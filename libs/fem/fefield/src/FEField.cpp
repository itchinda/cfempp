//************************************************************************
// --- CFEM++ Library
// ---
// --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.
// --- ALL RIGHTS RESERVED
// ---
// --- This software is protected by international copyright laws.
// --- Use, distribution, or modification of this software in any form,
// --- source or binary, for personal, academic, or non-commercial 
// --- purposes is permitted free of charge, provided that this 
// --- copyright notice and this permission notice appear in all 
// --- copies and supporting documentation.
// ---
// --- The authors and contributors provide this code "as is" without 
// --- any express or implied warranty. In no event shall they be 
// --- held liable for any damages arising from the use of this software.
//************************************************************************

#include "FEField.h"
#include "Mesh.h"
#include "DynamicVector.h"
#include "ReferenceElementFactory.h"

#include <stdexcept>
#include <iostream>

namespace cfem
{
    /**
     * @brief Constructor for FEField.
     * @param name Unique name of the field (e.g., "velocity", "temperature").
     * @param space Shared pointer to the FESpace (Lagrange, DG, etc.) providing the basis.
     */
    FEField::FEField(ConstructorToken, std::shared_ptr<FESpace> space, const std::string &name)
        : TrialFunction(space, name), FESpaceObserver(name),
         m_data(PETSC_DECIDE, space->getGlobalNumDofs())
    {
        m_data.setAll(0.0);
    }

    /**
     * @brief Resizes the underlying PETSc vector and clears local buffers.
     * Triggered by FESpace when the mesh or numbering changes.
     */
    void FEField::reallocate()
    {
        CfemInt newSize = m_space->getGlobalNumDofs();
        m_data = PetscVector(PETSC_DECIDE, newSize);
        m_data.setAll(0.0);
        m_localValues.clear();
        
        CFEM_INFO << "Field '" << m_name << "' reallocated to global size: " << newSize;
    }

    /**
     * @brief Synchronizes the global PETSc vector to a local buffer for fast access.
     * Uses a PETSc Scatter-to-All strategy to ensure every MPI process has a full copy 
     * of the field values for evaluations.
     */
    void FEField::updateLocalBuffer() const
    {
        CfemInt totalDofs = m_space->getGlobalNumDofs();
        
        // 1. Assurer la taille du buffer (allocation unique au début, réutilisée ensuite)
        if (m_localValues.size() != static_cast<size_t>(totalDofs))
        {
            m_localValues.assign(totalDofs, 0.0);
        }

        const Vec& petscVec = m_data.getRaw();
        
        // Obtenir le nombre de processus MPI
        int mpi_size;
        MPI_Comm_size(PETSC_COMM_WORLD, &mpi_size);

        if (mpi_size == 1) 
        {
            // Sur une seule machine, on lit directement le tableau interne de PETSc
            const PetscScalar *array;
            VecGetArrayRead(petscVec, &array);
            
            std::copy(array, array + totalDofs, m_localValues.begin());
            
            VecRestoreArrayRead(petscVec, &array);
        }
        else 
        {
            // pour ne pas être recréés à chaque fois, mais pour l'instant :
            VecScatter ctx;
            Vec localVec;

            VecCreateSeq(PETSC_COMM_SELF, totalDofs, &localVec);
            VecScatterCreateToAll(petscVec, &ctx, &localVec);

            VecScatterBegin(ctx, petscVec, localVec, INSERT_VALUES, SCATTER_FORWARD);
            VecScatterEnd(ctx, petscVec, localVec, INSERT_VALUES, SCATTER_FORWARD);

            const PetscScalar *array;
            VecGetArrayRead(localVec, &array);
            std::copy(array, array + totalDofs, m_localValues.begin());
            VecRestoreArrayRead(localVec, &array);

            VecDestroy(&localVec);
            VecScatterDestroy(&ctx);
        }
    }

    /**
     * @brief Retrieves the value at a specific mesh node for a given component.
     * Agnostic to how DOFs are stored (CG vs DG).
     */
    CfemReal FEField::getAtNode(CfemInt meshNodeId, CfemInt component) const
    {
        CfemInt globalIdx = m_space->getGlobalDofFromMeshNode(meshNodeId, component);
        
        if (globalIdx == -1 || globalIdx >= static_cast<CfemInt>(m_localValues.size()))
            CFEM_ERROR << "Mesh node ID out of bounds";
        
        return m_localValues[globalIdx];
    }

    /**
     * @brief Standard value evaluation at a local coordinate xi within a cell.
     */
    void FEField::evaluateWithoutShape(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        const Cell &cell = m_space->getMesh()->getCell(cellId);
        const auto &refEl = ReferenceElementFactory::getReferenceElement(cell.type, m_space->getOrderInt());

        refEl.evaluate(xi, ShapeEvaluationFlags::Value, ctx.scratchShape);
        this->evaluateWithShape(cellId, ctx.scratchShape, out);
    }

    /**
     * @brief Optimized evaluation using pre-computed shape functions.
     * Uses the FESpace abstraction to find global DOF indices.
     */
    void FEField::evaluateWithShape(CfemInt cellId, const ShapeInfo &shape, std::span<CfemReal> out) const
    {
        int nComp = this->getDimension();
        std::fill_n(out.data(), nComp, 0.0);
        
        const auto &shapeValues = shape.values;
        int nNodes = static_cast<int>(shapeValues.size());

        for (int i = 0; i < nNodes; ++i)
        {
            CfemReal s = shapeValues[i];
            for (int c = 0; c < nComp; ++c)
            {
                // Request the global index from the space abstraction
                CfemInt globalDofIdx = m_space->getGlobalDofIdx(cellId, i, c);
                
                if (globalDofIdx != -1) {
                    out[c] += s * m_localValues[globalDofIdx];
                }
            }
        }
    }

    void FEField::evaluateGradient(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        // ÉTAPE A : Est-ce qu'on est dans un intégrateur avec un mapping déjà prêt ?
        if (m_boundIndex != -1) {
            const MappingInfo* cachedMap = ctx.getMapping(m_boundIndex);
            if (cachedMap) {
                this->evaluateGradientWithCachedMapping(cellId, *cachedMap, out);
                return;
            }
        }

        // On appelle la version lourde qui recalcule la Jacobienne
        this->evaluateGradientWithoutCachedMapping(cellId, xi, out, ctx);
    }

    void FEField::evaluateGradientWithoutCachedMapping(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        const Cell &cell = m_space->getMesh()->getCell(cellId);
        const auto &refEl = ReferenceElementFactory::getReferenceElement(cell.type, m_space->getOrderInt());

        m_space->getMesh()->getCellVerticesCoords(cellId, ctx.scratchCoords);
        
        int expectedNodes = refEl.getNumNodes();
        if (ctx.scratchCoords.size() > (size_t)expectedNodes) {
            ctx.scratchCoords.resize(expectedNodes);
        }
        
        refEl.computeMapping(xi, ctx.scratchCoords, ShapeEvaluationFlags::Gradient, 
                                ctx.scratchShape, ctx.scratchMapping);
        // Évaluer le gradient physique
        this->evaluateGradientWithCachedMapping(cellId, ctx.scratchMapping, out);
    }

    /**
     * @brief Optimized gradient evaluation using pre-computed physical mapping.
     */
    void FEField::evaluateGradientWithCachedMapping(CfemInt cellId, const MappingInfo &mapping, std::span<CfemReal> out) const
    {
        int nComp = this->getDimension();
        
        std::fill_n(out.data(), nComp * 3, 0.0);

        const auto &grads = mapping.physicalShapeGradients;
        int nNodes = static_cast<int>(grads.size());

        for (int i = 0; i < nNodes; ++i)
        {
            const Vector3D &g = grads[i];
            for (int c = 0; c < nComp; ++c)
            {
                CfemInt globalDofIdx = m_space->getGlobalDofIdx(cellId, i, c);
                if (globalDofIdx != -1)
                {
                    CfemReal val = m_localValues[globalDofIdx];
                    out[c * 3 + 0] += val * g.x;
                    out[c * 3 + 1] += val * g.y;
                    out[c * 3 + 2] += val * g.z; 
                }
            }
        }
    }

    /**
     * @brief Implementation of the CfemExpression evaluation interface.
     */
    void FEField::evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        if (m_boundIndex >= 0) {
            const auto& data = ctx.getSpaceData(m_boundIndex);
            
            if (data.shape != nullptr) {
                this->evaluateWithShape(cellId, *(data.shape), out);
                return;
            }
        }

        this->evaluateWithoutShape(cellId, xi, out, ctx);
    }

    /**
     * @brief Projects a symbolic expression onto the FE Space (L2 Projection / Interpolation).
     * This handles the parallel MPI assembly and OpenMP node-visiting logic.
     */
    void FEField::project(const sym::CfemExpression &recipe)
    {
        CFEM_INFO << "projecting field into " << getName();

        
        recipe.prepare();

        const auto& mesh = m_space->getMesh();
        Vec& petscVec = m_data.getRaw();
        int dim = this->getDimension();

        const auto recipeFlags = recipe.requiredShapeFlags();

        // Determine the local ownership range for MPI
        CfemInt low, high;
        VecGetOwnershipRange(petscVec, &low, &high);
        CfemInt localSize = high - low;

        m_localValues.assign(localSize, 0.0);

        // Prepare for node visiting (CG: vertex-based, DG: cell-based)
        CfemInt numPositions = m_space->getNumUniquePositions();
        std::vector<uint8_t> visited(numPositions, 0);

        #pragma omp parallel
        {
            EvalContext evctx;
            std::vector<CfemReal> valBuffer(dim);
            std::span<CfemReal> valSpan(valBuffer);
            
            #pragma omp for schedule(guided)
            for (CfemInt cellId = 0; cellId < mesh->getNumCells(); ++cellId)
            {
                auto cell_type = mesh->getCell(cellId).type;
                const auto &refEl = ReferenceElementFactory::getReferenceElement(cell_type, m_space->getOrderInt());
                const auto &refElGeo = ReferenceElementFactory::getReferenceElement(cell_type, mesh->getOrderInt());

                mesh->getCellVerticesCoords(cellId, evctx.scratchCoords);

                for (int i = 0; i < (int)refEl.getNumNodes(); ++i)
                {
                    // Get the position ID for locking/visiting logic
                    CfemInt posId = m_space->getPositionId(cellId, i);
                    if (posId == -1) continue;

                    // Atomic check: Ensure each unique node is computed only once
                    if (visited[posId]) continue;
                    if (__atomic_test_and_set(&visited[posId], __ATOMIC_RELAXED)) continue;

                    // Evaluate the recipe at the physical location of the node
                    Vector3D xi = refEl.getNodeCoords(i);
                    // ShapeInfo nodeShape;
                    refElGeo.evaluate(xi, ShapeEvaluationFlags::Value | recipeFlags, evctx.scratchShape);
                    refElGeo.computeMapping(
                        xi, 
                        evctx.scratchCoords, 
                        ShapeEvaluationFlags::PhysicalPosition | recipeFlags, 
                        static_cast<const ShapeInfo&>(evctx.scratchShape), 
                        evctx.scratchMapping
                    );

                    evctx.physPos = evctx.scratchMapping.physicalPosition;
                    evctx.currentCellId = cellId;
                    evctx.currentXi = xi;
                    evctx.geometricMapping = &evctx.scratchMapping;

                    recipe.evaluate(cellId, xi, valSpan, evctx);

                    // Map calculated values to the global algebraic vector
                    for (int d = 0; d < dim; ++d)
                    {
                        CfemInt globalDofIdx = m_space->getGlobalDofIdx(cellId, i, d);

                        // MPI Check: Does this process own this specific DOF?
                        if (globalDofIdx >= low && globalDofIdx < high)
                        {
                            m_localValues[globalDofIdx - low] = valSpan[d];
                        }
                    }
                }
            }
        }

        // 3. Final PETSc Assembly
        PetscScalar *array;
        VecGetArray(petscVec, &array);
        std::copy(m_localValues.begin(), m_localValues.end(), array);
        VecRestoreArray(petscVec, &array);
        
        m_data.assemble();
        
        CFEM_INFO << "Projection of recipe onto field '" << m_name << "' completed.";
    }
}

namespace cfem::sym
{
    std::shared_ptr<FEField> Trial(std::shared_ptr<FESpace> V, const std::string& name)
    {
        sym::VariableRegistry::registerVariable(name);
        return FEField::create(name, V);
    }

    std::shared_ptr<const TestFunction> Test(std::shared_ptr<FESpace> V, const std::string& name) 
    {
        sym::VariableRegistry::registerVariable(name);
        
        return std::make_shared<TestFunction>(V, name);
    }
} // namespace cfem