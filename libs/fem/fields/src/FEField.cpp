#include "FEField.h"
#include "DynamicVector.h"
#include <stdexcept>
#include <iostream>

namespace cfem
{
    /**
     * @brief Constructor for FEField.
     * @param name Unique name of the field (e.g., "velocity", "temperature").
     * @param space Shared pointer to the FESpace (Lagrange, DG, etc.) providing the basis.
     */
    FEField::FEField(ConstructorToken, const std::string &name, std::shared_ptr<FESpace> space)
        : m_name(name),
          m_dimension(space->getDimension()),
          m_space(std::move(space)),
          m_data(PETSC_DECIDE, m_space->getGlobalNumDofs())
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
        
        // Ensure the local cache is sized correctly
        if (m_localValues.size() != static_cast<size_t>(totalDofs))
        {
            m_localValues.assign(totalDofs, 0.0);
        }

        Vec petscVec = m_data.getRaw();
        VecScatter ctx;
        Vec localVec;

        // Create a sequential vector on the local process to receive the full global data
        VecCreateSeq(PETSC_COMM_SELF, totalDofs, &localVec);
        VecScatterCreateToAll(petscVec, &ctx, &localVec);

        // Perform the communication
        VecScatterBegin(ctx, petscVec, localVec, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ctx, petscVec, localVec, INSERT_VALUES, SCATTER_FORWARD);

        // Extract values from the local PETSc vector to our std::vector buffer
        PetscScalar *array;
        VecGetArray(localVec, &array);
        std::copy(array, array + totalDofs, m_localValues.begin());
        VecRestoreArray(localVec, &array);

        // Clean up PETSc temporary objects
        VecDestroy(&localVec);
        VecScatterDestroy(&ctx);
    }

    /**
     * @brief Retrieves the value at a specific mesh node for a given component.
     * Agnostic to how DOFs are stored (CG vs DG).
     */
    CfemReal FEField::getAtNode(CfemInt meshNodeId, int component) const
    {
        CfemInt globalIdx = m_space->getGlobalDofFromMeshNode(meshNodeId, component);
        
        if (globalIdx == -1 || globalIdx >= static_cast<CfemInt>(m_localValues.size()))
            return 0.0;
        
        return m_localValues[globalIdx];
    }

    /**
     * @brief Standard value evaluation at a local coordinate xi within a cell.
     */
    void FEField::evaluateWithoutShape(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        const Cell &cell = m_space->getMesh()->getCell(cellId);
        CellType baseType = cell.type;

        if (m_space->getOrder() == Order::P1){
            baseType = Cell::getLinearEquivalent(baseType);
        }

        const auto &refEl = m_space->getRefElement(baseType);

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
        // 1. Déterminer le type de cellule correct pour l'espace (P1 ou P2)
        const Cell &cell = m_space->getMesh()->getCell(cellId);
        CellType baseType = cell.type;
        
        // Si l'espace est P1 mais le maillage est P2, on force le type linéaire (ex: TET10 -> TET4)
        if (m_space->getOrder() == Order::P1) {
            baseType = Cell::getLinearEquivalent(cell.type);
        }

        const auto &refEl = m_space->getRefElement(baseType);
        
        m_space->getMesh()->getCellVerticesCoords(cellId, ctx.scratchCoords);
        
        int expectedNodes = refEl.getNumNodes();
        if (ctx.scratchCoords.size() > (size_t)expectedNodes) {
            ctx.scratchCoords.resize(expectedNodes);
        }
        
        refEl.computeMapping(xi, ctx.scratchCoords, ShapeEvaluationFlags::Gradient, 
                                ctx.scratchShape, ctx.scratchMapping);
        // 4. Évaluer le gradient physique
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
        auto mesh = m_space->getMesh();
        Vec petscVec = m_data.getRaw();
        int dim = this->getDimension();

        // 1. Determine the local ownership range for MPI
        CfemInt low, high;
        VecGetOwnershipRange(petscVec, &low, &high);
        CfemInt localSize = high - low;

        m_localValues.assign(localSize, 0.0);

        // 2. Prepare for node visiting (CG: vertex-based, DG: cell-based)
        CfemInt numPositions = m_space->getNumUniquePositions();
        std::vector<uint8_t> visited(numPositions, 0);

        #pragma omp parallel
        {
            EvalContext evctx;
            std::vector<CfemReal> valBuffer(dim);
            std::span<CfemReal> valSpan(valBuffer);
            
            DynamicVector<Vector3D> physCoords;
            MappingInfo map;

            #pragma omp for schedule(guided)
            for (CfemInt cellId = 0; cellId < mesh->getNumCells(); ++cellId)
            {
                const auto &cell = mesh->getCell(cellId);
                const auto &refEl = m_space->getRefElement(cell.type);
                mesh->getCellVerticesCoords(cellId, physCoords);

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
                    ShapeInfo nodeShape;
                    refEl.evaluate(xi, ShapeEvaluationFlags::Value, nodeShape);
                    refEl.computeMapping(xi, physCoords, ShapeEvaluationFlags::None, nodeShape, map);

                    evctx.physPos = map.physicalPosition;
                    evctx.currentCellId = cellId;
                    evctx.currentXi = xi;

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

} // namespace cfem