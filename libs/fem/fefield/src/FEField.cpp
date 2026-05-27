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
        
        // Assurer la taille du buffer (allocation unique au début, réutilisée ensuite)
        if (m_localValues.size() != static_cast<size_t>(totalDofs))
        {
            m_localValues.assign(totalDofs, 0.0);
        }

        const Vec& petscVec = m_data.getRaw();
        
        // Obtenir le nombre de processus MPI
        CfemInt mpi_size;
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

    CfemReal FEField::getScalarFieldAtNode(CfemInt meshNodeId) const
    {
        CFEM_ASSERT(m_space->getFieldType() == FieldType::Scalar, "");

        CfemReal scalar = getAtNode(meshNodeId, 0);

        return scalar;
    }

    Vector2D FEField::get2DVectorFieldAtNode    (CfemInt meshNodeId) const
    {
        CFEM_ASSERT(m_space->getFieldType() == FieldType::Vector2D, "" );

        Vector2D vector;
        vector.x = getAtNode(meshNodeId, 0);
        vector.y = getAtNode(meshNodeId, 1);

        return vector;

    }

    Tensor2D FEField::get2DTensorFieldAtNode (CfemInt meshNodeId) const
    {
        CFEM_ASSERT(m_space->getFieldType() == FieldType::Tensor2D, "" );

        Tensor2D tensor;
        tensor.data[Tensor2D::idx_xx] = getAtNode(meshNodeId, Tensor2D::idx_xx);
        tensor.data[Tensor2D::idx_xy] = getAtNode(meshNodeId, Tensor2D::idx_xy);
        tensor.data[Tensor2D::idx_yx] = getAtNode(meshNodeId, Tensor2D::idx_yx);
        tensor.data[Tensor2D::idx_yy] = getAtNode(meshNodeId, Tensor2D::idx_yy);

        return tensor;

    }

    SymTensor2D FEField::get2DSymTensorFieldAtNode (CfemInt meshNodeId) const
    {
        CFEM_ASSERT(m_space->getFieldType() == FieldType::SymTensor2D, "" );

        SymTensor2D tensor;
        tensor.data[SymTensor2D::idx_xx] = getAtNode(meshNodeId, SymTensor2D::idx_xx);
        tensor.data[SymTensor2D::idx_xy] = getAtNode(meshNodeId, SymTensor2D::idx_xy);
        tensor.data[SymTensor2D::idx_yy] = getAtNode(meshNodeId, SymTensor2D::idx_yy);

        return tensor;
    }


    Vector3D FEField::get3DVectorFieldAtNode(CfemInt meshNodeId) const
    {
        CFEM_ASSERT(m_space->getFieldType() == FieldType::Vector3D, "" );

        Vector3D vector;
        vector.x = getAtNode(meshNodeId, 0);
        vector.y = getAtNode(meshNodeId, 1);
        vector.z = getAtNode(meshNodeId, 2);

        return vector;
    }

    Tensor3D FEField::get3DTensorFieldAtNode(CfemInt meshNodeId) const
    {
        CFEM_ASSERT(m_space->getFieldType() == FieldType::Tensor3D, "" );

        Tensor3D tensor;
        tensor.data[Tensor3D::idx_xx] = getAtNode(meshNodeId, Tensor3D::idx_xx);
        tensor.data[Tensor3D::idx_xy] = getAtNode(meshNodeId, Tensor3D::idx_xy);
        tensor.data[Tensor3D::idx_xz] = getAtNode(meshNodeId, Tensor3D::idx_xz);
        tensor.data[Tensor3D::idx_yx] = getAtNode(meshNodeId, Tensor3D::idx_yx);
        tensor.data[Tensor3D::idx_yy] = getAtNode(meshNodeId, Tensor3D::idx_yy);
        tensor.data[Tensor3D::idx_yz] = getAtNode(meshNodeId, Tensor3D::idx_yz);
        tensor.data[Tensor3D::idx_zx] = getAtNode(meshNodeId, Tensor3D::idx_zx);
        tensor.data[Tensor3D::idx_zy] = getAtNode(meshNodeId, Tensor3D::idx_zy);
        tensor.data[Tensor3D::idx_zz] = getAtNode(meshNodeId, Tensor3D::idx_zz);

        return tensor;
    }

    SymTensor3D FEField::get3DSymTensorFieldAtNode(CfemInt meshNodeId) const
    {
        CFEM_ASSERT(m_space->getFieldType() == FieldType::SymTensor3D, "" );

        SymTensor3D tensor;
        tensor.data[SymTensor3D::idx_xx] = getAtNode(meshNodeId, SymTensor3D::idx_xx); 
        tensor.data[SymTensor3D::idx_yy] = getAtNode(meshNodeId, SymTensor3D::idx_yy); 
        tensor.data[SymTensor3D::idx_zz] = getAtNode(meshNodeId, SymTensor3D::idx_zz); 
        tensor.data[SymTensor3D::idx_yz] = getAtNode(meshNodeId, SymTensor3D::idx_yz); 
        tensor.data[SymTensor3D::idx_xz] = getAtNode(meshNodeId, SymTensor3D::idx_xz); 
        tensor.data[SymTensor3D::idx_xy] = getAtNode(meshNodeId, SymTensor3D::idx_xy); 

        return tensor;
    }

    SkewSymTensor3D FEField::get3DSkewSymTensorFieldAtNode(CfemInt meshNodeId) const
    {
        CFEM_ASSERT(m_space->getFieldType() == FieldType::SkewSymTensor3D, "");

        SkewSymTensor3D tensor;
        tensor.data[SkewSymTensor3D::idx_yz] = getAtNode(meshNodeId, SkewSymTensor3D::idx_yz); 
        tensor.data[SkewSymTensor3D::idx_xz] = getAtNode(meshNodeId, SkewSymTensor3D::idx_xz); 
        tensor.data[SkewSymTensor3D::idx_xy] = getAtNode(meshNodeId, SkewSymTensor3D::idx_xy);

        return tensor;

    }

    /**
     * @brief Extracts cell degrees of freedom from the synchronized global local buffer into a contiguous matrix.
     * @param cellId Local index of the computational cell.
     * @param localDofs Reference to the scratch matrix [nNodes x fieldDimension] allocated in the memory arena.
     * @param ctx Unified high-performance evaluation context.
     */
    void FEField::extractCellDofs(CfemInt cellId, DynamicMatrix<CfemReal>& localDofs, EvalContext& ctx) const 
    {
        const CfemInt nNodes = m_space->getReferenceElement(cellId).getNumNodes();
        const CfemInt fDim   = this->getNumComponents(); // m_dim (e.g., 1 for scalar fields, 3 for 3D vectors)
        
        // Resize the scratch matrix inside the memory arena (zero allocation overhead if capacity suffices)
        localDofs.resize(nNodes, fDim);

        // Populate the dense local layout using direct global index lookup from FESpace
        for (CfemInt i = 0; i < nNodes; ++i) {
            for (CfemInt d = 0; d < fDim; ++d) {
                const CfemInt globalDofIdx = m_space->getGlobalDofIdx(cellId, i, d);
                
                // Safety check: Ensure the global index maps inside our synchronized cache buffer bounds
                CFEM_ASSERT(globalDofIdx >= 0 && globalDofIdx < static_cast<CfemInt>(m_localValues.size()),
                    "FEField::extractCellDofs - Resolved global DOF index is out of bounds of the synchronized m_localValues cache.");
                    
                localDofs[i][d] = m_localValues[globalDofIdx];
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
            
            if (data.shapeSingle != nullptr) {
                this->evaluateWithShape(cellId, *(data.shapeSingle), out);
                return;
            }
        }

        this->evaluateWithoutShape(cellId, xi, out, ctx);
    }


    /**
     * @brief Evaluates the standard field value at a single local coordinate xi without pre-computed shapes.
     * @param cellId Local index of the computational cell.
     * @param xi Local/Reference coordinates inside the reference element space.
     * @param out Destination memory span packing the physical field components [m_dim].
     * @param ctx Unified high-performance evaluation context.
     */
    void FEField::evaluateWithoutShape(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        const Cell &cell = m_space->getMesh()->getCell(cellId);
        const auto &refEl = ReferenceElementFactory::getReferenceElement(cell.type, m_space->getOrderInt());

        ctx.singleWS.shape.setup(refEl.getNumNodes(), false);
        ctx.singleWS.shape.validFlags |= ShapeEvaluationFlags::Value;
        refEl.evaluateShapesValues(xi, ctx.singleWS.shape);
        
        this->evaluateWithShape(cellId, ctx.singleWS.shape, out);
    }

    /**
     * @brief Interpolates field values at a single integration point using pre-allocated ShapeInfo.
     * @param cellId Local index of the computational cell.
     * @param shape Single-point reference shape container holding base function values.
     * @param out Destination memory span packing the physical field components [m_dim].
     */
    void FEField::evaluateWithShape(CfemInt cellId, const ShapeInfo &shape, std::span<CfemReal> out) const
    {
        const CfemInt nComp = this->getNumComponents();
        std::fill_n(out.data(), nComp, 0.0);
        
        const auto &shapeValues = shape.values();
        const CfemInt nNodes = shape.numNodes;

        // Contract nodal degrees of freedom linearly against element values
        for (CfemInt i = 0; i < nNodes; ++i)
        {
            const CfemReal s = shapeValues[i];
            for (CfemInt c = 0; c < nComp; ++c)
            {
                const CfemInt globalDofIdx = m_space->getGlobalDofIdx(cellId, i, c);
                
                // Bypass dead or unmapped DOFs safely (e.g., in specialized embedded spaces)
                if (globalDofIdx != -1) 
                {
                    CFEM_ASSERT(globalDofIdx >= 0 && globalDofIdx < static_cast<CfemInt>(m_localValues.size()),
                        "FEField::evaluateWithShape - Resolved global DOF index out of local cache bounds.");
                        
                    out[c] += s * m_localValues[globalDofIdx];
                }
            }
        }
    }

    /**
     * @brief Gateway method to evaluate the spatial physical gradient tensor at a single local coordinate.
     * @param cellId Local index of the computational cell.
     * @param xi Local/Reference coordinates inside the reference element space.
     * @param out Destination memory span packing the flattened gradient tensor [m_dim x spaceDim].
     * @param ctx Unified high-performance evaluation context.
     */
    void FEField::evaluateGradient(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        if (m_boundIndex != -1) 
        {
            const MappingInfo* cachedMap = ctx.getMapping(m_boundIndex);
            if (cachedMap) 
            {
                this->evaluateGradientWithCachedMapping(cellId, *cachedMap, out);
                return;
            }
        }

        this->evaluateGradientWithoutCachedMapping(cellId, xi, out, ctx);
    }

    /**
     * @brief Evaluates the physical spatial gradient tensor at an un-cached arbitrary local coordinate.
     * @param cellId Local index of the computational cell.
     * @param xi Local/Reference coordinates inside the reference element space.
     * @param out Destination memory span packing the flattened physical gradient tensor.
     * @param ctx Unified high-performance evaluation context.
     */
    void FEField::evaluateGradientWithoutCachedMapping(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        const auto &refEl = m_space->getReferenceElement(cellId);
        
        // Extract contiguous cell vertex coordinates into context scratchpad
        m_space->getMesh()->getCellVerticesCoords(cellId, ctx.scratchCoords);
        
        const CfemInt expectedNodes = refEl.getNumNodes();
        if (ctx.scratchCoords.size() > static_cast<size_t>(expectedNodes)) {
            ctx.scratchCoords.resize(expectedNodes);
        }
                
        ctx.singleWS.shape.shapeMustBeUpdated = false;
        
        refEl.computeMapping(
            xi, 
            ctx.scratchCoords,
            ShapeEvaluationFlags::FirstDerivatives,
            MappingEvaluationFlags::PhysicalFirstDerivatives,  // all internal dependency requirements (Jacobian, Determinant, InverseJacobian).
            ctx.singleWS.shape, 
            ctx.singleWS.mapping
        );
        
        // Forward the computed single metric point to our corrected tensor projection method
        this->evaluateGradientWithCachedMapping(cellId, ctx.singleWS.mapping, out);
    }

   /**
     * @brief Evaluates the physical gradient tensor at a single integration point using a cached MappingInfo.
     * @param cellId Local index of the computational cell.
     * @param mapping Single-point metric context containing precomputed physical shape gradients.
     * @param out Destination memory span packing the flattened gradient tensor.
     * Layout: [component * spatialDimension + spatialDirection]
     */
    void FEField::evaluateGradientWithCachedMapping(CfemInt cellId, const MappingInfo &mapping, std::span<CfemReal> out) const
    {
        const CfemInt nComp = this->getNumComponents();
        const CfemInt spatialDimension = m_space->getMesh()->getSpatialDimension();

        // Reset output buffer tracking the exact required tensor memory footprint
        std::fill_n(out.data(), nComp * spatialDimension, 0.0);

        const auto &gradsX = mapping.dN_dx();
        const auto &gradsY = mapping.dN_dy();
        const auto &gradsZ = mapping.dN_dz();
        const CfemInt nNodes = mapping.numNodes;  //  static_cast<CfemInt>(gradsX.size());

        // Loop over element basis functions (nodes)
        for (CfemInt i = 0; i < nNodes; ++i)
        {
            const CfemReal gx = gradsX[i];
            const CfemReal gy = (spatialDimension >= 2) ? gradsY[i] : 0.0;
            const CfemReal gz = (spatialDimension == 3) ? gradsZ[i] : 0.0;

            // Loop over physical field components
            for (CfemInt c = 0; c < nComp; ++c)
            {
                const CfemInt globalDofIdx = m_space->getGlobalDofIdx(cellId, i, c);
                
                CFEM_ASSERT(globalDofIdx >= 0 && globalDofIdx < static_cast<CfemInt>(m_localValues.size()),
                    "FEField::evaluateGradientWithCachedMapping - Resolved global DOF index is out of bounds.");

                const CfemReal val = m_localValues[globalDofIdx];
                const CfemInt rowOffset = c * spatialDimension;

                // Dynamically project derivatives based on the true dimensionality of the mesh
                out[rowOffset + 0] += val * gx;
                if (spatialDimension > 1) {
                    out[rowOffset + 1] += val * gy;
                    if (spatialDimension == 3) {
                        out[rowOffset + 2] += val * gz;
                    }
                }
            }
        }
    }


    // ========================================================================
    // BATCH EVALUATION (VALEURS)
    // ========================================================================

    void FEField::evaluateBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const
    {
        if (m_boundIndex >= 0) {
            const auto* shapeBatch = ctx.getShapeBatch(m_boundIndex);
            
            if (shapeBatch != nullptr) {
                this->evaluateBatchWithShape(cellId, nQp, *shapeBatch, outValues, ctx);
                return;
            }
        }

        this->evaluateBatchWithoutShape(cellId, nQp, outValues, ctx);
    }

    void FEField::evaluateBatchWithShape(CfemInt cellId, CfemInt nQp, const ShapeInfoBatch& shapeBatch, std::span<CfemReal> outValues, EvalContext& ctx) const
    {
        const CfemInt nComp = this->getNumComponents();
        const CfemInt nNodes = shapeBatch.numNodes;

        // On s'assure que le buffer de sortie a la bonne taille (SoA)
        CFEM_ASSERT(outValues.size() >= nComp * nQp, "FEField::evaluateBatchWithShape - Output buffer too small");
        std::fill_n(outValues.data(), nComp * nQp, 0.0);

        auto bookmark = ctx.getBookmark();
        this->extractCellDofs(cellId, ctx.scratchCellDofs, ctx); 

        for (CfemInt c = 0; c < nComp; ++c) {
            for (CfemInt q = 0; q < nQp; ++q) {
                const CfemReal* N = shapeBatch.values(q);
                
                CfemReal sum = 0.0;

                for (CfemInt i = 0; i < nNodes; ++i) {
                    sum += N[i] * ctx.scratchCellDofs[i][c];
                }
                outValues[c * nQp + q] = sum;
            }
        }

        ctx.releaseToBookmark(bookmark);
    }

    void FEField::evaluateBatchWithoutShape(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const
    {
        // Récupération sécurisée des coordonnées xiBatch depuis le contexte (doit être configuré par l'intégrateur)
        const auto* xiBatch = ctx.batchWS.xi;

        CFEM_ASSERT(xiBatch != nullptr && xiBatch->size() >= static_cast<size_t>(nQp), "FEField::evaluateBatchWithoutShape requires xiBatch to be populated in EvalContext.");

        const Cell &cell = m_space->getMesh()->getCell(cellId);
        const auto &refEl = ReferenceElementFactory::getReferenceElement(cell.type, m_space->getOrderInt());

        // Configuration temporaire du ShapeInfoBatch dans le contexte
        ctx.batchWS.shape->setup(nQp, refEl.getNumNodes(), false);
        ctx.batchWS.shape->validFlags |= ShapeEvaluationFlags::Value;
        
        refEl.evaluateShapesValuesBatch(*xiBatch, *(ctx.batchWS.shape));
        
        this->evaluateBatchWithShape(cellId, nQp, *(ctx.batchWS.shape), outValues, ctx);
    }


    // ========================================================================
    // BATCH EVALUATION (GRADIENTS)
    // ========================================================================

    void FEField::evaluateGradientBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const
    {
        if (m_boundIndex >= 0) {
            const auto* mappingBatch = ctx.getMappingBatch(m_boundIndex);
            
            if (mappingBatch != nullptr) {
                this->evaluateGradientBatchWithCachedMapping(cellId, nQp, *mappingBatch, outValues, ctx);
                return;
            }
            CFEM_INFO << "WARNING: Bound index present but no MappingBatch cached for FEField::evaluateGradientBatch";
        }

        this->evaluateGradientBatchWithoutCachedMapping(cellId, nQp, outValues, ctx);
    }

    void FEField::evaluateGradientBatchWithCachedMapping(CfemInt cellId, CfemInt nQp, const MappingInfoBatch& mappingBatch, std::span<CfemReal> outValues, EvalContext& ctx) const
    {
        const CfemInt nComp = this->getNumComponents();
        const CfemInt spatialDim = m_space->getMesh()->getSpatialDimension();
        const CfemInt nNodes = mappingBatch.numNodes;

        // Le format de sortie pour un tenseur (ou gradient de vecteur) en SoA:
        // Taille : (nComp * spatialDim) * nQp
        const CfemInt expectedSize = nComp * spatialDim * nQp;
        CFEM_ASSERT(outValues.size() >= static_cast<size_t>(expectedSize), "FEField::evaluateGradientBatchWithCachedMapping - Output buffer too small");
        std::fill_n(outValues.data(), expectedSize, 0.0);

        auto bookmark = ctx.getBookmark();
        this->extractCellDofs(cellId, ctx.scratchCellDofs, ctx); 

        // Boucle sur les composantes du champ (ex: u_x, u_y, u_z)
        for (CfemInt c = 0; c < nComp; ++c) 
        {
            const CfemInt baseRow = c * spatialDim; // Offset pour la composante c

            // Boucle interne sur les QP : parfait pour la prédiction de branchement et le chargement cache !
            for (CfemInt q = 0; q < nQp; ++q) 
            {
                const CfemReal* gx = mappingBatch.dN_dx(q);
                const CfemReal* gy = (spatialDim >= 2) ? mappingBatch.dN_dy(q) : nullptr;
                const CfemReal* gz = (spatialDim == 3) ? mappingBatch.dN_dz(q) : nullptr;

                CfemReal sumX = 0.0, sumY = 0.0, sumZ = 0.0;

                // Contraction locale (nœuds)
                for (CfemInt i = 0; i < nNodes; ++i) 
                {
                    const CfemReal dofVal = ctx.scratchCellDofs[i][c];
                    sumX += dofVal * gx[i];
                    if (spatialDim >= 2) sumY += dofVal * gy[i];
                    if (spatialDim == 3) sumZ += dofVal * gz[i];
                }

                // Écriture respectant le layout SoA : out[(composante_tensorielle) * nQp + q]
                outValues[(baseRow + 0) * nQp + q] = sumX;
                if (spatialDim >= 2) outValues[(baseRow + 1) * nQp + q] = sumY;
                if (spatialDim == 3) outValues[(baseRow + 2) * nQp + q] = sumZ;
            }
        }

        ctx.releaseToBookmark(bookmark);
    }

    void FEField::evaluateGradientBatchWithoutCachedMapping(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const
    {
        // Pareil que pour la valeur batch : on a besoin de xiBatch (via EvalContext)
        const auto& xiBatch = ctx.batchWS.xi;
        CFEM_ASSERT(xiBatch != nullptr && xiBatch->size() >= static_cast<size_t>(nQp), "FEField::evaluateGradientBatchWithoutCachedMapping requires xiBatch");

        const auto &refEl = m_space->getReferenceElement(cellId);
        
        m_space->getMesh()->getCellVerticesCoords(cellId, ctx.scratchCoords);
        const CfemInt expectedNodes = refEl.getNumNodes();
        if (ctx.scratchCoords.size() > static_cast<size_t>(expectedNodes)) {
            ctx.scratchCoords.resize(expectedNodes);
        }

        // On demande à l'élément de référence de construire la totale pour nQp points
        ctx.batchWS.shape->shapeMustBeUpdated = false;

        const bool isAffine = (m_space->getMesh()->getOrderInt() == 1) || m_space->getMesh()->isCellAffine(cellId);
        
        refEl.computeMappingBatch(
            *xiBatch, 
            ctx.scratchCoords,
            ShapeEvaluationFlags::FirstDerivatives,
            MappingEvaluationFlags::PhysicalFirstDerivatives,
            *(ctx.batchWS.shape), 
            *(ctx.batchWS.mapping)
        );
        
        this->evaluateGradientBatchWithCachedMapping(cellId, nQp, *(ctx.batchWS.mapping), outValues, ctx);
    }

    /**
     * @brief Projects a symbolic mathematical expression (recipe) onto the finite element field via nodal interpolation.
     * @param recipe The symbolic expression containing the mathematical definition to evaluate.
     */
    void FEField::interpolate(const sym::CfemExpression &recipe)
    {
        CFEM_INFO << "Interpolating symbolic expression onto field '" << getName() << "'.";

        // Initialize or bind any symbolic parameters or sub-expressions
        recipe.prepare();

        const auto& mesh = m_space->getMesh();
        Vec& petscVec = m_data.getRaw();
        const CfemInt dim = this->getNumComponents();

        // Query mathematical dependencies to determine which reference derivatives are needed
        const auto recipeFlags = recipe.requiredShapeFlags();
        const bool needsGeoNormals = recipe.dependsOnGeometricNormal();
        // Align evaluation flags to match the exact mathematical types required by the system
        const ShapeEvaluationFlags shapeFlags = ShapeEvaluationFlags::Value | recipeFlags |
                                                        (needsGeoNormals? ShapeEvaluationFlags::FirstDerivatives : ShapeEvaluationFlags::None);
        const MappingEvaluationFlags mappingFlags =  MappingEvaluationFlags::PhysicalPosition | transformShapeToMappingFlags(recipeFlags) |
                                                    (needsGeoNormals? MappingEvaluationFlags::Normal : MappingEvaluationFlags::None);

        // Determine the local process partition range owned within the distributed global PETSc system
        CfemInt low = 0;
        CfemInt high = 0;
        VecGetOwnershipRange(petscVec, &low, &high);
        const CfemInt localSize = high - low;

        // Reset and size the local write buffer matching this process's MPI ownership slice
        m_localValues.assign(localSize, 0.0);

        // Track visited nodes to guarantee a single evaluation pass for shared entities (CG/DG neutral)
        const CfemInt numPositions = m_space->getNumUniquePositions();
        std::vector<uint8_t> visited(numPositions, 0);

        #pragma omp parallel
        {
            // Allocation of thread-local unified context. Core structures (Arena, cache) remain local.
            EvalContext evctx;
            
            std::vector<CfemReal> valBuffer(dim);
            std::span<CfemReal> valSpan(valBuffer);
            
            #pragma omp for schedule(guided)
            for (CfemInt cellId = 0; cellId < mesh->getNumCells(); ++cellId)
            {
                const auto cell_type = mesh->getCell(cellId).type;
                const auto &refEl = ReferenceElementFactory::getReferenceElement(cell_type, m_space->getOrderInt());
                const auto &refElGeo = ReferenceElementFactory::getReferenceElement(cell_type, mesh->getOrderInt());

                // Collect node physical coordinates of the current cell into the context scratchpad
                mesh->getCellVerticesCoords(cellId, evctx.scratchCoords);

                const CfemInt nNodes = static_cast<CfemInt>(refEl.getNumNodes());
                for (CfemInt i = 0; i < nNodes; ++i)
                {
                    // Obtain unique global spatial position mapping ID
                    const CfemInt posId = m_space->getPositionId(cellId, i);
                    if (posId == -1) continue;

                    // Atomic test-and-set fence to prevent concurrent evaluations of identical shared nodes
                    if (visited[posId]) continue;
                    if (__atomic_test_and_set(&visited[posId], __ATOMIC_RELAXED)) continue;

                    // Map target reference element node coordinates
                    const Vector3D xi = refEl.getNodeCoords(i);

                    // Compute metric mapping properties at the specific isolated node location
                    refElGeo.computeMapping(
                        xi, 
                        evctx.scratchCoords, 
                        shapeFlags,
                        mappingFlags, 
                        evctx.singleWS.shape, 
                        evctx.singleWS.mapping
                    );

                    // Populate the single-point workspace with current state configurations
                    evctx.singleWS.physPos = evctx.singleWS.mapping.physicalPosition;
                    evctx.singleWS.xi      = xi;
                    evctx.geometricMapping = &evctx.singleWS.mapping;
                    evctx.currentCellId     = cellId;

                    // Bind the active mapping pointer to our newly established type-safe single slots
                    // evctx.setSpaceData(m_boundIndex, &evctx.singleWS.shape, &evctx.singleWS.mapping);

                    // Execute symbolic mathematical projection directly into thread-local buffer
                    recipe.evaluate(cellId, xi, valSpan, evctx);

                    // Map calculated values to the local algebraic partition slice
                    for (CfemInt d = 0; d < dim; ++d)
                    {
                        const CfemInt globalDofIdx = m_space->getGlobalDofIdx(cellId, i, d);

                        // MPI ownership range check: store only if it maps onto this local process memory segment
                        if (globalDofIdx >= low && globalDofIdx < high)
                        {
                            const size_t localProcessIdx = static_cast<size_t>(globalDofIdx - low);
                            m_localValues[localProcessIdx] = valSpan[d];
                        }
                    }
                }
            }
        }

        // 3. Complete distributed parallel synchronization into global PETSc architecture
        PetscScalar *array = nullptr;
        VecGetArray(petscVec, &array);
        std::copy(m_localValues.begin(), m_localValues.end(), array);
        VecRestoreArray(petscVec, &array);
        
        m_data.assemble();
        
        CFEM_INFO << "Nodal interpolation of recipe onto field '" << m_name << "' completed successfully.";
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