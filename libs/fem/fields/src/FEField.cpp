#include "FEField.h"
#include "ReferenceElementFactory.h"
#include "DynamicVector.h"

#include <stdexcept>
#include <iostream>
#include <malloc.h>

namespace cfem
{

    FEField::FEField(ConstructorToken, const std::string &name, std::shared_ptr<LagrangeSpace> space)
        : m_name(name),
          m_dimension(space->getDimension()),
          m_space(std::move(space)),
          m_data(PETSC_DECIDE, m_space->getGlobalNumDofs())
    {
        m_data.setAll(0.0);
    }

    void FEField::reallocate()
    {
        CfemInt newSize = m_space->getGlobalNumDofs();
        m_data = PetscVector(PETSC_DECIDE, newSize);
        m_data.setAll(0.0);
        m_localValues.clear();
        std::cout << "[FEField] Field '" << m_name << "' reallocated to size " << newSize << std::endl;
    }

    void FEField::updateLocalBuffer() const
    {
        CfemInt totalDofs = m_space->getGlobalNumDofs();
        if (m_localValues.size() != static_cast<size_t>(totalDofs))
        {
            m_localValues.assign(totalDofs, 0.0);
        }

        VecScatter ctx;
        Vec localVec;
        VecCreateSeq(PETSC_COMM_SELF, totalDofs, &localVec);
        VecScatterCreateToAll(m_data.getRaw(), &ctx, nullptr);

        VecScatterBegin(ctx, m_data.getRaw(), localVec, INSERT_VALUES, SCATTER_FORWARD);
        VecScatterEnd(ctx, m_data.getRaw(), localVec, INSERT_VALUES, SCATTER_FORWARD);

        PetscScalar *array;
        VecGetArray(localVec, &array);
        for (CfemInt i = 0; i < totalDofs; ++i)
        {
            m_localValues[i] = static_cast<CfemReal>(array[i]);
        }
        VecRestoreArray(localVec, &array);

        VecScatterDestroy(&ctx);
        VecDestroy(&localVec);
    }

    // void FEField::updateLocalBuffer() const {
    //     CfemInt totalDofs = m_space->getGlobalNumDofs();
        
    //     // 1. NE RIEN FAIRE si déjà à jour
    //     if (!m_isDirty && m_localValues.size() == (size_t)totalDofs) return;

    //     if (m_localValues.size() != (size_t)totalDofs) {
    //         m_localValues.assign(totalDofs, 0.0);
    //     }

    //     // 2. RÉUTILISER le contexte de communication (Indispensable pour la vitesse)
    //     if (!m_scatterCtx) {
    //         // Option A: Si tu es sur 1 seul process, utilise VecGetArray (Instant)
    //         // Option B: Si multi-process, garde le CreateToAll MAIS garde-le en membre de classe
    //         VecScatterCreateToAll(m_data.getRaw(), &m_scatterCtx, nullptr);
    //     }

    //     Vec tempLocalVec;
    //     VecCreateSeqWithArray(PETSC_COMM_SELF, 1, totalDofs, m_localValues.data(), &tempLocalVec);

    //     // 3. LE TRANSFERT (Maintenant ultra-rapide car le mapping est déjà prêt)
    //     VecScatterBegin(m_scatterCtx, m_data.getRaw(), tempLocalVec, INSERT_VALUES, SCATTER_FORWARD);
    //     VecScatterEnd(m_scatterCtx, m_data.getRaw(), tempLocalVec, INSERT_VALUES, SCATTER_FORWARD);

    //     VecDestroy(&tempLocalVec);
    //     m_isDirty = false;
    // }

    CfemReal FEField::getAtNode(CfemInt meshNodeId, int component) const
    {
        CfemInt nodeIdx = m_space->getLocalNodeIdx(meshNodeId);
        if (nodeIdx == -1)
            return 0.0;

        CfemInt numNodes = m_space->getNumActiveNodes();
        CfemInt globalIdx = nodeIdx + (component * numNodes);

        if (globalIdx >= static_cast<CfemInt>(m_localValues.size()))
            return 0.0;
        return m_localValues[globalIdx];
    }

    CfemReal FEField::getAt(CfemInt globalDofIdx, int component) const
    {
        if (globalDofIdx < 0)
            return 0.0;
        int numComponents = this->getDimension();
        return m_localValues[globalDofIdx * numComponents + component];
    }

    void FEField::evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        const ShapeInfo *shapeToUse = nullptr;
        if (m_space->getOrder() == Order::P1 && ctx.shapeP1)
        {
            shapeToUse = ctx.shapeP1;
        }
        else if (ctx.shapeP2)
        {
            shapeToUse = ctx.shapeP2;
        }

        ShapeInfo localShape;
        if (!shapeToUse)
        {
            CellType evalType = m_space->getMesh()->getCell(cellId).type;
            if (m_space->getOrder() == Order::P1)
            {
                evalType = Cell::getLinearEquivalent(evalType);
            }
            m_space->getRefElement(evalType).evaluate(xi, ShapeEvaluationFlags::Value, localShape);
            shapeToUse = &localShape;
        }

        if (shapeToUse)
        {
            this->evaluateWithShape(cellId, *shapeToUse, out);
        }
        else
        {
            this->evaluateAt(cellId, xi, out);
        }
    }

    void FEField::evaluateWithShape(CfemInt cellId, const ShapeInfo &shape, std::span<CfemReal> out) const
    {
        int nComp = this->getDimension();
        std::fill_n(out.data(), nComp, 0.0);

        const auto &cell_ids = m_space->getMesh()->getCellConnectivity(cellId);
        const auto &shapeValues = shape.values;
        CfemInt numNodes = m_space->getNumActiveNodes();

        for (size_t i = 0; i < shapeValues.size(); ++i)
        {
            CfemInt dofIdx = m_space->getLocalNodeIdx(cell_ids[i]);
            if (dofIdx != -1)
            {
                CfemReal s = shapeValues[i];
                for (int c = 0; c < nComp; ++c)
                {
                    out[c] += s * m_localValues[dofIdx + (c * numNodes)];
                }
            }
        }
    }

    void FEField::evaluateAt(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out) const
    {
        int nComp = this->getDimension();
        std::fill_n(out.data(), nComp, 0.0);

        const Cell &cell = m_space->getMesh()->getCell(cellId);
        CellType baseType = cell.type;

        // Logic for P1 enforcement (Fully Restored)
        if (m_space->getOrder() == Order::P1)
        {
            baseType = Cell::getLinearEquivalent(baseType);
        }

        const auto &refEl = m_space->getRefElement(baseType);
        ShapeInfo shape;
        refEl.evaluate(xi, ShapeEvaluationFlags::Value, shape);

        const auto &shapeValues = shape.values;
        const auto &cell_ids = m_space->getMesh()->getCellConnectivity(cell.id);

        CfemInt numNodes = m_space->getNumActiveNodes();
        for (size_t i = 0; i < shapeValues.size(); ++i)
        {
            CfemInt dofIdx = m_space->getLocalNodeIdx(cell_ids[i]);
            if (dofIdx != -1)
            {
                for (int c = 0; c < nComp; ++c)
                {
                    out[c] += shapeValues[i] * m_localValues[dofIdx + (c * numNodes)];
                }
            }
        }
    }

    void FEField::evaluateGradientAt(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        int nComp = this->getDimension(); 
        int outDim = nComp * 3;
        std::fill_n(out.data(), outDim, 0.0);

        const Cell &cell = m_space->getMesh()->getCell(cellId);
        CellType baseType = (m_space->getOrder() == Order::P1) ? 
                            Cell::getLinearEquivalent(cell.type) : cell.type;

        const auto &refEl = m_space->getRefElement(baseType);
        const auto &cell_ids = m_space->getMesh()->getCellConnectivity(cell.id);

        // Use Context scratchpads - NO HEAP ALLOC HERE
        m_space->getMesh()->getCellVerticesCoords(cellId, ctx.scratchCoords);

        // Compute mapping using the thread-local scratchpads
        refEl.computeMapping(xi, ctx.scratchCoords, ShapeEvaluationFlags::Gradient, 
                            ctx.scratchShape, ctx.scratchMapping);

        CfemInt numNodes = m_space->getNumActiveNodes();

        for (size_t i = 0; i < ctx.scratchMapping.physicalShapeGradients.size(); ++i)
        {
            CfemInt dofIdx = m_space->getLocalNodeIdx(cell_ids[i]);
            if (dofIdx != -1)
            {
                const Vector3D &grad_phys = ctx.scratchMapping.physicalShapeGradients[i];
                for (int c = 0; c < nComp; ++c)
                {
                    CfemReal nodalVal = m_localValues[dofIdx + (c * numNodes)];
                    // Layout: [du/dx, du/dy, du/dz, dv/dx, dv/dy, dv/dz, ...]
                    out[c * 3 + 0] += grad_phys.x * nodalVal;
                    out[c * 3 + 1] += grad_phys.y * nodalVal;
                    out[c * 3 + 2] += grad_phys.z * nodalVal;
                }
            }
        }

        // For debugging
        // for (int i=0; i < out.size()/int(3); ++i){
        //     std::cout << "Column "  <<i << ": " <<out[3*i + 0] << ", " << out[3*i + 1] << ", " << out[3*i + 2] << std::endl;
        // }
    }

    void FEField::evaluateGradientWithCachedMapping(CfemInt cellId, const MappingInfo &mapping, std::span<CfemReal> out) const
    {
        int nComp = this->getDimension();
        const auto &cell_ids = m_space->getMesh()->getCellConnectivity(cellId);
        CfemInt numNodes = m_space->getNumActiveNodes();

        //std::cout << out.size() << std::endl;

        std::fill_n(out.data(), nComp * 3, 0.0);

        for (size_t i = 0; i < mapping.physicalShapeGradients.size(); ++i)
        {
            CfemInt dofIdx = m_space->getLocalNodeIdx(cell_ids[i]);
            if (dofIdx != -1)
            {
                const Vector3D &grad_phys = mapping.physicalShapeGradients[i];
                //std::cout << grad_phys << std::endl;
                for (int c = 0; c < nComp; ++c)
                {
                    CfemReal val = m_localValues[dofIdx + (c * numNodes)];
                    out[c * 3 + 0] += grad_phys.x * val;
                    out[c * 3 + 1] += grad_phys.y * val;
                    out[c * 3 + 2] += grad_phys.z * val;
                }
            }
        }

        // For debugging
        // for (int i=0; i < out.size()/int(3); ++i){
        //     std::cout << "Column "  <<i << ": " <<out[3*i + 0] << ", " << out[3*i + 1] << ", " << out[3*i + 2] << std::endl;
        // }
    }

//     void FEField::project(const sym::CfemExpression &recipe)
//     {
//         recipe.prepare();

//         auto mesh = m_space->getMesh();
//         const int dim = getDimension();
//         const CfemInt numNodes = m_space->getNumActiveNodes();
//         auto lagrange = std::static_pointer_cast<LagrangeSpace>(m_space);

//         // Standard component-major buffer
//         std::vector<CfemReal> cpuBuffer(numNodes * dim, 0.0);

//         // --- PATH A: SUPER-FAST POINTWISE PROJECTION ---
//         // Used for analytical lambdas. Bypasses cells, mappings, and ref-elements.
//         if (recipe.isPointwise())
//         {
// #pragma omp parallel
//             {
//                 EvalContext evctx;

                
//                 CfemReal valStorage[3] = {0};
//                 std::span<CfemReal> valSpan(valStorage, dim);

// #pragma omp for
//                 for (CfemInt i = 0; i < numNodes; ++i)
//                 {
//                     evctx.resetBuffers();
//                     // We use the cached map to find the correct Mesh ID for DOF 'i'
//                     CfemInt vId = lagrange->getMeshNodeIdFromLocal(i);

//                     // Get coordinates directly from the Mesh
//                     evctx.physPos = mesh->getVertex(vId).coords;

//                     // Pointwise evaluation (Analytical)
//                     recipe.evaluate(-1, Vector3D{0, 0, 0}, valSpan, evctx);

//                     for (int d = 0; d < dim; ++d)
//                     {
//                         cpuBuffer[i + (d * numNodes)] = valSpan[d];
//                     }
//                 }
//             }
//         }
//         // --- PATH B: ROBUST CELL-BASED PROJECTION ---
//         // Used when the recipe contains other Fields. Uses 'visited' to skip shared nodes.
//         else
//         {

//             // Atomic-ready visited array (uint8_t is faster than vector<bool>)
//             std::vector<uint8_t> visited(numNodes, 0);

// #pragma omp parallel
//             {
//                 EvalContext evctx;
//                 MappingInfo map;
//                 DynamicVector<Vector3D> physCoords;
//                 // CfemReal valStorage[3] = {0};
//                 // std::span<CfemReal> valSpan(valStorage, dim);

// #pragma omp for
//                 for (CfemInt cellId = 0; cellId < mesh->getNumCells(); ++cellId)
//                 {
//                     evctx.resetBuffers();
//                     const auto &cell = mesh->getCell(cellId);
//                     const auto &refEl = m_space->getRefElement(cell.type);
//                     const auto &cellNodes = mesh->getCellConnectivity(cellId);
//                     mesh->getCellVerticesCoords(cellId, physCoords);

//                     for (size_t i = 0; i < refEl.getNumNodes(); ++i)
//                     {
//                         auto qpBookmark = evctx.getBookmark();
//                         std::span<CfemReal> spanOut = evctx.getBuffer(9);

//                         CfemInt localIdx = m_space->getLocalNodeIdx(cellNodes[i]);
//                         if (localIdx == -1)
//                             continue;

//                         // Atomic test-and-set: Only the first thread to reach this node processes it
//                         if (__atomic_test_and_set(&visited[localIdx], __ATOMIC_RELAXED))
//                         {
//                             continue;
//                         }

//                         Vector3D xi = refEl.getNodeCoords(i);
//                         ShapeInfo nodeShape;
//                         refEl.evaluate(xi, ShapeEvaluationFlags::Value, nodeShape);
//                         refEl.computeMapping(xi, physCoords, ShapeEvaluationFlags::None, nodeShape, map);

//                         evctx.physPos = map.physicalPosition;
//                         recipe.evaluate(cellId, xi, spanOut, evctx);

//                         for (int d = 0; d < dim; ++d)
//                         {
//                             cpuBuffer[localIdx + d * numNodes] = spanOut[d];
//                         }
//                         evctx.releaseToBookmark(qpBookmark);
//                     }
//                 }
//             }
//         }

//         // --- UPLOAD TO PETSC ---
//         // Instead of copying, use a pointer-based upload if your wrapper supports it

//         std::vector<CfemInt> indices(cpuBuffer.size());
//         for (CfemInt i = 0; i < indices.size(); ++i)
//             indices[i] = i;

//         VecSetValues(m_data.getRaw(), cpuBuffer.size(), indices.data(), cpuBuffer.data(), INSERT_VALUES);
//         m_data.assemble();

//         this->updateLocalBuffer();
//     }

void FEField::project(const sym::CfemExpression &recipe)
{
    CFEM_INFO << "[FEField] Projecting into '" << m_name << "' (RAM-optimized)...";
    recipe.prepare();

    auto mesh = m_space->getMesh();
    const CfemInt numNodes = m_space->getNumActiveNodes();
    const int dim = m_dimension;

    // 1. ON UTILISE DIRECTEMENT m_localValues (On supprime cpuBuffer)
    // On s'assure que m_localValues a la bonne taille
    m_localValues.assign(numNodes * dim, 0.0); 
    
    // On garde 'visited' car c'est un bitset/bool, c'est léger (30 Mo pour 30M de noeuds)
    std::vector<uint8_t> visited(numNodes, 0);

    // 2. Calcul parallèle
    #pragma omp parallel
    {
        EvalContext evctx;
        MappingInfo map;
        DynamicVector<Vector3D> physCoords;
        CfemReal valStorage[9]; 
        std::span<CfemReal> valSpan(valStorage, dim);

        #pragma omp for schedule(guided)
        for (CfemInt cellId = 0; cellId < mesh->getNumCells(); ++cellId) {
            const auto &cell = mesh->getCell(cellId);
            const auto &refEl = m_space->getRefElement(cell.type);
            const auto &cellNodes = mesh->getCellConnectivity(cellId);
            mesh->getCellVerticesCoords(cellId, physCoords);

            for (size_t i = 0; i < refEl.getNumNodes(); ++i) {
                CfemInt localIdx = m_space->getLocalNodeIdx(cellNodes[i]);
                
                if (localIdx == -1 || visited[localIdx]) continue;
                
                // Spin-lock léger pour les noeuds partagés
                if (__atomic_test_and_set(&visited[localIdx], __ATOMIC_RELAXED)) continue;

                evctx.resetBuffers();
                Vector3D xi = refEl.getNodeCoords(i);
                
                // (Calcul shape/mapping...)
                ShapeInfo nodeShape;
                refEl.evaluate(xi, ShapeEvaluationFlags::Value, nodeShape);
                refEl.computeMapping(xi, physCoords, ShapeEvaluationFlags::None, nodeShape, map);

                evctx.physPos = map.physicalPosition;
                evctx.currentCellId = cellId;
                evctx.currentXi = xi;

                recipe.evaluate(cellId, xi, valSpan, evctx);

                // ÉCRITURE DIRECTE DANS LE BUFFER FINAL
                for (int d = 0; d < dim; ++d) {
                    m_localValues[localIdx + (d * numNodes)] = valSpan[d];
                }
            }
        }
    } // Sortie OMP : TOUS les EvalContext sont détruits ici.

    // 3. LIBÉRATION DU TABLEAU VISITED
    visited.clear();
    visited.shrink_to_fit();

    // 4. SYNCHRONISATION PETSc (Depuis m_localValues)
    // Au lieu de copier cpuBuffer -> PETSc -> localBuffer,
    // on fait localBuffer -> PETSc directement.
    PetscScalar *array;
    VecGetArray(m_data.getRaw(), &array);
    std::copy(m_localValues.begin(), m_localValues.end(), array);
    VecRestoreArray(m_data.getRaw(), &array);
    m_data.assemble();

    // 5. PLUS BESOIN DE updateLocalBuffer() ! 
    // Il est déjà à jour, on économise une grosse copie et une allocation PETSc interne.

#ifdef __linux__
    malloc_trim(0); 
#endif

    CFEM_INFO << "[FEField] Projection completed. Memory flushed.";
}

} // namespace cfem