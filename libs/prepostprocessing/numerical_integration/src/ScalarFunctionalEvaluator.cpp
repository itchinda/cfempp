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

#include <map>
#include <set>

#include "EvalContext.h"
#include "LagrangeSpace.h"
#include "MappingInfo.h"
#include "MappingInfoBatch.h"
#include "MappingInfoBatchView.h"
#include "MathOperators.h"
#include "Mesh.h"
#include "MeshEntity.h"
#include "ReferenceElementFactory.h"
#include "ScalarFunctionalEvaluator.h"
#include "ShapeInfo.h"
#include "ShapeInfoBatch.h"

namespace cfem
{
    ScalarFunctionalEvaluator::ScalarFunctionalEvaluator(std::shared_ptr<const Mesh>    mesh,
                           const IntegrationSchemeParams &schemes_params) : 
        m_mesh(std::move(mesh)),
        m_schemes_params(schemes_params)
    {
        if (!m_mesh)
            CFEM_THROW("ScalarFunctionalEvaluator requires a valid mesh.");
    }

    void ScalarFunctionalEvaluator::L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr,
                            CfemReal                                         *l2Out,
                            CfemReal                                         *h1SemiOut,
                            const MeshEntity                                 &mesh_entity) 
    {
        CfemReal l2_sqr = 0.0, h1semi_sqr = 0.0;
        std::vector<EvaluationTarget> targets;
        if (l2Out) targets.push_back({cfem::sym::mag_sqr(expr), &l2_sqr});
        if (h1SemiOut) targets.push_back({cfem::sym::mag_sqr(cfem::sym::grad(expr)), &h1semi_sqr});
        if (targets.empty()) return;
        dispatchIntegration(targets, mesh_entity, m_mesh->getEntityConst(UNIVERSE_ID));
        if (l2Out) *l2Out = std::sqrt(std::max<CfemReal>(0.0, l2_sqr));
        if (h1SemiOut) *h1SemiOut = std::sqrt(std::max<CfemReal>(0.0, h1semi_sqr));
    }

    void ScalarFunctionalEvaluator::L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr,
                            CfemReal                                         *l2Out, 
                            CfemReal                                         *h1SemiOut, 
                            const std::string                                &entity_name) {
        L2Norm(expr, l2Out, h1SemiOut, m_mesh->getEntityConst(entity_name));
    }

    void ScalarFunctionalEvaluator::L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr,
                                    CfemReal                                         *l2Out,
                                    CfemReal                                         *h1SemiOut,
                                    CfemInt                                           entity_id) {
        L2Norm(expr, l2Out, h1SemiOut, m_mesh->getEntityConst(entity_id));
    }

    void ScalarFunctionalEvaluator::L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr,
                                    CfemReal                                         *l2Out, 
                                    CfemReal                                         *h1SemiOut) {
        L2Norm(expr, l2Out, h1SemiOut, m_mesh->getEntityConst(UNIVERSE_ID));
    }

    CfemReal ScalarFunctionalEvaluator::computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr, 
                                         const MeshEntity &entity) {
        CfemReal result = 0.0;
        dispatchIntegration({{expr, &result}}, entity, m_mesh->getEntityConst(UNIVERSE_ID));
        return result;
    }

    CfemReal ScalarFunctionalEvaluator::computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr) {
        return computeIntegral(expr, m_mesh->getEntityConst(UNIVERSE_ID));
    }

    CfemReal ScalarFunctionalEvaluator::computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr,
                                                 const std::string                                &entity_name) {
        return computeIntegral(expr, m_mesh->getEntityConst(entity_name));
    }

    CfemReal ScalarFunctionalEvaluator::computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr,
                                                 CfemInt                                           entity_id) {
        return computeIntegral(expr, m_mesh->getEntityConst(entity_id));
    }

    // ********************************************************************************************** 

    CfemReal ScalarFunctionalEvaluator::computeFlux(const std::shared_ptr<const sym::CfemExpression>& expr, 
                                             CfemInt boundaryId) 
    {
        return computeFlux(expr, m_mesh->getEntityConst(boundaryId), m_mesh->getEntityConst(UNIVERSE_ID));
    }

    CfemReal ScalarFunctionalEvaluator::computeFlux(const std::shared_ptr<const sym::CfemExpression>& expr, 
                                             const std::string& boundary_name) 
    {
        return computeFlux(expr, m_mesh->getEntityConst(boundary_name), m_mesh->getEntityConst(UNIVERSE_ID));
    }

    CfemReal ScalarFunctionalEvaluator::computeFlux(const std::shared_ptr<const sym::CfemExpression>& expr, 
                                             const std::string& boundary_name,
                                             const std::string& volume_name) 
    {
        return computeFlux(expr, m_mesh->getEntityConst(boundary_name), m_mesh->getEntityConst(volume_name));
    }

    CfemReal ScalarFunctionalEvaluator::computeFlux(const std::shared_ptr<const sym::CfemExpression>& expr, 
                                             CfemInt boundary_id,
                                             CfemInt volume_id) 
    {
        return computeFlux(expr, m_mesh->getEntityConst(boundary_id), m_mesh->getEntityConst(volume_id));
    }

    CfemReal ScalarFunctionalEvaluator::computeFlux(const std::shared_ptr<const sym::CfemExpression>& expr, 
                                             const MeshEntity& boundary_entity) 
    {
        return computeFlux(expr, boundary_entity, m_mesh->getEntityConst(UNIVERSE_ID));
    }

    CfemReal ScalarFunctionalEvaluator::computeFlux(const std::shared_ptr<const sym::CfemExpression>& expr, 
                                            const MeshEntity& boundary_entity, 
                                            const MeshEntity& volume_entity)
    {
        if( ! m_mesh->isPartOfTheMeshEntities(boundary_entity.getId()) ){
            CFEM_THROW(boundary_entity.getName() << " is not part of the mesh.");
        }

        CfemInt entityDim = boundary_entity.getTopologicalDimension();
        CfemInt meshDim   = m_mesh->getTopologicalDimension();
        
        bool isBoundary = (entityDim == meshDim - 1);
        bool isManifoldSurface = (entityDim == 2 && meshDim == 2);

        if (!isBoundary && !isManifoldSurface) {
            CFEM_THROW("computeFlux: Flux can only be computed on a boundary (dim " 
                       << meshDim - 1 << ") or a surface manifold (dim 2). Target entity '" 
                       << boundary_entity.getName() << "' has dimension " << entityDim);
        }

        if (isBoundary) {
            if( ! m_mesh->isPartOfTheMeshEntities(volume_entity.getId()) ){
                CFEM_THROW(volume_entity.getName() << " is not part of the mesh.");
            }
            if( volume_entity.getTopologicalDimension() != meshDim) {
                CFEM_THROW("The body must have the same dimension as the mesh (" << meshDim << ").");
            }
        }

        if (expr->getNumComponents() != 3) {
            CFEM_THROW("computeFlux: The expression must be a 3D vector. "
                       << "Got dimension " << expr->getNumComponents());
        }
    
        auto fluxExpr = sym::dot(expr, sym::geometric_normal());
        CfemReal totalFlux = 0.0;
        
        if (isBoundary) {
            integrateBoundary({ EvaluationTarget{fluxExpr, &totalFlux} }, boundary_entity, volume_entity);
        } else {
            // Cas Manifold (isManifoldSurface == true) : 
            // La frontière EST le maillage (ex: la sphère), on intègre directement sur les cellules.
            integrate({ EvaluationTarget{fluxExpr, &totalFlux} }, boundary_entity);
        }

        return totalFlux;
    }

    void ScalarFunctionalEvaluator::dispatchIntegration(const std::vector<EvaluationTarget> &targets,
                                                 const MeshEntity                    &entity,
                                                 const MeshEntity                    &volume) 
    {
        if (entity.getCellIds().empty() && entity.getFacetIds().empty()) {
            for (const auto &target : targets) if (target.second) *(target.second) = 0.0;
            return;
        }

        const CfemInt meshDim = m_mesh->getTopologicalDimension();
        const CfemInt entityDim = entity.getTopologicalDimension();

        if (entityDim == meshDim) integrate(targets, entity);
        else if (entityDim == meshDim - 1) integrateBoundary(targets, entity, volume);
        else for (const auto &target : targets) if (target.second) *(target.second) = 0.0;
    }

    void ScalarFunctionalEvaluator::integrate(const std::vector<EvaluationTarget> &targets,
                                       const MeshEntity                    &entity)
    {
        if (targets.empty()) return;
        CFEM_INFO << "Integrating over entity " << entity.getName();

        // ============================================================
        // PHASE 1 : ANALYSE ET BINDING O(1)
        // ============================================================
        std::vector<size_t> activeIdx;
        ShapeEvaluationFlags fieldShapesFlags = ShapeEvaluationFlags::None;
        bool fieldsDependsOnGeoNormal = false;
        std::set<std::shared_ptr<const FESpace>> requiredSpaces;
        CfemInt geoOrder = m_mesh->getOrderInt();
        CfemInt spatialDim =m_mesh->getSpatialDimension();

        for (size_t i = 0; i < targets.size(); ++i) {
            if (!targets[i].first) continue;

            if (targets[i].first->getNumComponents() != 1) {
                CFEM_THROW("ScalarFunctionalEvaluator only accepts scalar expressions. "
                        << "Please pass each component separately.");
            }

            targets[i].first->prepare();

            if (targets[i].second) *(targets[i].second) = 0.0;
            activeIdx.push_back(i);
            fieldShapesFlags |= targets[i].first->requiredShapeFlags();
            fieldsDependsOnGeoNormal = fieldsDependsOnGeoNormal || targets[i].first->dependsOnGeometricNormal();

            targets[i].first->collectRequiredSpaces(requiredSpaces);
        }

        if (activeIdx.empty()) return;

        // Association des espaces à des IDs uniques pour le contexte
        std::map<const FESpace*, int> spaceToId;
        std::vector<std::shared_ptr<const FESpace>> idToSpace;
        CfemInt nextId = 0;
        for (const auto& sp : requiredSpaces) {
            spaceToId[sp.get()] = nextId++;
            idToSpace.push_back(sp);
        }

        for (size_t idx : activeIdx) {
            targets[idx].first->bind(spaceToId);
        }

        const bool needsFieldFirstDerivatives = hasFlag(fieldShapesFlags, ShapeEvaluationFlags::FirstDerivatives);
        const bool needsFieldSecondDerivatives = hasFlag(fieldShapesFlags, ShapeEvaluationFlags::SecondDerivatives);
        
        // Configuration des flags géométriques requis
        MappingEvaluationFlags geoMappingFlags = MappingEvaluationFlags::PhysicalPosition |
                                                MappingEvaluationFlags::Determinant |
                                                (fieldsDependsOnGeoNormal ? MappingEvaluationFlags::Normal : MappingEvaluationFlags::None);
        ShapeEvaluationFlags geoShapeFlags = transformMappingToShapeFlags(geoMappingFlags); 

        if (needsFieldFirstDerivatives) {
            geoShapeFlags |= ShapeEvaluationFlags::FirstDerivatives;
            geoMappingFlags |= MappingEvaluationFlags::InverseJacobian; 
        }

        IntegrationScheme scheme;

        // ============================================================
        // PHASE 2 : CACHING DE REFERENCE BATCH (Read-Only)
        // ============================================================
        struct SpaceCacheBatch {
            CfemInt contextIndex;
            ShapeInfoBatch shapesBatch;
            CfemInt paddedNodeStride;
        };

        struct CellCacheEntryBatch {
            const QuadratureRule* quad = nullptr;
            DynamicVector<Vector3D> xiPoints; 
            ShapeInfoBatch geoShapesBatch; 
            std::vector<SpaceCacheBatch> spaceCaches; 
        };

        std::vector<CellCacheEntryBatch> fastCache(static_cast<int>(CellType::COUNT));

        for (CellType cell_type : m_mesh->getCellTypesPresent())
        {
            CfemInt idx = static_cast<int>(cell_type);
            const auto &quad = scheme.getRule(cell_type, m_schemes_params.getNumIntgPt(cell_type));
            const auto &lQuadPts = quad.getPoints();

            fastCache[idx].quad = &quad;
            
            CfemInt nQp = static_cast<CfemInt>(lQuadPts.size());
            fastCache[idx].xiPoints.resize(nQp);
            for (CfemInt q = 0; q < nQp; ++q) {
                fastCache[idx].xiPoints[q] = lQuadPts[q].xi;
            }

            const auto &refElGeo = ReferenceElementFactory::getReferenceElement(cell_type, geoOrder);
            
            refElGeo.evaluateShapesBatch(fastCache[idx].xiPoints, geoShapeFlags, fastCache[idx].geoShapesBatch);
            fastCache[idx].geoShapesBatch.shapeMustBeUpdated = false;

            for (const auto& space : idToSpace) {
                SpaceCacheBatch spCache;
                spCache.contextIndex = spaceToId[space.get()];
                
                CfemInt spaceOrder = space->getOrderInt();
                const auto& refElSpace = ReferenceElementFactory::getReferenceElement(cell_type, spaceOrder);
                
                refElSpace.evaluateShapesBatch(fastCache[idx].xiPoints, fieldShapesFlags, spCache.shapesBatch);
                spCache.shapesBatch.shapeMustBeUpdated = false;

                spCache.paddedNodeStride = kernel_utils::computePaddedStride(spCache.shapesBatch.numNodes);
                
                fastCache[idx].spaceCaches.push_back(std::move(spCache));
            }
        }

        const std::vector<CfemInt> &cells_vIds = entity.getCellIds();

        // ============================================================
        // PHASE 3 : BOUCLE D'INTÉGRATION PARALLÈLE (100% BATCH)
        // ============================================================
#pragma omp parallel
        {
            std::vector<CfemReal> threadTotals(targets.size(), 0.0);
            EvalContext ctx;
            std::vector<CfemReal> batchOutVals; 

            // Le SEUL Owner de la mémoire géométrique
            MappingInfoBatch mapGeoBatch; 
            
            // L'instance locale de la vue géométrique (pour que son adresse survive à la boucle)
            MappingInfoBatchView geoView;

            // NOUVEAU : Un vecteur de Vues (Zéro allocation dynamique)
            std::vector<MappingInfoBatchView> spaceMappingsViews(idToSpace.size()); 
            
            DynamicVector<Vector3D> physCoords;
            CellType lastType = CellType::COUNT; 
            const ReferenceElement* refElGeo = nullptr;
            const CellCacheEntryBatch* cache = nullptr;

#pragma omp for schedule(static)
            for (CfemInt i = 0; i < entity.getNumCells(); ++i)
            {
                CfemInt cellid = cells_vIds[i];
                const auto &cell = m_mesh->getCell(cellid);

                if (cell.type != lastType) {
                    lastType = cell.type;
                    refElGeo = &ReferenceElementFactory::getReferenceElement(lastType, geoOrder);
                    cache = &fastCache[static_cast<int>(lastType)];

                    CfemInt nQp = cache->xiPoints.size();
                    
                    // Seul le Owner géométrique nécessite un redimensionnement (setup)
                    mapGeoBatch.setup(nQp, refElGeo->getNumNodes(), refElGeo->getDimension(), false);
                    
                    if (batchOutVals.size() < static_cast<size_t>(nQp)) {
                        batchOutVals.resize(nQp);
                    }
                }

                m_mesh->getCellVerticesCoords(cellid, physCoords);
                const bool isAffine = (geoOrder == 1) || m_mesh->isCellAffine(cellid);

                // -------------------------------------------------------------
                // ÉTAPE 3.1 : CALCULS BATCH GÉOMÉTRIQUES (Le Owner)
                // -------------------------------------------------------------
                mapGeoBatch.isCellAffine = isAffine; 
                
                refElGeo->computeMappingBatch(
                            cache->xiPoints, 
                            physCoords, 
                            geoMappingFlags, 
                            cache->geoShapesBatch, 
                            mapGeoBatch
                        );

                // On génère la vue géométrique pure (O(1))
                geoView = mapGeoBatch.createView();

                // -------------------------------------------------------------
                // ÉTAPE 3.2 : CONFIGURATION DES VUES METIERS ET DE L'ARENA
                // -------------------------------------------------------------
                // CRITIQUE HPC : On remet l'Arena à zéro au DÉBUT de la cellule.
                // Cela libère la mémoire des dérivées physiques de la cellule précédente.
                ctx.resetBuffers(); 
                
                CfemInt nQp = cache->xiPoints.size();

                for (size_t sIdx = 0; sIdx < cache->spaceCaches.size(); ++sIdx) {
                    const auto& spCache = cache->spaceCaches[sIdx];
                    
                    spaceMappingsViews[sIdx] = geoView;

                    if (needsFieldFirstDerivatives) {
                        
                        // L'intégrateur demande la mémoire continue à l'Arena
                        CfemInt paddedNodes = spCache.paddedNodeStride;
                        std::span<CfemReal> arenaStorage = ctx.getBuffer(spatialDim * nQp * paddedNodes);
                        
                        spaceMappingsViews[sIdx].setupPhysicalDerivativesMemory(arenaStorage, spCache.shapesBatch.numNodes, spatialDim, nQp);
                        refElGeo->computePhysicalFirstDerivativesBatch(spCache.shapesBatch, spaceMappingsViews[sIdx]);
                    }
                }

                // -------------------------------------------------------------
                // ÉTAPE 3.3 : BINDING AU CONTEXTE O(1)
                // -------------------------------------------------------------
                ctx.currentCellId = cellid;
                ctx.geometricMappingBatch = &geoView;
                ctx.batchWS.physPos = &mapGeoBatch.physicalPositions;

                // On lie les Vues au contexte
                for (size_t sIdx = 0; sIdx < cache->spaceCaches.size(); ++sIdx) {
                    ctx.setSpaceDataBatch(
                        cache->spaceCaches[sIdx].contextIndex, 
                        &cache->spaceCaches[sIdx].shapesBatch, 
                        &spaceMappingsViews[sIdx]
                    );
                }

                // -------------------------------------------------------------
                // ÉTAPE 3.4 : ÉVALUATION DES EXPRESSIONS
                // -------------------------------------------------------------
                const auto &lQuadPts = cache->quad->getPoints();
                std::span<CfemReal> spanOutBatch(batchOutVals.data(), nQp);

                for (size_t idx : activeIdx) {
                    // On pose un signet POUR PRÉSERVER les dérivées physiques déjà allouées !
                    auto bookmark = ctx.getBookmark();

                    // Évaluation de l'arbre d'expression complet
                    targets[idx].first->evaluateBatch(cellid, nQp, spanOutBatch, ctx);

                    // Intégration numérique (Réduction vectorisable)
                    CfemReal localSum = 0.0;
                    for (CfemInt q = 0; q < nQp; ++q) {
                        localSum += spanOutBatch[q] * mapGeoBatch.detJs[q] * lQuadPts[q].weight;
                    }
                    threadTotals[idx] += localSum;

                    // On libère uniquement la mémoire temporaire allouée par evaluateBatch (ex: tenseurs intermédiaires)
                    ctx.releaseToBookmark(bookmark); 
                }
            }

#pragma omp critical
            {
                for (size_t idx : activeIdx) {
                    if (targets[idx].second) *(targets[idx].second) += threadTotals[idx];
                }
            }
        }
    }

    void ScalarFunctionalEvaluator::integrateBoundary(const std::vector<EvaluationTarget> &targets,
                                           const MeshEntity                    &boundary,
                                           const MeshEntity                    &volume)
    {
        if (targets.empty() || !boundary.containsFacets() || !volume.containsCells()) return;

        // ============================================================
        // BINDING
        // ============================================================
        std::vector<size_t> activeIdx;
        ShapeEvaluationFlags fieldShapesFlags = ShapeEvaluationFlags::None;
        std::set<std::shared_ptr<const FESpace>> requiredSpaces;
        bool fieldsDependsOnGeoNormal = false;
        
        CfemInt geoOrder = m_mesh->getOrderInt();
        CfemInt spatialDim = m_mesh->getSpatialDimension();

        for (size_t i = 0; i < targets.size(); ++i) {
            if (!targets[i].first) continue;

            if (targets[i].first->getNumComponents() != 1) {
                CFEM_THROW( "ScalarFunctionalEvaluator only accepts scalar expressions. "
                           << "Please pass each component separately.");
            }

            targets[i].first->prepare();

            if (targets[i].second) *(targets[i].second) = 0.0;
            activeIdx.push_back(i);
            fieldShapesFlags |= targets[i].first->requiredShapeFlags();
            targets[i].first->collectRequiredSpaces(requiredSpaces);

            fieldsDependsOnGeoNormal = fieldsDependsOnGeoNormal || targets[i].first->dependsOnGeometricNormal();
        }
        if (activeIdx.empty()) return;

        std::map<const FESpace*, int> spaceToId;
        std::vector<std::shared_ptr<const FESpace>> idToSpace;
        CfemInt nextId = 0;
        for (const auto& sp : requiredSpaces) {
            spaceToId[sp.get()] = nextId++;
            idToSpace.push_back(sp);
        }

        for (size_t idx : activeIdx) {
            targets[idx].first->bind(spaceToId);
        }

        const bool needsFieldFirstDerivatives = hasFlag(fieldShapesFlags, ShapeEvaluationFlags::FirstDerivatives);
        const bool needsFieldSecondDerivatives = hasFlag(fieldShapesFlags, ShapeEvaluationFlags::SecondDerivatives);
        
        // Facet flags (for dS, Normal et Position) ---
        const MappingEvaluationFlags facetGeoMappingFlags = MappingEvaluationFlags::PhysicalPosition | 
                                                            MappingEvaluationFlags::Determinant |
                                                            (fieldsDependsOnGeoNormal ? MappingEvaluationFlags::Normal : MappingEvaluationFlags::None);
        const ShapeEvaluationFlags facetGeoShapeFlags = transformMappingToShapeFlags(facetGeoMappingFlags);
        
        // Parent flags (Parent geo mapping inverse jacobian for derivatives if required)
        const ShapeEvaluationFlags parentGeoShapeFlags = fieldShapesFlags;
        const MappingEvaluationFlags parentGeoMappingFlags = transformShapeToMappingFlags(fieldShapesFlags);


        IntegrationScheme scheme;

        // ============================================================
        // HIBRID BATCH CACHING (Face + Volume)
        // ============================================================
        struct BoundarySpaceCacheBatch {
            CfemInt contextIndex; 
            ShapeInfoBatch shapesBatch; // Shapes evaluated at the facet nodes projected on the parent cell
            CfemInt paddedNodeStride;   // This will serve in the MappingInfoBatchView for HPC flat data storage
        };
        
        struct BoundaryCacheEntryBatch {
            const QuadratureRule* quad = nullptr;
            DynamicVector<Vector3D> xiFacetPoints;  // Reference facet coords
            DynamicVector<Vector3D> xiParentPoints; // Projected corrds on the parent cell
            
            ShapeInfoBatch geoShapesParentBatch; 
            ShapeInfoBatch geoShapesFacetBatch;
            std::vector<BoundarySpaceCacheBatch> spaceCaches;
        };
        
        // Indexed by: [CellType][LocalFacetId]
        std::vector<std::vector<BoundaryCacheEntryBatch>> parentFaceCache(static_cast<int>(CellType::COUNT));

        for (CellType cType : m_mesh->getCellTypesPresent()) 
        {
            CfemInt cIdx = static_cast<int>(cType);
            CfemInt numFaces = Cell::getFaceCount(cType); 
            parentFaceCache[cIdx].resize(numFaces);

            const auto& refElParentGeo = ReferenceElementFactory::getReferenceElement( cType, geoOrder );

            for (CfemInt f_local = 0; f_local < numFaces; ++f_local) 
            {
                CellType fType = Cell::getFacetType(cType, f_local);
                const auto &refElFacetGeo = ReferenceElementFactory::getReferenceElement(fType, geoOrder);
                const auto &quad = scheme.getRule(fType, m_schemes_params.getNumIntgPt(fType));
                const auto &lQuadPts = quad.getPoints();

                auto& entry = parentFaceCache[cIdx][f_local];
                entry.quad = &quad;
                
                CfemInt nQp = static_cast<CfemInt>(lQuadPts.size());
                entry.xiFacetPoints.resize(nQp);
                entry.xiParentPoints.resize(nQp);
                
                // Quadrature points batch pre-projection
                for (CfemInt q = 0; q < nQp; ++q) {
                    entry.xiFacetPoints[q] = lQuadPts[q].xi;
                    entry.xiParentPoints[q] = refElParentGeo.projectFacetRefCoordToCellRefCoord(f_local, lQuadPts[q].xi);
                }

                // Batch shapes evaluations
                refElFacetGeo.evaluateShapesBatch(entry.xiFacetPoints, facetGeoShapeFlags, entry.geoShapesFacetBatch);
                entry.geoShapesFacetBatch.shapeMustBeUpdated = false;
                
                // Batch shape evaluation on the parent cell using the projected reference facet coords
                refElParentGeo.evaluateShapesBatch(entry.xiParentPoints, parentGeoShapeFlags, entry.geoShapesParentBatch);
                entry.geoShapesParentBatch.shapeMustBeUpdated = false;

                // Caching des espaces métiers
                for (size_t sIdx = 0; sIdx < idToSpace.size(); ++sIdx) {
                    BoundarySpaceCacheBatch spCache;
                    spCache.contextIndex = spaceToId[idToSpace[sIdx].get()];
                    
                    CfemInt spaceOrder = idToSpace[sIdx]->getOrderInt();
                    const auto& refElSpace = ReferenceElementFactory::getReferenceElement(cType, spaceOrder);
                    
                    // Fields/expressions are evaluated at the PROJECTED coordinates (xiParentPoints)
                    refElSpace.evaluateShapesBatch(entry.xiParentPoints, fieldShapesFlags, spCache.shapesBatch);
                    spCache.shapesBatch.shapeMustBeUpdated = false;
                    
                    spCache.paddedNodeStride = kernel_utils::computePaddedStride(spCache.shapesBatch.numNodes);
                    
                    entry.spaceCaches.push_back(std::move(spCache));
                }
            }
        }

        const std::vector<CfemInt> &facet_ids = boundary.getFacetIds();

        // ============================================================
        // INTEGRATION LOOP USING THE BATCH ARCHITECTURE
        // ============================================================
    #pragma omp parallel
        {
            std::vector<CfemReal> threadTotals(targets.size(), 0.0);
            EvalContext ctx;
            std::vector<CfemReal> batchOutVals; 

            // Two geometrics owners :one for the surface, another for the volume
            MappingInfoBatch mapFacetGeoBatch; 
            MappingInfoBatch mapParentGeoBatch; 
            
            MappingInfoBatchView geoFacetView;
            std::vector<MappingInfoBatchView> spaceMappingsViews(idToSpace.size()); 
            
            DynamicVector<Vector3D> facetNodes;
            DynamicVector<Vector3D> parentNodes;
            
            CellType lastParentType = CellType::COUNT; 
            CfemInt lastLocalFacet = -1;
            
            const ReferenceElement* refElFacetGeo = nullptr;
            const ReferenceElement* refElParentGeo = nullptr;
            const BoundaryCacheEntryBatch* cache = nullptr;

    #pragma omp for schedule(static)
            for (CfemInt i = 0; i < boundary.getNumFacets(); ++i)
            {
                CfemInt fid = facet_ids[i];
                const auto &facet = m_mesh->getFacet(fid);

                CfemInt parentId = -1;
                CfemInt localFacetId = -1; 
                
                if (facet.left_cell != -1 && volume.containsCell(facet.left_cell)) {
                    parentId = facet.left_cell;
                    localFacetId = facet.left_cell_facet_id;
                } else if (facet.right_cell != -1 && volume.containsCell(facet.right_cell)) {
                    parentId = facet.right_cell;
                    localFacetId = facet.right_cell_facet_id;
                }
                
                if (parentId == -1) continue;

                const auto &parentCell = m_mesh->getCell(parentId);
                const bool isAffine = (geoOrder == 1) || m_mesh->isCellAffine(parentId);

                // Cache and allocationb management
                if (parentCell.type != lastParentType || localFacetId != lastLocalFacet) {
                    lastParentType = parentCell.type;
                    lastLocalFacet = localFacetId;
                    
                    refElFacetGeo = &ReferenceElementFactory::getReferenceElement(facet.type, geoOrder);
                    refElParentGeo = &ReferenceElementFactory::getReferenceElement(parentCell.type, geoOrder);
                    cache = &parentFaceCache[static_cast<int>(lastParentType)][lastLocalFacet];

                    CfemInt nQp = cache->xiFacetPoints.size();
                    
                    mapFacetGeoBatch.setup(nQp, refElFacetGeo->getNumNodes(), refElFacetGeo->getDimension(), false);
                    mapParentGeoBatch.setup(nQp, refElParentGeo->getNumNodes(), refElParentGeo->getDimension(), false);
                    
                    if (batchOutVals.size() < static_cast<size_t>(nQp)) {
                        batchOutVals.resize(nQp);
                    }
                }

                m_mesh->getFacetVerticesCoords(fid, facetNodes);
                m_mesh->getCellVerticesCoords(parentId, parentNodes);

                // -------------------------------------------------------------
                // GEOMETRY BATCH EVALUATION (Surface et Volume)
                // -------------------------------------------------------------
                mapFacetGeoBatch.isCellAffine = isAffine;
                mapParentGeoBatch.isCellAffine = isAffine;

                refElFacetGeo->computeMappingBatch(cache->xiFacetPoints, facetNodes, facetGeoMappingFlags, 
                                                   cache->geoShapesFacetBatch, mapFacetGeoBatch);

                if (needsFieldFirstDerivatives || needsFieldSecondDerivatives) {
                    refElParentGeo->computeMappingBatch(cache->xiParentPoints, parentNodes, parentGeoMappingFlags, 
                                                        cache->geoShapesParentBatch, mapParentGeoBatch);
                }

                geoFacetView = mapFacetGeoBatch.createView();

                // -------------------------------------------------------------
                // HYBRID VIEW ON THE MAPPING AND ARENA
                // -------------------------------------------------------------
                ctx.resetBuffers(); 
                CfemInt nQp = cache->xiFacetPoints.size();

                for (size_t sIdx = 0; sIdx < cache->spaceCaches.size(); ++sIdx) {
                    const auto& spCache = cache->spaceCaches[sIdx];
                    
                    // Lightweight view
                    spaceMappingsViews[sIdx] = geoFacetView; 
                    
                    if (needsFieldFirstDerivatives) {
                        // Here we need the parent inverse mapping jacobian
                        spaceMappingsViews[sIdx].invJacobians = mapParentGeoBatch.invJacobians;
                        spaceMappingsViews[sIdx].validFlags |= MappingEvaluationFlags::InverseJacobian;

                        CfemInt paddedNodes = spCache.paddedNodeStride;
                        std::span<CfemReal> arenaStorage = ctx.getBuffer(spatialDim * nQp * paddedNodes);
                        
                        spaceMappingsViews[sIdx].setupPhysicalDerivativesMemory(arenaStorage, spCache.shapesBatch.numNodes, spatialDim, nQp);
                        
                        // CRITICAL: Derivatives are computed using the parent geo mapping
                        refElParentGeo->computePhysicalFirstDerivativesBatch(spCache.shapesBatch, spaceMappingsViews[sIdx]);
                    }
                }

                // -------------------------------------------------------------
                // BINDING IN THE CONTEXT
                // -------------------------------------------------------------
                ctx.currentCellId = parentId;
                ctx.geometricMappingBatch = &geoFacetView;
                ctx.batchWS.physPos = &mapFacetGeoBatch.physicalPositions;

                for (size_t sIdx = 0; sIdx < cache->spaceCaches.size(); ++sIdx) {
                    ctx.setSpaceDataBatch(
                        cache->spaceCaches[sIdx].contextIndex, 
                        &cache->spaceCaches[sIdx].shapesBatch, 
                        &spaceMappingsViews[sIdx]
                    );
                }

                // -------------------------------------------------------------
                // FIELD EVALUATION ON THE BOUNDARY
                // -------------------------------------------------------------
                const auto &lQuadPts = cache->quad->getPoints();
                std::span<CfemReal> spanOutBatch(batchOutVals.data(), nQp);

                for (size_t idx : activeIdx) {
                    auto bookmark = ctx.getBookmark();

                    targets[idx].first->evaluateBatch(parentId, nQp, spanOutBatch, ctx);

                    CfemReal localSum = 0.0;
                    for (CfemInt q = 0; q < nQp; ++q) {
                        // Le mapFacetGeoBatch.detJs[q] contient bien le 'dS' surfacique.
                        localSum += spanOutBatch[q] * mapFacetGeoBatch.detJs[q] * lQuadPts[q].weight;
                    }
                    threadTotals[idx] += localSum;

                    ctx.releaseToBookmark(bookmark); 
                }
            }

    #pragma omp critical
            {
                for (size_t idx : activeIdx) {
                    if (targets[idx].second) *(targets[idx].second) += threadTotals[idx];
                }
            }
        }
    }

} // namespace cfem
