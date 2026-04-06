#include "DomainIntegrator.h"
#include "ReferenceElementFactory.h"
#include "LagrangeSpace.h"
#include "MathOperators.h"
#include "Mesh.h"
#include "MeshEntity.h"

#include <set>

namespace cfem
{

    DomainIntegrator::DomainIntegrator(std::shared_ptr<const Mesh> mesh,
                                       const IntegrationSchemeParams &schemes_params) : m_mesh(std::move(mesh)), m_schemes_params(schemes_params)
    {
        if (!m_mesh)
            CFEM_ERROR << "DomainIntegrator requires a valid mesh.";
    }

    void DomainIntegrator::L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr,
                                  CfemReal *l2Out,
                                  CfemReal *h1SemiOut,
                                  const MeshEntity &mesh_entity)
    {
        CfemReal l2_sqr = 0.0, h1semi_sqr = 0.0;
        std::vector<EvaluationTarget> targets;

        if (l2Out)
            targets.push_back({cfem::sym::mag_sqr(expr), &l2_sqr});
        if (h1SemiOut)
            targets.push_back({cfem::sym::mag_sqr(cfem::sym::grad(expr)), &h1semi_sqr});

        if (targets.empty())
            return;

        dispatchIntegration(targets, mesh_entity, m_mesh->getEntityConst(UNIVERSE_ID));

        if (l2Out)
            *l2Out = std::sqrt(std::max<CfemReal>(0.0, l2_sqr));
        if (h1SemiOut)
            *h1SemiOut = std::sqrt(std::max<CfemReal>(0.0, h1semi_sqr));
    }

    void DomainIntegrator::L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr,
                                  CfemReal *l2Out,
                                  CfemReal *h1SemiOut,
                                  const std::string &entity_name)
    {
        L2Norm(expr, l2Out, h1SemiOut, m_mesh->getEntityConst(entity_name));
    }

    void DomainIntegrator::L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr,
                                  CfemReal *l2Out,
                                  CfemReal *h1SemiOut,
                                  CfemInt entity_id)
    {
        L2Norm(expr, l2Out, h1SemiOut, m_mesh->getEntityConst(entity_id));
    }

    void DomainIntegrator::L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr,
                                  CfemReal *l2Out,
                                  CfemReal *h1SemiOut)
    {
        L2Norm(expr, l2Out, h1SemiOut, m_mesh->getEntityConst(UNIVERSE_ID));
    }

    CfemReal DomainIntegrator::computeIntegral(
        const std::shared_ptr<const sym::CfemExpression> &expr,
        const MeshEntity &entity)
    {
        CfemReal result = 0.0;
        dispatchIntegration({{expr, &result}}, entity, m_mesh->getEntityConst(UNIVERSE_ID));
        return result;
    }

    CfemReal DomainIntegrator::computeIntegral(
        const std::shared_ptr<const sym::CfemExpression> &expr)
    {
        return computeIntegral(expr, m_mesh->getEntityConst(UNIVERSE_ID));
    }

    CfemReal DomainIntegrator::computeIntegral(
        const std::shared_ptr<const sym::CfemExpression> &expr,
        const std::string &entity_name)
    {
        return computeIntegral(expr, m_mesh->getEntityConst(entity_name));
    }

    CfemReal DomainIntegrator::computeIntegral(
        const std::shared_ptr<const sym::CfemExpression> &expr, 
        CfemInt entity_id)
    {
        return computeIntegral(expr, m_mesh->getEntityConst(entity_id));
    }

    void DomainIntegrator::dispatchIntegration(const std::vector<EvaluationTarget> &targets,
                                               const MeshEntity &entity,
                                               const MeshEntity &volume)
    {
        // 1. If the entity is completely empty, the integral is exactly 0.
        if (entity.getCellIds().empty() && entity.getFacetIds().empty())
        {
            for (const auto &target : targets)
            {
                if (target.second)
                {
                    *(target.second) = 0.0;
                }
            }
            return;
        }

        const CfemInt meshDim = m_mesh->getDimension();
        const CfemInt entityDim = entity.getDimension();

        // 2. Handle Volume Integration (3D cells in a 3D mesh, or 2D in 2D)
        if (entityDim == meshDim)
        {
            integrate(targets, entity);
        }
        // 3. Handle Boundary Integration (2D facets in a 3D mesh, or 1D in 2D)
        else if (entityDim == meshDim - 1)
        {
            integrateBoundary(targets, entity, volume);
        }
        // 4. Fallback for lower dimensions (Vertices, edges in 3D, etc.)
        else
        {
            // These entities have a geometric measure (volume/area/length) of zero.
            // Therefore, any continuous integral over them yields zero.
            for (const auto &target : targets)
            {
                if (target.second)
                {
                    *(target.second) = 0.0;
                }
            }
        }
    }

    /**
     * @brief Performs numerical integration of multiple expressions over the domain.
     *
     * This function handles the core volume integration engine. It loops over all cells
     * in the mesh, dynamically detects the required functional spaces (P1, P2, etc.) from
     * the target expressions, caches quadrature rules, and evaluates the integrals in parallel.
     *
     * @param mesh The computational mesh representing the domain.
     * @param targets A list of pairs containing the target symbolic expressions to evaluate
     * and the memory addresses where the results will be accumulated.
     * @param params Parameters controlling the integration scheme (e.g., number of integration points).
     *
     * @note This function :
     * - processes multiple targets simultaneously to reduce mesh traversal overhead.
     * - utilizes OpenMP for multithreaded domain parallelization.
     */
    void DomainIntegrator::integrate(
        const std::vector<EvaluationTarget> &targets,
        const MeshEntity &entity
    )
    {
        if (targets.empty())
            return;

        // ============================================================
        // EXPRESSION AND SPACE NEEDS ANALYSIS
        // ============================================================
        std::vector<size_t> activeIdx;
        ShapeEvaluationFlags fieldFlags = ShapeEvaluationFlags::None;
        std::set<std::shared_ptr<const LagrangeSpace>> requiredSpaces;

        for (size_t i = 0; i < targets.size(); ++i)
        {
            if (!targets[i].first)
                continue;

            targets[i].first->prepare();
            if (targets[i].second)
                *(targets[i].second) = 0.0;

            activeIdx.push_back(i);
            fieldFlags |= targets[i].first->requiredShapeFlags();

            targets[i].first->collectRequiredSpaces(requiredSpaces);
        }

        if (activeIdx.empty())
            return;

        // Analyze detected element orders
        bool needsP1 = false;
        bool needsP2 = false;
        for (const auto &space : requiredSpaces)
        {
            if (space->getOrder() == Order::P1)
                needsP1 = true;
            if (space->getOrder() == Order::P2)
                needsP2 = true;
        }

        // Do the fields need physical gradients (e.g., H1-semi norm)?
        const bool needsFieldGradients = hasFlag(fieldFlags, ShapeEvaluationFlags::Gradient);

        // Reference values and gradients are always required to compute the Jacobian determinant.
        const ShapeEvaluationFlags geoFlags = ShapeEvaluationFlags::Value | ShapeEvaluationFlags::Gradient;

        IntegrationScheme scheme;

        // ============================================================
        // QUADRATURE AND REFERENCE VALUE CACHING
        // ============================================================
        DynamicVector<CacheEntry> fastCache(static_cast<int>(CellType::COUNT));

        for (CellType type : m_mesh->getCellTypesPresent())
        {
            int idx = static_cast<int>(type);
            const auto &quad = scheme.getRule(type, m_schemes_params.getNumIntgPt(type));
            const auto &lQuadPts = quad.getPoints();

            fastCache[idx].quad = &quad;

            // P2 Pre-computations (Default geometry and/or P2 fields)
            // Note: If the mesh is linear, we force the use of linear geometry.
            const CellType geoType = (m_mesh->getOrder() == MeshOrder::Quadratic) ? type : Cell::getLinearEquivalent(type);
            const auto &refElGeo = ReferenceElementFactory::getReferenceElement(geoType);

            fastCache[idx].v_shapesP2.resize(lQuadPts.size());

            // If P1 is required by expressions (or if the mesh is purely linear)
            if (needsP1 || geoType != type)
            {
                fastCache[idx].v_shapesP1.resize(lQuadPts.size());
            }

            for (size_t q = 0; q < lQuadPts.size(); ++q){
                ShapeEvaluationFlags evalFlags = geoFlags | fieldFlags;
                refElGeo.evaluate(lQuadPts[q].xi, evalFlags, fastCache[idx].v_shapesP2[q]);
                fastCache[idx].v_shapesP2[q].shapeMustBeUpdated = false;

                // If a P1 field has been detected, we evaluate it on the P1 reference element
                if (needsP1 || geoType != type)
                {
                    const auto &refElP1 = ReferenceElementFactory::getReferenceElement(Cell::getLinearEquivalent(type));
                    refElP1.evaluate(lQuadPts[q].xi, evalFlags, fastCache[idx].v_shapesP1[q]);
                    fastCache[idx].v_shapesP1[q].shapeMustBeUpdated = false;
                }
            }
        }

        // get a const reference to the cells' ids list
        const std::vector<CfemInt> &cells_vIds = entity.getCellIds();

        // ============================================================
        // PARALLEL INTEGRATION LOOP
        // ============================================================
#pragma omp parallel
        {
            std::vector<CfemReal> threadTotals(targets.size(), 0.0);
            EvalContext ctx;

            CfemReal outVal[9] = {0.0};
            std::span<CfemReal> spanOut(outVal, 1);

            MappingInfo mapP2; // Master Mapping (Geometry)
            MappingInfo mapP1; // Slave Mapping (P1 field)
            DynamicVector<Vector3D> physCoords;

#pragma omp for schedule(static)
            for (CfemInt i = 0; i < entity.getNumCells(); ++i)
            {
                CfemInt cellid = cells_vIds[i];
                const auto &cell = m_mesh->getCell(cellid);

                const CellType geoType = (m_mesh->getOrder() == MeshOrder::Quadratic) ? cell.type : Cell::getLinearEquivalent(cell.type);
                const auto &refElGeo = ReferenceElementFactory::getReferenceElement(geoType);

                const bool cellIsLinear = (geoType == Cell::getLinearEquivalent(cell.type));

                auto &cache = fastCache[static_cast<int>(cell.type)];
                const auto &lQuadPts = cache.quad->getPoints();

                m_mesh->getCellVerticesCoords(cellid, physCoords);

                // Pre-allocation
                mapP2.setup(static_cast<int>(physCoords.size()));

                if (needsP1 && geoType != Cell::getLinearEquivalent(cell.type))
                {
                    mapP1.setup(static_cast<int>(cache.v_shapesP1[0].values.size()));
                }

                for (size_t q = 0; q < lQuadPts.size(); ++q)
                {
                    const auto &qp = lQuadPts[q];

                    refElGeo.computeMapping(
                        qp.xi,
                        physCoords,
                        needsFieldGradients ? ShapeEvaluationFlags::Gradient : ShapeEvaluationFlags::None,
                        static_cast<const ShapeInfo &>(cache.v_shapesP2[q]),
                        mapP2);

                    const CfemReal wDetJ = std::abs(mapP2.detJ) * qp.weight;

                    // CONTEXT CONFIGURATION
                    ctx.currentCellId = cellid;
                    ctx.currentXi = qp.xi;
                    ctx.physPos = mapP2.physicalPosition;

                    ctx.shapeP2 = &cache.v_shapesP2[q];
                    ctx.mappingP2 = &mapP2;

                    // SLAVE MAPPING: ONLY IF NECESSARY
                    if (needsP1)
                    {
                        if (cellIsLinear)
                        {
                            // If the geometry is already P1, no transformation is needed
                            ctx.shapeP1 = ctx.shapeP2;
                            ctx.mappingP1 = ctx.mappingP2;
                        }
                        else
                        {
                            // Light copy without vector reallocation
                            mapP1.detJ = mapP2.detJ;
                            mapP1.invJacobian = mapP2.invJacobian;
                            mapP1.physicalPosition = mapP2.physicalPosition;

                            if (needsFieldGradients)
                            {
                                // Only perform matrix transformation (J^-T * dN) if required
                                mapP1.computePhysicalGradientsOnly(cache.v_shapesP1[q]);
                            }

                            ctx.shapeP1 = &cache.v_shapesP1[q];
                            ctx.mappingP1 = &mapP1;
                        }
                    }
                    else
                    {
                        // Fallback: point to P2 to avoid nullptrs
                        ctx.shapeP1 = ctx.shapeP2;
                        ctx.mappingP1 = ctx.mappingP2;
                    }

                    // EVALUATION AND SUMMATION
                    for (size_t idx : activeIdx)
                    {
                        targets[idx].first->evaluate(i, qp.xi, spanOut, ctx);
                        threadTotals[idx] += spanOut[0] * wDetJ;
                    }
                }
            }

            // CRITICAL REDUCTION OF LOCAL THREAD TOTALS
#pragma omp critical
            {
                for (size_t idx : activeIdx)
                {
                    if (targets[idx].second)
                    {
                        *(targets[idx].second) += threadTotals[idx];
                    }
                }
            }
        }
    }


    void DomainIntegrator::integrateBoundary(const std::vector<EvaluationTarget> &targets,
                                         const MeshEntity                        &boundary,
                                         const MeshEntity                        &volume)
    {
        if (targets.empty() || !boundary.containsFacets() || !volume.containsCells()) return;

        // ============================================================
        // 1. EXPRESSION AND SPACE NEEDS ANALYSIS
        // ============================================================
        std::vector<size_t> activeIdx;
        ShapeEvaluationFlags fieldFlags = ShapeEvaluationFlags::None;
        std::set<std::shared_ptr<const LagrangeSpace>> requiredSpaces;

        for (size_t i = 0; i < targets.size(); ++i) {
            if (!targets[i].first) continue;
            targets[i].first->prepare();
            if (targets[i].second) *(targets[i].second) = 0.0;
            activeIdx.push_back(i);
            fieldFlags |= targets[i].first->requiredShapeFlags();
            targets[i].first->collectRequiredSpaces(requiredSpaces);
        }
        if (activeIdx.empty()) return;

        bool needsP1 = false;
        for (const auto &space : requiredSpaces) {
            if (space->getOrder() == Order::P1) needsP1 = true;
        }

        const bool needsFieldGradients = hasFlag(fieldFlags, ShapeEvaluationFlags::Gradient);
        
        // Le parent a besoin du Gradient SI l'expression le demande.
        // La géométrie du parent a TOUJOURS besoin de la valeur pour interpoler la position.
        const ShapeEvaluationFlags geoFlags = ShapeEvaluationFlags::Value | ShapeEvaluationFlags::Gradient;
        IntegrationScheme scheme;

        // ============================================================
        // 2. CACHING PAR FACE LOCALE DU PARENT
        // ============================================================
        struct BoundaryCacheEntry {
            const QuadratureRule* quad = nullptr;
            std::vector<ShapeInfo> v_shapesP2;
            std::vector<ShapeInfo> v_shapesP1;
        };
        
        std::vector<std::vector<BoundaryCacheEntry>> parentFaceCache(static_cast<int>(CellType::COUNT));
        
        for (CellType cType : m_mesh->getCellTypesPresent()) 
        {
            int cIdx = static_cast<int>(cType);
            CfemInt numFaces = Cell::getLocalFacesCount(cType); 
            parentFaceCache[cIdx].resize(numFaces);

            const auto& refElParent = ReferenceElementFactory::getReferenceElement(cType);
            const auto& refElP1Parent = ReferenceElementFactory::getReferenceElement(Cell::getLinearEquivalent(cType));

            for (int f_local = 0; f_local < numFaces; ++f_local) 
            {
                CellType fType = Cell::getFacetType(cType, f_local);
                
                const auto &quad = scheme.getRule(fType, m_schemes_params.getNumIntgPt(fType));
                const auto &lQuadPts = quad.getPoints();

                auto& entry = parentFaceCache[cIdx][f_local];
                entry.quad = &quad;
                entry.v_shapesP2.resize(lQuadPts.size());
                if (needsP1) entry.v_shapesP1.resize(lQuadPts.size());

                for (size_t q = 0; q < lQuadPts.size(); ++q) 
                {
                    Vector3D xi_parent = refElParent.projectFacetRefCoordToCellRefCoord(f_local, lQuadPts[q].xi);
                    
                    // On active le drapeau de gradient sur le parent SI ET SEULEMENT SI c'est requis
                    ShapeEvaluationFlags evalFlags = geoFlags | fieldFlags;
                    
                    refElParent.evaluate(xi_parent, evalFlags, entry.v_shapesP2[q]);
                    entry.v_shapesP2[q].shapeMustBeUpdated = false;

                    if (needsP1) {
                        refElP1Parent.evaluate(xi_parent, evalFlags, entry.v_shapesP1[q]);
                        entry.v_shapesP1[q].shapeMustBeUpdated = false;
                    }
                }
            }
        }

        const std::vector<CfemInt> &facet_ids = boundary.getFacetIds();

        // ============================================================
        // 3. BOUCLE DE CALCUL PARALLÈLE
        // ============================================================
    #pragma omp parallel
        {
            std::vector<CfemReal> threadTotals(targets.size(), 0.0);
            EvalContext ctx;
            CfemReal outVal[9] = {0.0};
            std::span<CfemReal> spanOut(outVal, 1);

            MappingInfo mapSurf;   
            MappingInfo mapParent; 
            DynamicVector<Vector3D> facetNodes;
            DynamicVector<Vector3D> parentNodes;

            // Stockage local pour les fonctions de forme de la facette
            ShapeInfo facetShape; 

    #pragma omp for schedule(static)
            for (CfemInt i = 0; i < boundary.getNumFacets(); ++i)
            {
                CfemInt fid = facet_ids[i];
                const auto &facet = m_mesh->getFacet(fid);

                CfemInt parentId = -1;
                int localFacetId = -1; 
                
                if (facet.left_cell != -1 && volume.containsCell(facet.left_cell)) {
                    parentId = facet.left_cell;
                    localFacetId = facet.left_cell_facet_id;
                } else if (facet.right_cell != -1 && volume.containsCell(facet.right_cell)) {
                    parentId = facet.right_cell;
                    localFacetId = facet.right_cell_facet_id;
                }
                
                if (parentId == -1) continue;

                const auto &parentCell = m_mesh->getCell(parentId);
                const CellType geoTypeParent = (m_mesh->getOrder() == MeshOrder::Quadratic) ? parentCell.type : Cell::getLinearEquivalent(parentCell.type);
                const bool parentIsLinear = (geoTypeParent == Cell::getLinearEquivalent(parentCell.type));

                const auto& cache = parentFaceCache[static_cast<int>(parentCell.type)][localFacetId];
                const auto &lQuadPts = cache.quad->getPoints();

                m_mesh->getFacetVerticesCoords(fid, facetNodes);
                m_mesh->getCellVerticesCoords(parentId, parentNodes);

                mapSurf.setup(static_cast<int>(facetNodes.size()));
                mapParent.setup(static_cast<int>(parentNodes.size()));

                const auto &refElFacet = ReferenceElementFactory::getReferenceElement(facet.type);
                const auto &refElParent = ReferenceElementFactory::getReferenceElement(parentCell.type);

                for (size_t q = 0; q < lQuadPts.size(); ++q)
                {
                    const auto &qp = lQuadPts[q];

                    // --------------------------------------------------------
                    // A. MAPPING DE SURFACE
                    // --------------------------------------------------------
                    refElFacet.evaluate(qp.xi, ShapeEvaluationFlags::Value | ShapeEvaluationFlags::Gradient, facetShape);
                    
                    // On délègue tout le calcul du Jacobien et du detJ (via dim=2) à ton moteur
                    refElFacet.computeMapping(qp.xi, facetNodes, ShapeEvaluationFlags::Gradient, 
                                            facetShape, mapSurf);
                                        
                    const CfemReal dS = mapSurf.detJ * qp.weight;

                    // --------------------------------------------------------
                    // B. MAPPING DU PARENT (C'est lui qui fournit le gradient 3D !)
                    // --------------------------------------------------------
                    Vector3D xi_parent = refElParent.projectFacetRefCoordToCellRefCoord(localFacetId, qp.xi);

                    refElParent.computeMapping(
                        xi_parent, 
                        parentNodes, 
                        needsFieldGradients? ShapeEvaluationFlags::Gradient : ShapeEvaluationFlags::None,
                        cache.v_shapesP2[q], 
                        mapParent
                    );
                    
                    // --------------------------------------------------------
                    // C. CONFIGURATION DU CONTEXTE
                    // --------------------------------------------------------
                    ctx.currentCellId = parentId;
                    ctx.currentXi = xi_parent; 
                    ctx.physPos = mapSurf.physicalPosition; // Position physique lue sur la face

                    ctx.shapeP2 = &cache.v_shapesP2[q];
                    ctx.mappingP2 = &mapParent; // Mapping contenant les dérivées 3D complètes

                    // Gestion P1
                    if (needsP1) {
                        if (parentIsLinear) {
                            ctx.shapeP1 = ctx.shapeP2;
                            ctx.mappingP1 = ctx.mappingP2;
                        } else {
                            mapParent.detJ = mapParent.detJ;
                            mapParent.invJacobian = mapParent.invJacobian;
                            mapParent.physicalPosition = mapParent.physicalPosition;

                            if (needsFieldGradients) {
                                mapParent.computePhysicalGradientsOnly(cache.v_shapesP1[q]);
                            }
                            ctx.shapeP1 = &cache.v_shapesP1[q];
                            ctx.mappingP1 = &mapParent;
                        }
                    } else {
                        ctx.shapeP1 = ctx.shapeP2;
                        ctx.mappingP1 = ctx.mappingP2;
                    }

                    // --------------------------------------------------------
                    // D. EVALUATION
                    // --------------------------------------------------------
                    for (size_t idx : activeIdx) {
                        targets[idx].first->evaluate(parentId, xi_parent, spanOut, ctx);
                        threadTotals[idx] += spanOut[0] * dS;
                    }
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