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

    void DomainIntegrator::L2Norm(const MeshEntity &mesh_entity,
                                  const std::shared_ptr<const sym::CfemExpression> &expr,
                                  CfemReal *l2Out,
                                  CfemReal *h1SemiOut)
    {
        CfemReal l2_sqr = 0.0, h1semi_sqr = 0.0;
        std::vector<EvaluationTarget> targets;

        if (l2Out)
            targets.push_back({cfem::sym::mag_sqr(expr), &l2_sqr});
        if (h1SemiOut)
            targets.push_back({cfem::sym::mag_sqr(cfem::sym::grad(expr)), &h1semi_sqr});

        if (targets.empty())
            return;

        dispatchIntegration(mesh_entity, targets);

        if (l2Out)
            *l2Out = std::sqrt(std::max<CfemReal>(0.0, l2_sqr));
        if (h1SemiOut)
            *h1SemiOut = std::sqrt(std::max<CfemReal>(0.0, h1semi_sqr));
    }

    void DomainIntegrator::L2Norm(const std::string &entity_name,
                                  const std::shared_ptr<const sym::CfemExpression> &expr,
                                  CfemReal *l2Out,
                                  CfemReal *h1SemiOut)
    {
        L2Norm(m_mesh->getEntityConst(entity_name), expr, l2Out, h1SemiOut);
    }

    void DomainIntegrator::L2Norm(CfemInt entity_id,
                                  const std::shared_ptr<const sym::CfemExpression> &expr,
                                  CfemReal *l2Out,
                                  CfemReal *h1SemiOut)
    {
        L2Norm(m_mesh->getEntityConst(entity_id), expr, l2Out, h1SemiOut);
    }

    void DomainIntegrator::L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr,
                                  CfemReal *l2Out,
                                  CfemReal *h1SemiOut)
    {
        L2Norm(m_mesh->getEntityConst("Universe"), expr, l2Out, h1SemiOut);
    }

    CfemReal DomainIntegrator::computeIntegral(
        const MeshEntity &entity,
        const std::shared_ptr<const sym::CfemExpression> &expr)
    {
        CfemReal result = 0.0;
        dispatchIntegration(entity, {{expr, &result}});
        return result;
    }

    CfemReal DomainIntegrator::computeIntegral(
        const std::shared_ptr<const sym::CfemExpression> &expr)
    {
        return computeIntegral(m_mesh->getEntityConst("Universe"), expr);
    }

    CfemReal DomainIntegrator::computeIntegral(
        const std::string &entity_name,
        const std::shared_ptr<const sym::CfemExpression> &expr)
    {
        return computeIntegral(m_mesh->getEntityConst(entity_name), expr);
    }

    CfemReal DomainIntegrator::computeIntegral(
        CfemInt entity_id,
        const std::shared_ptr<const sym::CfemExpression> &expr)
    {
        return computeIntegral(m_mesh->getEntityConst(entity_id), expr);
    }

    void DomainIntegrator::dispatchIntegration(const MeshEntity &entity,
                                               const std::vector<EvaluationTarget> &targets)
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
            integrate(entity, targets);
        }
        // 3. Handle Boundary Integration (2D facets in a 3D mesh, or 1D in 2D)
        else if (entityDim == meshDim - 1)
        {
            integrateBoundary(entity, targets);
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
        const MeshEntity &entity,
        const std::vector<EvaluationTarget> &targets)
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
                // "Geometry" and P2 are often merged if using isoparametric elements,
                // but here we ensure we evaluate on the correct reference element.
                refElGeo.evaluate(lQuadPts[q].xi, geoFlags | fieldFlags, fastCache[idx].v_shapesP2[q]);
                fastCache[idx].v_shapesP2[q].shapeMustBeUpdated = false;

                // If a P1 field has been detected, we evaluate it on the P1 reference element
                if (needsP1 || geoType != type)
                {
                    const auto &refElP1 = ReferenceElementFactory::getReferenceElement(Cell::getLinearEquivalent(type));
                    refElP1.evaluate(lQuadPts[q].xi, geoFlags | fieldFlags, fastCache[idx].v_shapesP1[q]);
                    fastCache[idx].v_shapesP1[q].shapeMustBeUpdated = false;
                }
            }
        }

        // CFEM_INFO << fastCache.size()  << " HHHH";

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
                    ctx.currentCellId = i;
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
        // fastCache.clear();
    }

} // namespace cfem