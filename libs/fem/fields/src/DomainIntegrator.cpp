#include "DomainIntegrator.h"
#include "ReferenceElementFactory.h"

namespace cfem
{

    void DomainIntegrator::integrate(
        const std::shared_ptr<LagrangeSpace> &space,
        const std::vector<EvaluationTarget> &targets,
        const IntegrationSchemeParams &params)
    {
        if (targets.empty())
            return;

        // 1. Préparation des targets
        std::vector<size_t> pointwiseIdx;
        std::vector<size_t> fieldIdx;

        for (size_t i = 0; i < targets.size(); ++i)
        {
            if (!targets[i].first)
                continue;
            targets[i].first->prepare();
            if (targets[i].second)
                *(targets[i].second) = 0.0;

            if (targets[i].first->isPointwise())
                pointwiseIdx.push_back(i);
            else
                fieldIdx.push_back(i);
        }

        auto mesh = space->getMesh();
        IntegrationScheme scheme;
        const size_t numTargets = targets.size();
        auto space_order_equivalence = space->getOrder() == Order::P1 ? MeshOrder::Linear : MeshOrder::Quadratic;

        // 2. Pré-caching des quadratures et ShapeInfo
        DynamicVector<CacheEntry> fastCache(static_cast<int>(CellType::COUNT));
        ShapeEvaluationFlags internalFields = ShapeEvaluationFlags::Value | ShapeEvaluationFlags::Gradient;

        for (CellType type : mesh->getCellTypesPresent())
        {
            int idx = static_cast<int>(type);
            const auto &quad = scheme.getRule(type, params.getNumIntgPt(type));
            const auto &refElP2 = space->getRefElement(type);
            const auto &refElP1 = space->getRefElement(Cell::getLinearEquivalent(type));
            const auto &lQuadPts = quad.getPoints();

            fastCache[idx].quad = &quad;
            fastCache[idx].v_shapesP2.resize(lQuadPts.size());
            fastCache[idx].v_shapesP1.resize(lQuadPts.size());

            for (size_t q = 0; q < lQuadPts.size(); ++q)
            {
                refElP2.evaluate(lQuadPts[q].xi, internalFields, fastCache[idx].v_shapesP2[q]);
                refElP1.evaluate(lQuadPts[q].xi, internalFields, fastCache[idx].v_shapesP1[q]);
                fastCache[idx].v_shapesP2[q].shapeMustBeUpdated = false;
                fastCache[idx].v_shapesP1[q].shapeMustBeUpdated = false;
            }
        }

        // 3. Boucle parallèle
#pragma omp parallel
        {
            std::vector<CfemReal> threadTotals(numTargets, 0.0);
            EvalContext ctx;
            CfemReal outVal[1] = {0.0};
            std::span<CfemReal> spanOut(outVal, 1);

            MappingInfo threadMapP2; // Le mapping maître (Géométrie)
            MappingInfo threadMapP1; // Le mapping esclave (Champ P1)
            DynamicVector<Vector3D> threadPhysCoords;

#pragma omp for schedule(static)
            for (CfemInt i = 0; i < mesh->getNumCells(); ++i)
            {
                // ctx.resetBuffers();
                const auto &cell = mesh->getCell(i);
                const auto &refElGeo = ReferenceElementFactory::getReferenceElement(cell.type);
                auto &cache = fastCache[static_cast<int>(cell.type)];
                const auto &lQuadPts = cache.quad->getPoints();

                mesh->getCellVerticesCoords(i, threadPhysCoords);

                for (size_t q = 0; q < lQuadPts.size(); ++q)
                {

                    // auto qpBookmark = ctx.getBookmark();
                    // std::span<CfemReal> spanOut = ctx.getBuffer(9);

                    const auto &qp = lQuadPts[q];

                    // 1. CALCUL MAÎTRE (Coûteux : Position, detJ, invJ et gradients P2)
                    threadMapP2.setup(static_cast<int>(threadPhysCoords.size()));
                    refElGeo.computeMapping(qp.xi, threadPhysCoords,
                                            ShapeEvaluationFlags::Gradient,
                                            cache.v_shapesP2[q], threadMapP2);

                    // 2. Préparation du contexte commun
                    ctx.physPos = threadMapP2.physicalPosition;
                    ctx.currentCellId = i;
                    ctx.currentXi = qp.xi;
                    ctx.shapeP2 = &cache.v_shapesP2[q];
                    ctx.mappingP2 = &threadMapP2;

                    const CfemReal wDetJ = std::abs(threadMapP2.detJ) * qp.weight;

                    // 3. OPTIMISATION : Transformation rapide pour le P1
                    if (space_order_equivalence != mesh->getOrder())
                    {
                        // Copier les propriétés géométriques (detJ, invJ, position) depuis P2
                        threadMapP1 = threadMapP2;

                        // Mettre à jour UNIQUEMENT les gradients physiques avec J^{-T}
                        threadMapP1.computePhysicalGradientsOnly(cache.v_shapesP1[q]);

                        ctx.shapeP1 = &cache.v_shapesP1[q];
                        ctx.mappingP1 = &threadMapP1;
                    }
                    else
                    {
                        // Si le maillage et le champ sont du même ordre, P1 = P2
                        ctx.shapeP1 = ctx.shapeP2;
                        ctx.mappingP1 = ctx.mappingP2;
                    }

                    // --- Évaluation (Pointwise & Fields) ---
                    for (size_t idx : pointwiseIdx)
                    {
                        targets[idx].first->evaluate(i, qp.xi, spanOut, ctx);
                        threadTotals[idx] += spanOut[0] * wDetJ;
                    }

                    for (size_t idx : fieldIdx)
                    {
                        targets[idx].first->evaluate(i, qp.xi, spanOut, ctx);
                        threadTotals[idx] += spanOut[0] * wDetJ;
                    }

                    // ctx.releaseToBookmark(qpBookmark);
                }
            }

#pragma omp critical
            {
                for (size_t t = 0; t < numTargets; ++t)
                    if (targets[t].second)
                        *(targets[t].second) += threadTotals[t];
            }
        }

        fastCache.clear();
    }

} // namespace cfem
