#ifndef CFEM_DOMAIN_INTEGRATOR_H
#define CFEM_DOMAIN_INTEGRATOR_H

#include "CfemExpression.h"
#include "MathOperators.h"
#include "LagrangeSpace.h"
#include "Quadrature.h"
#include "cfem_topology.h" 
#include <memory>
#include <cmath>
#include <stdexcept>

namespace cfem
{
    class DomainIntegrator
    {
        struct CacheEntry
        {
            const QuadratureRule *quad = nullptr;
            DynamicVector<ShapeInfo> v_shapesP1; 
            DynamicVector<ShapeInfo> v_shapesP2; 
        };

        using EvaluationTarget = std::pair<std::shared_ptr<const sym::CfemExpression>, CfemReal*>;
        
        static void integrate(const std::shared_ptr<LagrangeSpace> &space,
                       const std::vector<EvaluationTarget>  &targets,
                       const IntegrationSchemeParams        &params = {});

    public:
            
        static void L2Norm(
                    const std::shared_ptr<LagrangeSpace> &space,
                    const std::shared_ptr<const sym::CfemExpression> &expr,
                    CfemReal* l2Out,
                    CfemReal* h1SemiOut,
                    const IntegrationSchemeParams& params = {}
                )
        {
            CfemReal l2_sqr = 0.0, h1semi_sqr = 0.0;
            std::vector<EvaluationTarget> targets;

            if (l2Out)     targets.push_back({ cfem::sym::mag_sqr(expr), &l2_sqr });
            if (h1SemiOut) targets.push_back({ cfem::sym::mag_sqr(cfem::sym::grad(expr)), &h1semi_sqr });

            if (targets.empty()) return;

            integrate(space, targets, params);

            if (l2Out)     *l2Out = std::sqrt(std::max<CfemReal>(0.0, l2_sqr));
            if (h1SemiOut) *h1SemiOut = std::sqrt(std::max<CfemReal>(0.0, h1semi_sqr));


        }
    };
}
#endif