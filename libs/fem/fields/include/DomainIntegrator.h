#ifndef CFEM_DOMAIN_INTEGRATOR_H
#define CFEM_DOMAIN_INTEGRATOR_H

#include "CfemExpression.h"
#include "Quadrature.h"
#include "cfem_topology.h"
#include <memory>
#include <cmath>
#include <stdexcept>

namespace cfem
{
    class Mesh;
    class MeshEntity;
    class DomainIntegrator
    {
    private:
        struct CacheEntry
        {
            const QuadratureRule *quad = nullptr;
            DynamicVector<ShapeInfo> v_shapesP1; 
            DynamicVector<ShapeInfo> v_shapesP2; 
        };

        std::shared_ptr<const Mesh> m_mesh;
        IntegrationSchemeParams m_schemes_params;
        
        // This assumes the mesh doesn't change between calls.
        bool m_cacheInitialized = false;
        

        using EvaluationTarget = std::pair<std::shared_ptr<const sym::CfemExpression>, CfemReal*>;
        
        void integrate(const MeshEntity                    &domain,
                       const std::vector<EvaluationTarget> &targets);
        /**
         * @brief Integrates over facets (Boundary Integration).
         * Used when domain consists of facets (dim < mesh_dim).
         */
        void integrateBoundary(const MeshEntity                    &domain,
                               const std::vector<EvaluationTarget>& targets) {};

        void dispatchIntegration(const MeshEntity &entity,
                                 const std::vector<EvaluationTarget> &targets);

        inline bool needsMapping(ShapeEvaluationFlags flags)
        {
            return hasFlag(flags, ShapeEvaluationFlags::Gradient) ||
                hasFlag(flags, ShapeEvaluationFlags::Hessian);
        }

    public:

        explicit DomainIntegrator(std::shared_ptr<const Mesh> mesh,
                                  const IntegrationSchemeParams& params = {});
        
        void setIntegrationSchemeParams(const IntegrationSchemeParams& params){
            if (m_schemes_params != params) {
                m_schemes_params = params;
                m_cacheInitialized = false;
            }
        }

        void L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr,
                    CfemReal* l2Out,
                    CfemReal* h1SemiOut);

        void L2Norm(const std::string& entity_name,
                    const std::shared_ptr<const sym::CfemExpression> &expr,
                    CfemReal* l2Out,
                    CfemReal* h1SemiOut);

        void L2Norm(CfemInt entity_id,
                    const std::shared_ptr<const sym::CfemExpression> &expr,
                    CfemReal* l2Out,
                    CfemReal* h1SemiOut);

        void L2Norm(const MeshEntity& entity,
                    const std::shared_ptr<const sym::CfemExpression> &expr,
                    CfemReal* l2Out,
                    CfemReal* h1SemiOut);
        /**
         * @brief Exemple d'une autre fonction simplifiée rendue possible
         * Calcule l'intégrale d'une fonction (ex: le volume si expr = 1.0)
         */
        CfemReal computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr);

        CfemReal computeIntegral(const std::string& entity_name,
                                 const std::shared_ptr<const sym::CfemExpression> &expr);

        CfemReal computeIntegral(CfemInt entity_id,
                                 const std::shared_ptr<const sym::CfemExpression> &expr);

        CfemReal computeIntegral(const MeshEntity& entity,
                                 const std::shared_ptr<const sym::CfemExpression> &expr);
    };
}
#endif