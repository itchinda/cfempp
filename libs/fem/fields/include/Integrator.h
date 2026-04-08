#ifndef CFEM_INTEGRATOR_H
#define CFEM_INTEGRATOR_H

#include "CfemExpression.h"
#include "Quadrature.h"
#include "cfem_topology.h"
#include "FESpace.h" // <-- NOUVEAU: On utilise l'abstraction FESpace
#include <memory>
#include <cmath>
#include <stdexcept>
#include <set>

namespace cfem
{
    class Mesh;
    class MeshEntity;
    
    class Integrator
    {
    private:
        std::shared_ptr<const Mesh> m_mesh;
        IntegrationSchemeParams m_schemes_params;
        bool m_cacheInitialized = false;

        using EvaluationTarget = std::pair<std::shared_ptr<const sym::CfemExpression>, CfemReal*>;
        
        void integrate(const std::vector<EvaluationTarget> &targets, const MeshEntity &domain);
        
        void integrateBoundary(const std::vector<EvaluationTarget> &targets,
                               const MeshEntity                    &boundary,
                               const MeshEntity                    &volume);

        void dispatchIntegration(const std::vector<EvaluationTarget> &targets,
                                 const MeshEntity &entity,
                                 const MeshEntity &);

        inline bool needsMapping(ShapeEvaluationFlags flags)
        {
            return hasFlag(flags, ShapeEvaluationFlags::Gradient) ||
                   hasFlag(flags, ShapeEvaluationFlags::Hessian);
        }

    public:
        explicit Integrator(std::shared_ptr<const Mesh> mesh,
                            const IntegrationSchemeParams& params = {});
        
        void setIntegrationSchemeParams(const IntegrationSchemeParams& params){
            if (m_schemes_params != params) {
                m_schemes_params = params;
                m_cacheInitialized = false;
            }
        }

        void L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr, CfemReal* l2Out, CfemReal* h1SemiOut);
        void L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr, CfemReal* l2Out, CfemReal* h1SemiOut, const std::string& entity_name);
        void L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr, CfemReal* l2Out, CfemReal* h1SemiOut, CfemInt entity_id);
        void L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr, CfemReal* l2Out, CfemReal* h1SemiOut, const MeshEntity& entity);

        CfemReal computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr);
        CfemReal computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr, const std::string& entity_name);
        CfemReal computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr, CfemInt entity_id);
        CfemReal computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr, const MeshEntity& entity);
    };
}
#endif