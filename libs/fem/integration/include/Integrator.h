#ifndef CFEM_INTEGRATOR_H
#define CFEM_INTEGRATOR_H

#include "cfemutils.h"
#include "LagrangeSpace.h"
#include "CfemExpression.h"
#include "SymbolicExpression.h"
#include "MeshEntity.h"
#include "Quadrature.h"

#include <memory>
#include <vector>

namespace cfem
{
    class Integrator
    {
        using L2NormTargets = std::tuple<const std::shared_ptr<const CfemExpression>&, CfemReal*, CfemReal*>;
        using IntegralTargets = std::tuple<const std::shared_ptr<const CfemExpression>&, CfemReal*>;
    public:
        static void L2Norm(
            const std::shared_ptr<Mesh> &mesh,
            std::vector<L2NormTargets>& targets,
            const IntegrationSchemeParams &schemes_params = {},
            const MeshEntity* entity = nullptr
        );

        static void L2Norm(
            const std::shared_ptr<Mesh> &mesh,
            const MeshEntity* entity,
            const std::shared_ptr<const CfemExpression> &expr,
            CfemReal* l2norm,
            CfemReal* h1seminnorm,
            const IntegrationSchemeParams &schemes_params = {}
        );

        static void Integrate(); // for a general scalar function

    private:

        void integrate_internal(
            const std::shared_ptr<Mesh> &mesh,
            const std::shared_ptr<const CfemExpression> &expr,
            const IntegrationSchemeParams &schemes_params = {},
            const MeshEntity &entity);
    };

}
#endif