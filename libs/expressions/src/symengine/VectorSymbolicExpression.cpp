#include "symengine/VectorSymbolicExpression.h"
#include "symengine/TensorSymbolicExpression.h"
#include "symengine/SymEngineSymbols.h"

namespace cfem::sym
{
    VectorSymbolicExpression::VectorSymbolicExpression(
        std::shared_ptr<const CfemExpression> dx,
        std::shared_ptr<const CfemExpression> dy,
        std::shared_ptr<const CfemExpression> dz)
    {
        // Use dynamic_pointer_cast to safely convert from CfemExpression to ScalarSymbolicExpression
        m_components[0] = std::dynamic_pointer_cast<const ScalarSymbolicExpression>(dx);
        m_components[1] = std::dynamic_pointer_cast<const ScalarSymbolicExpression>(dy);
        m_components[2] = std::dynamic_pointer_cast<const ScalarSymbolicExpression>(dz);

        CFEM_ASSERT(m_components[0] && m_components[1] && m_components[2], 
            "VectorSymbolicExpression: All components must be of type ScalarSymbolicExpression.");
    }

    std::shared_ptr<const CfemExpression> VectorSymbolicExpression::grad() const {
        std::shared_ptr<const ScalarSymbolicExpression> comps[9];
        for(int i=0; i<3; ++i) {
            comps[i*3 + 0] = m_components[i]->diff(sym::x_rcp);
            comps[i*3 + 1] = m_components[i]->diff(sym::y_rcp);
            comps[i*3 + 2] = m_components[i]->diff(sym::z_rcp);
        }
        return std::make_shared<TensorSymbolicExpression>(comps[0], comps[1], comps[2], 
                                                          comps[3], comps[4], comps[5], 
                                                          comps[6], comps[7], comps[8]);
    }
}