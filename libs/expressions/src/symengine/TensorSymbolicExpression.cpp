#include "symengine/TensorSymbolicExpression.h"

namespace cfem::sym
{
    TensorSymbolicExpression::TensorSymbolicExpression(
        std::shared_ptr<const ScalarSymbolicExpression> dxx, std::shared_ptr<const ScalarSymbolicExpression> dxy, std::shared_ptr<const ScalarSymbolicExpression> dxz,
        std::shared_ptr<const ScalarSymbolicExpression> dyx, std::shared_ptr<const ScalarSymbolicExpression> dyy, std::shared_ptr<const ScalarSymbolicExpression> dyz,
        std::shared_ptr<const ScalarSymbolicExpression> dzx, std::shared_ptr<const ScalarSymbolicExpression> dzy, std::shared_ptr<const ScalarSymbolicExpression> dzz)
    {
        m_components[0] = std::move(dxx);
        m_components[1] = std::move(dxy);
        m_components[2] = std::move(dxz);
        m_components[3] = std::move(dyx);
        m_components[4] = std::move(dyy);
        m_components[5] = std::move(dyz);
        m_components[6] = std::move(dzx);
        m_components[7] = std::move(dzy);
        m_components[8] = std::move(dzz);
    }

    // std::shared_ptr<const CfemExpression> TensorSymbolicExpression::grad() const  //TODO:  grad implementation goes here.
}