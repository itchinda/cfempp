#include "symengine/TensorSymbolicExpression.h"
#include "symengine/SymEngineSymbols.h"

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

    std::shared_ptr<const CfemExpression> TensorSymbolicExpression::diff(CfemInt var_idx) const {
        auto var = x_rcp;
        if (var_idx == 1){
            var = y_rcp;
        }
        else if (var_idx == 2){
            var = z_rcp;
        }
        else if (var_idx == 3){
            var = t_rcp;
        }else {
            if(var_idx != 0) CFEM_ERROR << "Index out of bounds.";
        }

        return std::make_shared<const TensorSymbolicExpression>(m_components[0]->diff(var), m_components[1]->diff(var), m_components[2]->diff(var),
                                                                m_components[3]->diff(var), m_components[4]->diff(var), m_components[5]->diff(var),
                                                                m_components[6]->diff(var), m_components[7]->diff(var), m_components[8]->diff(var));
    }

    std::shared_ptr<const CfemExpression> TensorSymbolicExpression::div() const {
        // Row 1 (X component of the resulting vector)
        auto div_x = m_components[0]->get_expression().diff(x_rcp) +
                     m_components[1]->get_expression().diff(y_rcp) +
                     m_components[2]->get_expression().diff(z_rcp);

        // Row 2 (Y component of the resulting vector)
        auto div_y = m_components[3]->get_expression().diff(x_rcp) +
                     m_components[4]->get_expression().diff(y_rcp) +
                     m_components[5]->get_expression().diff(z_rcp);

        // Row 3 (Z component of the resulting vector)
        auto div_z = m_components[6]->get_expression().diff(x_rcp) +
                     m_components[7]->get_expression().diff(y_rcp) +
                     m_components[8]->get_expression().diff(z_rcp);

        return std::make_shared<VectorSymbolicExpression>( 
            std::make_shared<ScalarSymbolicExpression>(div_x, m_components[0]->get_vars()), 
            std::make_shared<ScalarSymbolicExpression>(div_y, m_components[0]->get_vars()), 
            std::make_shared<ScalarSymbolicExpression>(div_z, m_components[0]->get_vars())
        );
    }

    // std::shared_ptr<const CfemExpression> TensorSymbolicExpression::grad() const  //TODO:  grad implementation goes here.

    std::shared_ptr<const CfemExpression> TensorSymbolicExpression::timeDerivative() const {
        std::shared_ptr<const ScalarSymbolicExpression> comps[9];
        for (CfemInt i=0; i<9; ++i){
            comps[i] = m_components[i]->diff(t_rcp);
        }
        
        return std::make_shared<TensorSymbolicExpression>(comps[0], comps[1], comps[2],
                                                          comps[3], comps[4], comps[5],
                                                          comps[6], comps[7], comps[8]);
    }
}