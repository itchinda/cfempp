//************************************************************************
// --- CFEM++ Library
// ---
// --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.
// --- ALL RIGHTS RESERVED
// ---
// --- This software is protected by international copyright laws.
// --- Use, distribution, or modification of this software in any form,
// --- source or binary, for personal, academic, or non-commercial 
// --- purposes is permitted free of charge, provided that this 
// --- copyright notice and this permission notice appear in all 
// --- copies and supporting documentation.
// ---
// --- The authors and contributors provide this code "as is" without 
// --- any express or implied warranty. In no event shall they be 
// --- held liable for any damages arising from the use of this software.
//************************************************************************

#include "symengine/VectorSymbolicExpression.h"
#include "symengine/TensorSymbolicExpression.h"
#include "symengine/SymEngineSymbols.h"
#include "BinaryExpression.h"

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

    std::shared_ptr<const CfemExpression> VectorSymbolicExpression::diff(CfemInt var_idx) const {
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
        return std::make_shared<const VectorSymbolicExpression>(m_components[0]->diff(var), m_components[1]->diff(var), m_components[2]->diff(var));
    }

    std::shared_ptr<const CfemExpression> VectorSymbolicExpression::createGradient() const {
        std::shared_ptr<const ScalarSymbolicExpression> comps[9];
        for(CfemInt i=0; i<3; ++i) {
            comps[i*3 + 0] = m_components[i]->diff(sym::x_rcp);
            comps[i*3 + 1] = m_components[i]->diff(sym::y_rcp);
            comps[i*3 + 2] = m_components[i]->diff(sym::z_rcp);
        }
        return std::make_shared<TensorSymbolicExpression>(comps[0], comps[1], comps[2], 
                                                          comps[3], comps[4], comps[5], 
                                                          comps[6], comps[7], comps[8]);
    }

    std::shared_ptr<const CfemExpression> VectorSymbolicExpression::div() const {
        auto dux_dx = m_components[0]->get_expression().diff(x_rcp); 
        auto duy_dy = m_components[1]->get_expression().diff(y_rcp);
        auto duz_dz = m_components[2]->get_expression().diff(z_rcp);

        return std::make_shared<ScalarSymbolicExpression>(dux_dx + duy_dy + duz_dz, m_components[0]->get_vars());
    }

    std::shared_ptr<const CfemExpression> VectorSymbolicExpression::createTimeDerivative() const {
        std::shared_ptr<const ScalarSymbolicExpression> comps[3];
        for (CfemInt i=0; i<3; ++i){
            comps[i] = m_components[i]->diff(sym::t_rcp);
        }
        
        return std::make_shared<VectorSymbolicExpression>(comps[0], comps[1], comps[2]);
    }
}