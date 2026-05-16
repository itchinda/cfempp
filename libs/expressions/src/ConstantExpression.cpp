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

#include "ConstantExpression.h"

namespace cfem::sym
{

    std::shared_ptr<const CfemExpression> ConstantExpression::div() const {
        CfemInt dim = getNumComponents();
            
        switch (dim) {
            case 1:
                CFEM_ERROR << "Divergence operator cannot be applied to a scalar constant expression: " << toString();
            case 3:
                return zeros::Zero1D; 
            case 9:
                return zeros::Zero3D;
            default:
                CFEM_ERROR << "Divergence operator not supported for ConstantExpression of dimension " << dim;
        }
        return nullptr;
    }

    std::string ConstantExpression::toString() const { 
        CfemInt dim = getNumComponents();
        if (dim==1) return formatCoefficient(m_value[0]);

        std::ostringstream out;
        CfemInt i = 0;
        out << "[";
        for (auto val : m_value){
            if (i < dim-1) out << val << ", ";
            ++i;
        }
        out << m_value[dim-1] << "]";
        return out.str();
    }

}