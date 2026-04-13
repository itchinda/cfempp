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

#include "ZeroExpression.h"

namespace cfem::sym
{

    std::shared_ptr<const CfemExpression> ZeroExpression::div() const {
        switch (m_dim) {
            case 1:
                CFEM_ERROR << "Divergence operator cannot be applied to a scalar zero expression.";
            case 3:
                return zeros::Zero1D; 
            case 9:
                return zeros::Zero3D;
            default:
                CFEM_ERROR << "Divergence operator not supported for ZeroExpression of dimension " << m_dim;
        }
        return nullptr;
    }

    std::string ZeroExpression::toString() const { 
        CfemInt dim = getDimension();
        if (dim==1) return "0";
        
        std::ostringstream out;
        out << "[";
        for (CfemInt i = 0; i< m_dim-1; ++i){
            out << "0" << ", ";
        }
        out << "0" << "]";
        return out.str();
    }
}

