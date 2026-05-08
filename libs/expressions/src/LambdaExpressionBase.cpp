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

#include "LambdaExpressionBase.h"

namespace cfem::sym
{
    LambdaExpressionBase::LambdaExpressionBase(CfemInt dimension,
                                               std::shared_ptr<const CfemExpression> gradient,
                                               std::shared_ptr<const CfemExpression> timeDeriv): 
        m_dimension(dimension)
    {
        
        m_gradient = std::move(gradient);
        m_timeDeriv = std::move(timeDeriv);
        if (m_gradient){
            if (m_dimension == 1 && m_gradient->getDimension() != 3){
                CFEM_THROW("The gradient of a scalar lambda expression must be a 3D vector.");
            }
            if (m_dimension == 3 && m_gradient->getDimension() != 9){
                CFEM_THROW("The gradient of a vector lambda expression must be a 2nd order tensor (dimension 9).");
            }
        }

        if (m_timeDeriv){
            if (m_timeDeriv->getDimension() != m_dimension){
                CFEM_THROW( "The time derivative must have the same dimension ("
                           << m_dimension << ") as the lambda expression.");
            }
        }
    }
}