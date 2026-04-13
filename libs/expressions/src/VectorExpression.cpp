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

#include "VectorExpression.h"
#include "TensorExpression.h"
#include "symengine/ScalarSymbolicExpression.h"
#include "symengine/VectorSymbolicExpression.h"
#include "symengine/TensorSymbolicExpression.h"
#include "BinaryExpression.h"
#include "UnaryExpression.h"
#include "ScaledExpression.h"
#include "LambdaExpression.h"
#include "GradExpression.h"

#include "MathOperators.h"

namespace cfem::sym
{
    using namespace utils;

    std::shared_ptr<const CfemExpression> VectorExpression::createGradient() const
    {
        auto g_x = m_components[0]->grad(); 
        auto g_y = m_components[1]->grad(); 
        auto g_z = m_components[2]->grad(); 

        auto C00 = extractComponent(g_x, 0); auto C01 = extractComponent(g_x, 1); auto C02 = extractComponent(g_x, 2);
        auto C10 = extractComponent(g_y, 0); auto C11 = extractComponent(g_y, 1); auto C12 = extractComponent(g_y, 2);
        auto C20 = extractComponent(g_z, 0); auto C21 = extractComponent(g_z, 1); auto C22 = extractComponent(g_z, 2);

        return std::make_shared<TensorExpression>(
            C00, C01, C02,
            C10, C11, C12,
            C20, C21, C22
        );
    }

    std::shared_ptr<const CfemExpression> VectorExpression::createTimeDerivative() const{
        return std::make_shared<const VectorExpression>( m_components[0]->timeDerivative(), m_components[1]->timeDerivative(), m_components[2]->timeDerivative() );
    }

    std::shared_ptr<const CfemExpression> VectorExpression::div() const {
        auto dux_dx = m_components[0]->diff(0); 
        auto duy_dy = m_components[1]->diff(1); 
        auto duz_dz = m_components[2]->diff(2); 

        return dux_dx + duy_dy + duz_dz;
    }
}