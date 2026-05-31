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

#pragma once

#include <functional>
#include <symengine/basic.h>
#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/lambda_double.h>
#include <symengine/symbol.h>
#include <vector>

#include "ScalarSymbolicExpression.h"
#include "TensorSymbolicExpression.h"
#include "VectorSymbolicExpression.h"


namespace cfem::sym
{

    /**
     * @namespace cfem::sym
     * @brief Symbolic utilities using SymEngine for CFEM expressions.
     */
    // using SymEngine::Expression;
    // using SymEngine::RCP;
    // using SymEngine::Symbol;

    /** @name Predefined spatial symbols
     *  @{
     */
    inline const auto x_rcp = SymEngine::symbol("x");
    inline const auto y_rcp = SymEngine::symbol("y");
    inline const auto z_rcp = SymEngine::symbol("z");
    inline const auto t_rcp = SymEngine::symbol("t");

    namespace symengine{
        inline const SymEngine::Expression x(x_rcp);
        inline const SymEngine::Expression y(y_rcp);
        inline const SymEngine::Expression z(z_rcp);
        inline const SymEngine::Expression t(t_rcp);
    }
    


    /**
     * @brief Extracts all symbolic variables from a SymEngine expression.
     * @param ex The expression to analyze
     * @return Vector of unique symbols (RCP<const Symbol>)
     */
    std::vector<SymEngine::RCP<const SymEngine::Symbol>> extract_symbols(const SymEngine::Expression &ex);
}