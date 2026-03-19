#pragma once

#include "ScalarSymbolicExpression.h"
#include "VectorSymbolicExpression.h"
#include "TensorSymbolicExpression.h"

#include <vector>
#include <functional>

#include <symengine/expression.h>
#include <symengine/lambda_double.h>
#include <symengine/symbol.h>
#include <symengine/functions.h>
#include <symengine/basic.h>


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