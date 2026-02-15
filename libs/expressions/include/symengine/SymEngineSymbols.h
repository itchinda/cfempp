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
    using SymEngine::Expression;
    using SymEngine::RCP;
    using SymEngine::Symbol;

    /** @name Predefined spatial symbols
     *  @{
     */
    inline const auto x_rcp = SymEngine::symbol("x");
    inline const auto y_rcp = SymEngine::symbol("y");
    inline const auto z_rcp = SymEngine::symbol("z");
    inline const auto t_rcp = SymEngine::symbol("t");

    inline const SymEngine::Expression x(x_rcp);
    inline const SymEngine::Expression y(y_rcp);
    inline const SymEngine::Expression z(z_rcp);
    inline const SymEngine::Expression t(t_rcp);


    // // On définit les constantes de base une fois pour toutes
    // auto __gradX__ = std::make_shared<ConstantExpression>(std::vector({1.0, 0.0, 0.0}));
    // auto __gradY__ = std::make_shared<ConstantExpression>(std::vector({1.0, 0.0, 0.0}));
    // auto __gradZ__ = std::make_shared<ConstantExpression>(std::vector({1.0, 0.0, 0.0}));

    // // On crée nos variables "Live" branchées sur le contexte
    // auto x_ = std::make_shared<VariableExpression>(&ctx.physPos.x, 1, __gradX__);
    // auto y_ = std::make_shared<VariableExpression>(&ctx.physPos.y, 1, __gradY__);
    // auto z_ = std::make_shared<VariableExpression>(&ctx.physPos.z, 1, __gradZ__);
    // auto t_ = std::make_shared<VariableExpression>(&ctx.currentTime, 1, zeros::Zero3D);
    /** @} */



    /**
     * @brief Extracts all symbolic variables from a SymEngine expression.
     * @param ex The expression to analyze
     * @return Vector of unique symbols (RCP<const Symbol>)
     */
    std::vector<RCP<const Symbol>> extract_symbols(const Expression &ex);
}