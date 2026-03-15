
#include "symengine/ScalarSymbolicExpression.h"
#include "symengine/VectorSymbolicExpression.h"
#include "symengine/SymEngineSymbols.h"

namespace cfem::sym{

    ScalarSymbolicExpression::ScalarSymbolicExpression(
        const SymEngine::Expression                                &ex, 
        const std::vector<SymEngine::RCP<const SymEngine::Symbol>> &vars)
        : m_expression(ex), m_vars(vars)
    {
        std::vector<SymEngine::RCP<const SymEngine::Basic>> inputs;
        inputs.reserve(m_vars.size());

        for (const auto &v : m_vars)
            inputs.push_back(v);

        m_evaluator = std::make_unique<SymbolicEvaluator>(m_expression.get_basic(), inputs);

    }

    std::shared_ptr<const CfemExpression> ScalarSymbolicExpression::grad() const {
        return std::make_shared<VectorSymbolicExpression>(diff(sym::x_rcp), diff(sym::y_rcp), diff(sym::z_rcp));
    }

}
