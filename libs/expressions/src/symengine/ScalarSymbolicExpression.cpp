
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

    std::shared_ptr<const CfemExpression> ScalarSymbolicExpression::diff(CfemInt var_idx ) const{
        if (var_idx == 0) { return diff(x_rcp); }
        else if (var_idx == 1){ return diff(y_rcp); }
        else if (var_idx == 2){ return diff(z_rcp); }
        else if (var_idx == 3){ return diff(t_rcp); }
        else{
            CFEM_ERROR << "index out of bounds. Use 0, 1, 2, and 3, respectively for the derivative with respect to x, y, z, and t.";
            return nullptr;
        }
    }

    std::shared_ptr<const CfemExpression> ScalarSymbolicExpression::createGradient() const {
        return std::make_shared<VectorSymbolicExpression>(diff(sym::x_rcp), diff(sym::y_rcp), diff(sym::z_rcp));
    }

    std::shared_ptr<const CfemExpression> ScalarSymbolicExpression::div() const{
        CFEM_ERROR << "Divergence operator cannot be applied to a scalar field: " << toString();
        return nullptr;
    }

    std::shared_ptr<const CfemExpression> ScalarSymbolicExpression::createTimeDerivative() const {
        auto dExprDt = m_expression.diff(t_rcp);
        
        return std::make_shared<ScalarSymbolicExpression>(dExprDt, m_vars);
    }
}
