#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
#include "CfemExpression.h"
#include "EvalContext.h"
#include <vector>
#include <span>
#include <memory>

#include <symengine/expression.h>
#include <symengine/lambda_double.h>
#include <symengine/symbol.h>
#include <symengine/functions.h>
#include <symengine/basic.h>

#include "SymbolicEvaluator.h"


namespace cfem::sym
{

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class ScalarSymbolicExpression : public CfemExpression
    {
    public:

        ScalarSymbolicExpression(const SymEngine::Expression                               &ex, 
                                const std::vector<SymEngine::RCP<const SymEngine::Symbol>> &vars);
        CfemInt getDimension() const override { return 1; }
        bool isPointwise() const override { return true; }

        const SymEngine::Expression &get_expression() const { return m_expression; }
        const std::vector<SymEngine::RCP<const SymEngine::Symbol>> &get_vars() const { return m_vars; }

        void evaluate(CfemInt, const Vector3D &, std::span<CfemReal> out, EvalContext &ctx) const{
            const CfemReal input[3] = { ctx.physPos.x, ctx.physPos.y, ctx.physPos.z };
            out[0] = static_cast<CfemReal>(m_evaluator->evaluate(input));
        }

        std::shared_ptr<const ScalarSymbolicExpression> diff(
            const SymEngine::RCP<const SymEngine::Symbol> &var) const {
            return std::make_shared<ScalarSymbolicExpression>(m_expression.diff(var), m_vars);
        }

        std::shared_ptr<const CfemExpression> grad() const override;

    private:
        SymEngine::Expression                                m_expression;
        std::vector<SymEngine::RCP<const SymEngine::Symbol>> m_vars;
        std::unique_ptr<SymbolicEvaluator>                   m_evaluator;   //! Compiled evaluator owned by this expression instance
    };

}