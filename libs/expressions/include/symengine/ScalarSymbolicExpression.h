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
#include <regex>

namespace cfem::sym
{

    inline std::string finalizeString(std::string str) {
        // This regex looks for a decimal point followed by zeros at the end of a number
        // It turns 1.0 -> 1, 4.0 -> 4, and 3.14000 -> 3.14
        str = std::regex_replace(str, std::regex(R"(\.0+(?![0-9]))"), "");
        str = std::regex_replace(str, std::regex(R"((\.[0-9]*[1-9])0+(?![0-9]))"), "$1");
        return str;
    }
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
            const CfemReal input[4] = { ctx.physPos.x, ctx.physPos.y, ctx.physPos.z, ctx.currentTime};
            out[0] = static_cast<CfemReal>(m_evaluator->evaluate(input));
        }

        std::shared_ptr<const ScalarSymbolicExpression> diff( const SymEngine::RCP<const SymEngine::Symbol> &var) const {
            return std::make_shared<const ScalarSymbolicExpression>(m_expression.diff(var), m_vars);
        }

        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx ) const override;

        std::shared_ptr<const CfemExpression> timeDerivative() const override;

        std::shared_ptr<const CfemExpression> grad() const override;
        std::shared_ptr<const CfemExpression> div() const override;

        std::string toString() const override { return finalizeString( m_expression.get_basic()->__str__() ); } // e.g., "x"
        ExpressionType getType() const override { return ExpressionType::ScalarSymbolic; }


    private:
        SymEngine::Expression                                m_expression;
        std::vector<SymEngine::RCP<const SymEngine::Symbol>> m_vars;
        std::unique_ptr<SymbolicEvaluator>                   m_evaluator;   //! Compiled evaluator owned by this expression instance
    };




}