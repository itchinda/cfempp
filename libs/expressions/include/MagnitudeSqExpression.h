#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
#include "TensorO2.h"
#include "CfemExpression.h"
#include "EvalContext.h"
#include <span>
#include <memory>

namespace cfem::sym
{


    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class MagnitudeSqExpression : public CfemExpression
    {
        std::shared_ptr<const CfemExpression> m_expr;

    public:
        explicit MagnitudeSqExpression(std::shared_ptr<const CfemExpression> expr)
            : m_expr(std::move(expr)) 
        {
        }

        CfemInt getDimension() const override { return 1; }

        void prepare() const override { m_expr->prepare(); }

        void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override;

        std::shared_ptr<const CfemExpression> getInner() const { return m_expr; }
        std::shared_ptr<const CfemExpression> grad() const override;
        std::shared_ptr<const CfemExpression> timeDerivative() const override;


        CfemInt precedence() const override { return 5; }
        ExpressionType getType() const override { return ExpressionType::MagnitudeSq; }
        std::string toString() const override { return "grad(" + m_expr->toString() + ")"; }
    };

}

#include "MagnitudeSqExpression.tpp"