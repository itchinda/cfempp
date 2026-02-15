#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
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
    class UnaryExpression : public CfemExpression {
    public:
        UnaryExpression(std::shared_ptr<const CfemExpression> expr, MathOp op);

        CfemInt getDimension() const override { return m_dimOut; }
        MathOp getOp() const { return m_op; }
        std::shared_ptr<const CfemExpression> getInner() const { return m_expr; }
        bool isPointwise() const override { return m_expr->isPointwise(); }

        void evaluate(CfemInt cellIdx, const Vector3D& point, 
                      std::span<CfemReal> out, EvalContext& ctx) const override;

        std::shared_ptr<const CfemExpression> grad() const override;

    private:
        std::shared_ptr<const CfemExpression> m_expr;
        MathOp m_op;
        CfemInt m_dimIn;
        CfemInt m_dimOut;
    };
}

#include "UnaryExpression.tpp"