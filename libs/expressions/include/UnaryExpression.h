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

        bool isContinuous() const override{ 
            return m_expr->isContinuous();
        }

        bool dependsOnFEField() const override { return m_expr->dependsOnFEField(); }
        bool dependsOnShape() const override { return m_expr->dependsOnShape(); }
        
        ShapeEvaluationFlags requiredShapeFlags() const override {
            return m_expr->requiredShapeFlags();
        }

        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            if (m_expr) m_expr->bind(spaceToId);
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            if (m_expr) m_expr->collectRequiredSpaces(spaces);
        }

        void evaluate(CfemInt cellIdx, const Vector3D& point, 
                      std::span<CfemReal> out, EvalContext& ctx) const override;

        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const override;
        std::shared_ptr<const CfemExpression> grad() const override;
        std::shared_ptr<const CfemExpression> div() const override;
        std::shared_ptr<const CfemExpression> timeDerivative() const override;

        std::string toString() const override;
        CfemInt precedence() const override;
        ExpressionType getType() const override { return ExpressionType::Unary; }

    private:
        std::shared_ptr<const CfemExpression> m_expr;
        MathOp m_op;
        CfemInt m_dimIn;
        CfemInt m_dimOut;

        std::shared_ptr<const CfemExpression> getDerivativeWrtInner() const;
    };
}

#include "UnaryExpression.tpp"