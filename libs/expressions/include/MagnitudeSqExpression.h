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
        CfemInt m_innerDim = 0;
        std::shared_ptr<const CfemExpression> m_expr;
        

    public:
        explicit MagnitudeSqExpression(std::shared_ptr<const CfemExpression> expr)
            : m_innerDim(expr->getDimension()), m_expr(std::move(expr))
        {
            if (!m_expr) CFEM_ERROR << "MagnitudeSqExpression requires a valid inner expression.";
        }

        CfemInt getDimension() const override { return 1; }
        bool dependsOnFEField() const override { return m_expr->dependsOnFEField(); }
        bool dependsOnShape() const override { return m_expr->dependsOnShape(); }
        ShapeEvaluationFlags requiredShapeFlags() const override {
            return m_expr->requiredShapeFlags();
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            if (m_expr) m_expr->collectRequiredSpaces(spaces);
        }

        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            m_expr->bind(spaceToId);
        }

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