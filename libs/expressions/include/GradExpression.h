#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
#include "CfemExpression.h"
#include "LambdaExpression.h"
#include <span>
#include <memory>


namespace cfem::sym
{

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class GradExpression : public CfemExpression
    {
        std::shared_ptr<const CfemExpression> m_expr;

    public:
        explicit GradExpression(std::shared_ptr<const CfemExpression> expr) : m_expr(std::move(expr)) {
            if (auto lambda = std::dynamic_pointer_cast<const LambdaExpressionBase>(expr) ){
                if (!(lambda->grad()))
                    CFEM_ERROR << "LambdaExpressions are currentely not supported by the automatic differentation.";
            }
        }

        CfemInt getDimension() const override { return m_expr->getDimension() * 3; }
        void prepare() const override { m_expr->prepare(); }
        void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override;
        ExpressionType getType() const override { return ExpressionType::Gradient; }
        bool dependsOnFEField() const override { return m_expr->dependsOnFEField(); }
        bool dependsOnShape() const override { return m_expr->dependsOnShape(); }

        ShapeEvaluationFlags requiredShapeFlags() const override {
            if (!dependsOnShape())
                return ShapeEvaluationFlags::None;

            return m_expr->requiredShapeFlags() | ShapeEvaluationFlags::Gradient;
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const LagrangeSpace>>& spaces) const override {
            if (m_expr) m_expr->collectRequiredSpaces(spaces);
        }

        std::shared_ptr<const CfemExpression> timeDerivative() const override{
            return std::make_shared<GradExpression>(m_expr->timeDerivative());
        }

        std::string toString() const override{ return "grad(" + m_expr->toString() +")"; }

        std::shared_ptr<const CfemExpression> grad() const override {
            CFEM_ERROR << "GradExpression does not yet support grad(). Cannot compute spatial derivatives for GradExpression.";
            return nullptr; 
        }


    
    };


}