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
            if (auto lambda = std::dynamic_pointer_cast<const LambdaExpression>(expr) ){
                if (!(lambda->grad()))
                    CFEM_ERROR("LambdaExpressions are currentely not supported by the automatic differentation.");
            }
        }

        CfemInt getDimension() const override { return m_expr->getDimension() * 3; }
        void prepare() const override { m_expr->prepare(); }
        void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override;
    };


}