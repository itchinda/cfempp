#pragma once

#include "BinaryExpression.h"
#include "ConstantExpression.h"

namespace cfem::sym
{

    inline void MagnitudeSqExpression::evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        auto bookmark = ctx.getBookmark(); 

        int dim = m_expr->getDimension();
        
        // On utilise le scratchpad au lieu de allouer sur la pile/tas
        std::span<CfemReal> spanBuf = ctx.getBuffer(dim);

        m_expr->evaluate(cellId, xi, spanBuf, ctx);

        CfemReal sumSq = 0.0;
        for (int i = 0; i < dim; ++i){
            sumSq += spanBuf[i] * spanBuf[i];
        }
        out[0] = sumSq;

        ctx.releaseToBookmark(bookmark); 
    }

        // --- MagnitudeSq Analytical Gradient (Chain Rule) ---
    inline std::shared_ptr<const CfemExpression> MagnitudeSqExpression::grad() const {
        auto g = m_expr->grad();
        if (!g) {CFEM_ERROR("You requested the gradient, but the object is in wrong state. nullptr."); return nullptr;}
        auto two = std::make_shared<ConstantExpression>(std::vector<CfemReal>{2.0});
        // grad(f^2) = 2 * f * grad(f)
        // return two * m_expr * g;
        return std::make_shared<BinaryExpression>(std::make_shared<BinaryExpression>(two, m_expr, MathOp::Mul), g, MathOp::Mul);
    }

}