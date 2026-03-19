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

     /**
     * @brief Computes 2 * u * grad(u).
     * Since u is dim N and grad(u) is dim N*3, the dot product results in dim 3.
     */
    inline std::shared_ptr<const CfemExpression> MagnitudeSqExpression::grad() const 
    {
        auto two = std::make_shared<ConstantExpression>(2.0);
        return std::make_shared<BinaryExpression>(
            two, 
            std::make_shared<BinaryExpression>(m_expr, m_expr->grad(), MathOp::Dot),
            MathOp::Mul
        );
    }

    /**
     * @brief Computes 2 * u * u_dot.
     * Scalar result (dim 1).
     */
    inline std::shared_ptr<const CfemExpression> MagnitudeSqExpression::timeDerivative() const 
    {
        auto two = std::make_shared<ConstantExpression>(2.0);
        
        return std::make_shared<BinaryExpression>(
            two,
            std::make_shared<BinaryExpression>(m_expr, m_expr->timeDerivative(), MathOp::Dot),
            MathOp::Mul
        );
    }

}