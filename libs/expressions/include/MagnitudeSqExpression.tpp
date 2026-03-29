#pragma once

#include "BinaryExpression.h"
#include "ConstantExpression.h"

namespace cfem::sym
{

    // inline void MagnitudeSqExpression::evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    // {
    //     auto bookmark = ctx.getBookmark(); 

    //     int dim = m_expr->getDimension();
        
    //     // On utilise le scratchpad au lieu de allouer sur la pile/tas
    //     std::span<CfemReal> spanBuf = ctx.getBuffer(dim);

    //     m_expr->evaluate(cellId, xi, spanBuf, ctx);

    //     CfemReal sumSq = 0.0;
    //     for (int i = 0; i < dim; ++i){
    //         sumSq += spanBuf[i] * spanBuf[i];
    //     }
    //     out[0] = sumSq;

    //     ctx.releaseToBookmark(bookmark); 
    // }

    inline void MagnitudeSqExpression::evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        // 1. Utilisation d'un buffer sur la pile pour les cas courants (dim <= 9)
        // 9 couvre les gradients de vecteurs 3D (3x3).
        CfemReal stackBuf[9];
        std::span<CfemReal> targetSpan;
        bool usesScratchpad = false;
        EvalContext::Bookmark bookmark;

        if (m_innerDim <= 9) {
            targetSpan = {stackBuf, static_cast<size_t>(m_innerDim)};
        } else {
            // Fallback pour les cas exotiques (très hautes dimensions)
            bookmark = ctx.getBookmark();
            targetSpan = ctx.getBuffer(m_innerDim);
            usesScratchpad = true;
        }

        // 2. Évaluation de l'expression interne
        m_expr->evaluate(cellId, xi, targetSpan, ctx);

        // 3. Calcul de la somme des carrés avec un switch pour aider le compilateur
        CfemReal sumSq = 0.0;
        switch (m_innerDim) {
            case 1: 
                sumSq = targetSpan[0] * targetSpan[0]; 
                break;
            case 3: // Cas ultra-fréquent : Gradient de scalaire ou Vecteur 3D
                sumSq = targetSpan[0]*targetSpan[0] + targetSpan[1]*targetSpan[1] + targetSpan[2]*targetSpan[2];
                break;
            case 9: // Cas fréquent : Gradient de vecteur 3D
                sumSq = targetSpan[0]*targetSpan[0] + targetSpan[1]*targetSpan[1] + targetSpan[2]*targetSpan[2] +
                        targetSpan[3]*targetSpan[3] + targetSpan[4]*targetSpan[4] + targetSpan[5]*targetSpan[5] +
                        targetSpan[6]*targetSpan[6] + targetSpan[7]*targetSpan[7] + targetSpan[8]*targetSpan[8];
                //for (int i = 0; i < 9; ++i) sumSq += targetSpan[i] * targetSpan[i];
                break;
            default:
                for (int i = 0; i < m_innerDim; ++i) sumSq += targetSpan[i] * targetSpan[i];
                break;
        }

        out[0] = sumSq;

        // 4. Nettoyage uniquement si nécessaire
        if (usesScratchpad) {
            ctx.releaseToBookmark(bookmark);
        }
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