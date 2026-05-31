//************************************************************************
// --- CFEM++ Library
// ---
// --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.
// --- ALL RIGHTS RESERVED
// ---
// --- This software is protected by international copyright laws.
// --- Use, distribution, or modification of this software in any form,
// --- source or binary, for personal, academic, or non-commercial 
// --- purposes is permitted free of charge, provided that this 
// --- copyright notice and this permission notice appear in all 
// --- copies and supporting documentation.
// ---
// --- The authors and contributors provide this code "as is" without 
// --- any express or implied warranty. In no event shall they be 
// --- held liable for any damages arising from the use of this software.
//************************************************************************

#pragma once

#include "BinaryExpression.h"
#include "ConstantExpression.h"

namespace cfem::sym
{

    inline void MagnitudeSqExpression::evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        CfemReal stackBuf[9];
        std::span<CfemReal> targetSpan;
        bool usesScratchpad = false;
        EvalContext::Bookmark bookmark;

        if (m_innerDim <= 9) {
            targetSpan = {stackBuf, static_cast<size_t>(m_innerDim)};
        } else {
            bookmark = ctx.getBookmark();
            targetSpan = ctx.getBuffer(m_innerDim);
            usesScratchpad = true;
        }

        m_expr->evaluate(cellId, xi, targetSpan, ctx);

        CfemReal sumSq = 0.0;
        switch (m_innerDim) {
            case 1: 
                sumSq = targetSpan[0] * targetSpan[0]; 
                break;
            case 3:
                sumSq = targetSpan[0]*targetSpan[0] + targetSpan[1]*targetSpan[1] + targetSpan[2]*targetSpan[2];
                break;
            case 9:
                sumSq = targetSpan[0]*targetSpan[0] + targetSpan[1]*targetSpan[1] + targetSpan[2]*targetSpan[2] +
                        targetSpan[3]*targetSpan[3] + targetSpan[4]*targetSpan[4] + targetSpan[5]*targetSpan[5] +
                        targetSpan[6]*targetSpan[6] + targetSpan[7]*targetSpan[7] + targetSpan[8]*targetSpan[8];
                break;
            default:
                for (int i = 0; i < m_innerDim; ++i) sumSq += targetSpan[i] * targetSpan[i];
                break;
        }

        out[0] = sumSq;

        if (usesScratchpad) {
            ctx.releaseToBookmark(bookmark);
        }
    }

    inline void MagnitudeSqExpression::evaluateBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const
    {
        // La magnitude au carré produit un scalaire. Le buffer de sortie doit faire au moins nQp.
        CFEM_ASSERT(outValues.size() >= static_cast<size_t>(nQp), "Batch evaluation output buffer too small for scalar result");

        // Initialisation de la somme à zéro pour tous les points de quadrature
        for (CfemInt q = 0; q < nQp; ++q) {
            outValues[q] = 0.0;
        }

        // Sauvegarde de l'Arena et allocation du buffer SoA pour l'expression interne
        auto bookmark = ctx.getBookmark();
        std::span<CfemReal> in_buf = ctx.getBuffer(m_innerDim * nQp);

        m_expr->evaluateBatch(cellId, nQp, in_buf, ctx);

        // =========================================================================
        // FORMAT SoA: in_buf[composante * nQp + q]
        // MagnitudeSq = Somme des carrés des composantes.
        // La boucle externe sur les composantes (i) assure la vectorisation SIMD
        // de la boucle interne (q).
        // =========================================================================
        for (CfemInt i = 0; i < m_innerDim; ++i) {
            CfemInt i_x_nQp = i * nQp;
            for (CfemInt q = 0; q < nQp; ++q) {
                CfemReal val = in_buf[i_x_nQp + q];
                outValues[q] += val * val;
            }
        }

        // Libération de la mémoire
        ctx.releaseToBookmark(bookmark);
    }

     /**
     * @brief Computes 2 * u * grad(u).
     * Since u is dim N and grad(u) is dim N*3, the dot product results in dim 3.
     */
    inline std::shared_ptr<const CfemExpression> MagnitudeSqExpression::createGradient() const 
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
    inline std::shared_ptr<const CfemExpression> MagnitudeSqExpression::createTimeDerivative() const 
    {
        auto two = std::make_shared<ConstantExpression>(2.0);
        
        return std::make_shared<BinaryExpression>(
            two,
            std::make_shared<BinaryExpression>(m_expr, m_expr->timeDerivative(), MathOp::Dot),
            MathOp::Mul
        );
    }

}