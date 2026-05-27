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

#include <string>

namespace cfem::sym
{
    inline void BinaryExpression::evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        // auto bookmark = ctx.getBookmark();

        // std::span<CfemReal> spanL = ctx.getBuffer(m_dimL);
        // std::span<CfemReal> spanR = ctx.getBuffer(m_dimR);

        CfemReal stackL[9], stackR[9];
        std::span<CfemReal> spanL, spanR;
        
        bool needsBookmark = (m_dimL > 9 || m_dimR > 9);
        EvalContext::Bookmark bookmark;

        if (needsBookmark) {
            bookmark = ctx.getBookmark();
            spanL = ctx.getBuffer(m_dimL);
            spanR = ctx.getBuffer(m_dimR);
        } else {
            spanL = {stackL, static_cast<size_t>(m_dimL)};
            spanR = {stackR, static_cast<size_t>(m_dimR)};
        }

        m_lhs->evaluate(cellId, xi, spanL, ctx);
        m_rhs->evaluate(cellId, xi, spanR, ctx);

        switch (m_op)
        {
            case MathOp::Add:{
                for (CfemInt i = 0; i < m_dimOut; ++i) out[i] = spanL[i] + spanR[i];
                break;
            }
            case MathOp::Sub:{
                for (CfemInt i = 0; i < m_dimOut; ++i) out[i] = spanL[i] - spanR[i];
                break;
            }
            case MathOp::Mul: {
                if (m_dimL == 1) { // Scalaire * (Scalaire ou Vecteur ou Tenseur)
                    for (CfemInt i = 0; i < m_dimOut; ++i) out[i] = spanL[0] * spanR[i];
                }
                else if (m_dimR == 1) { // (Vecteur ou Tenseur) * Scalaire
                    for (CfemInt i = 0; i < m_dimOut; ++i) out[i] = spanL[i] * spanR[0];
                }
                else if (m_dimL == 9 && m_dimR == 9) { // TENSEUR * TENSEUR (MatMul)
                    for (CfemInt i = 0; i < 3; ++i) {
                        for (CfemInt j = 0; j < 3; ++j) {
                            out[i * 3 + j] = spanL[i * 3 + 0] * spanR[0 * 3 + j] +
                                             spanL[i * 3 + 1] * spanR[1 * 3 + j] +
                                             spanL[i * 3 + 2] * spanR[2 * 3 + j];
                        }
                    }
                }
                else if (m_dimL == 9 && m_dimR == 3) { // Tenseur * Vecteur -> Vecteur
                    out[0] = spanL[0]*spanR[0] + spanL[1]*spanR[1] + spanL[2]*spanR[2];
                    out[1] = spanL[3]*spanR[0] + spanL[4]*spanR[1] + spanL[5]*spanR[2];
                    out[2] = spanL[6]*spanR[0] + spanL[7]*spanR[1] + spanL[8]*spanR[2];
                }
                else { // Par défaut : Multiplication élément par élément (Hadamard)
                    for (CfemInt i = 0; i < m_dimOut; ++i) out[i] = spanL[i] * spanR[i];
                }
                break;
            }
            case MathOp::Div:{
                if (m_dimR == 1) // (Vecteur ou Tenseur) / Scalaire
                    for (CfemInt i = 0; i < m_dimOut; ++i) out[i] = spanL[i] / spanR[0];
                else // Division élément par élément
                    for (CfemInt i = 0; i < m_dimOut; ++i) out[i] = spanL[i] / spanR[i];
                break;
            }
            case MathOp::Dot:{ 
                if (m_dimL == 3 && m_dimR == 3) { // Vecteur . Vecteur -> Scalaire
                    out[0] = spanL[0]*spanR[0] + spanL[1]*spanR[1] + spanL[2]*spanR[2];
                }
                else if (m_dimL == 9 && m_dimR == 3) { // Tenseur . Vecteur -> Vecteur
                    out[0] = spanL[0]*spanR[0] + spanL[1]*spanR[1] + spanL[2]*spanR[2];
                    out[1] = spanL[3]*spanR[0] + spanL[4]*spanR[1] + spanL[5]*spanR[2];
                    out[2] = spanL[6]*spanR[0] + spanL[7]*spanR[1] + spanL[8]*spanR[2];
                }
                else if (m_dimL == 3 && m_dimR == 9) { // Vecteur . Tenseur -> Vecteur
                    out[0] = spanL[0]*spanR[0] + spanL[1]*spanR[3] + spanL[2]*spanR[6];
                    out[1] = spanL[0]*spanR[1] + spanL[1]*spanR[4] + spanL[2]*spanR[7];
                    out[2] = spanL[0]*spanR[2] + spanL[1]*spanR[5] + spanL[2]*spanR[8];
                }
                else if (m_dimL == 9 && m_dimR == 9) { // Tenseur . Tenseur -> Tenseur (MatMul)
                    for (CfemInt i = 0; i < 3; ++i) {
                        for (CfemInt j = 0; j < 3; ++j) {
                            out[i*3 + j] = spanL[i*3 + 0] * spanR[0*3 + j] +
                                        spanL[i*3 + 1] * spanR[1*3 + j] +
                                        spanL[i*3 + 2] * spanR[2*3 + j];
                        }
                    }
                }
                else if (m_dimL == 1) { // Scalaire . (Vecteur/Tenseur)
                    for (CfemInt i = 0; i < m_dimOut; ++i) out[i] = spanL[0] * spanR[i];
                }
                else if (m_dimR == 1) { // (Vecteur/Tenseur) . Scalaire
                    for (CfemInt i = 0; i < m_dimOut; ++i) out[i] = spanL[i] * spanR[0];
                }
                break;
            }
            case MathOp::Inner: { // Double contraction sigma : epsilon -> Scalaire
                CfemReal sum = 0.0;
                for (CfemInt i = 0; i < 9; ++i) sum += spanL[i] * spanR[i];
                out[0] = sum;
                break;
            }
            case MathOp::Cross:{ // Uniquement 3D
                out[0] = spanL[1] * spanR[2] - spanL[2] * spanR[1];
                out[1] = spanL[2] * spanR[0] - spanL[0] * spanR[2];
                out[2] = spanL[0] * spanR[1] - spanL[1] * spanR[0];
                break;
            }
            case MathOp::Outer:
                // spanL = [a0, a1, a2], spanR = [b0, b1, b2]
                // out = [a0*b0, a0*b1, a0*b2, a1*b0, a1*b1, a1*b2, a2*b0, a2*b1, a2*b2]
                for (CfemInt i = 0; i < 3; ++i) {
                    for (CfemInt j = 0; j < 3; ++j) {
                        out[i * 3 + j] = spanL[i] * spanR[j];
                    }
                }
                break;
            case MathOp::Pow:{
                for (CfemInt i = 0; i < m_dimOut; ++i) out[i] = std::pow(spanL[i], spanR[i]);
                break;
            }
            case MathOp::Max:{
                for (CfemInt i = 0; i < m_dimOut; ++i) out[i] = std::max(spanL[i], spanR[i]);
                break;
            }
            case MathOp::Min:{
                for (CfemInt i = 0; i < m_dimOut; ++i) out[i] = std::min(spanL[i], spanR[i]);
                break;
            }
            case MathOp::Atan2:{
                for (CfemInt i = 0; i < m_dimOut; ++i) out[i] = std::atan2(spanL[i], spanR[i]);
                break;
            }
            default:
                break;
        }
        // On libère la mémoire utilisée par spanL et spanR pour les prochains nœuds
        if (needsBookmark) {
            ctx.releaseToBookmark(bookmark);
        }
    }

    inline void BinaryExpression::evaluateBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> out, EvalContext &ctx) const
    {
        // On sauvegarde l'état de l'Arena avant d'allouer les buffers temporaires par lot
        auto bookmark = ctx.getBookmark();

        // Allocation des buffers SoA pour les nQp points : taille = dimension * nQp
        std::span<CfemReal> spanL = ctx.getBuffer(m_dimL * nQp);
        std::span<CfemReal> spanR = ctx.getBuffer(m_dimR * nQp);

        // Évaluation en cascade de l'arbre d'expression
        m_lhs->evaluateBatch(cellId, nQp, spanL, ctx);
        m_rhs->evaluateBatch(cellId, nQp, spanR, ctx);

        // =========================================================================
        // RAPPEL FORMAT SoA: span[comp * nQp + q]
        // Les boucles internes sur `q` assurent une lecture/écriture séquentielle pure
        // =========================================================================

        switch (m_op)
        {
            case MathOp::Add: {
                for (CfemInt i = 0; i < m_dimOut; ++i) {
                    for (CfemInt q = 0; q < nQp; ++q) {
                        out[i * nQp + q] = spanL[i * nQp + q] + spanR[i * nQp + q];
                    }
                }
                break;
            }
            case MathOp::Sub: {
                for (CfemInt i = 0; i < m_dimOut; ++i) {
                    for (CfemInt q = 0; q < nQp; ++q) {
                        out[i * nQp + q] = spanL[i * nQp + q] - spanR[i * nQp + q];
                    }
                }
                break;
            }
            case MathOp::Mul: {
                if (m_dimL == 1) { // Scalaire * (Scalaire/Vecteur/Tenseur)
                    for (CfemInt i = 0; i < m_dimOut; ++i) {
                        for (CfemInt q = 0; q < nQp; ++q) {
                            out[i * nQp + q] = spanL[q] * spanR[i * nQp + q]; // spanL[0*nQp + q] == spanL[q]
                        }
                    }
                }
                else if (m_dimR == 1) { // (Vecteur/Tenseur) * Scalaire
                    for (CfemInt i = 0; i < m_dimOut; ++i) {
                        for (CfemInt q = 0; q < nQp; ++q) {
                            out[i * nQp + q] = spanL[i * nQp + q] * spanR[q];
                        }
                    }
                }
                else if (m_dimL == 9 && m_dimR == 9) { // Tenseur * Tenseur
                    for (CfemInt i = 0; i < 3; ++i) {
                        for (CfemInt j = 0; j < 3; ++j) {
                            CfemInt outIdx = i * 3 + j;
                            CfemInt L0 = i * 3 + 0, L1 = i * 3 + 1, L2 = i * 3 + 2;
                            CfemInt R0 = 0 * 3 + j, R1 = 1 * 3 + j, R2 = 2 * 3 + j;
                            for (CfemInt q = 0; q < nQp; ++q) {
                                out[outIdx * nQp + q] = spanL[L0 * nQp + q] * spanR[R0 * nQp + q] +
                                                        spanL[L1 * nQp + q] * spanR[R1 * nQp + q] +
                                                        spanL[L2 * nQp + q] * spanR[R2 * nQp + q];
                            }
                        }
                    }
                }
                else if (m_dimL == 9 && m_dimR == 3) { // Tenseur * Vecteur
                    for (CfemInt i = 0; i < 3; ++i) {
                        CfemInt L0 = i * 3 + 0, L1 = i * 3 + 1, L2 = i * 3 + 2;
                        for (CfemInt q = 0; q < nQp; ++q) {
                            out[i * nQp + q] = spanL[L0 * nQp + q] * spanR[0 * nQp + q] +
                                               spanL[L1 * nQp + q] * spanR[1 * nQp + q] +
                                               spanL[L2 * nQp + q] * spanR[2 * nQp + q];
                        }
                    }
                }
                else { // Élément par élément (Hadamard)
                    for (CfemInt i = 0; i < m_dimOut; ++i) {
                        for (CfemInt q = 0; q < nQp; ++q) {
                            out[i * nQp + q] = spanL[i * nQp + q] * spanR[i * nQp + q];
                        }
                    }
                }
                break;
            }
            case MathOp::Div: {
                if (m_dimR == 1) { // (Vecteur/Tenseur) / Scalaire
                    for (CfemInt i = 0; i < m_dimOut; ++i) {
                        for (CfemInt q = 0; q < nQp; ++q) {
                            out[i * nQp + q] = spanL[i * nQp + q] / spanR[q];
                        }
                    }
                } else { // Élément par élément
                    for (CfemInt i = 0; i < m_dimOut; ++i) {
                        for (CfemInt q = 0; q < nQp; ++q) {
                            out[i * nQp + q] = spanL[i * nQp + q] / spanR[i * nQp + q];
                        }
                    }
                }
                break;
            }
            case MathOp::Dot: {
                if (m_dimL == 3 && m_dimR == 3) { // Vecteur . Vecteur -> Scalaire
                    for (CfemInt q = 0; q < nQp; ++q) {
                        out[q] = spanL[0 * nQp + q] * spanR[0 * nQp + q] + 
                                 spanL[1 * nQp + q] * spanR[1 * nQp + q] + 
                                 spanL[2 * nQp + q] * spanR[2 * nQp + q];
                    }
                }
                else if (m_dimL == 9 && m_dimR == 3) { // Tenseur . Vecteur -> Vecteur
                    for (CfemInt i = 0; i < 3; ++i) {
                        CfemInt L0 = i * 3 + 0, L1 = i * 3 + 1, L2 = i * 3 + 2;
                        for (CfemInt q = 0; q < nQp; ++q) {
                            out[i * nQp + q] = spanL[L0 * nQp + q] * spanR[0 * nQp + q] +
                                               spanL[L1 * nQp + q] * spanR[1 * nQp + q] +
                                               spanL[L2 * nQp + q] * spanR[2 * nQp + q];
                        }
                    }
                }
                else if (m_dimL == 3 && m_dimR == 9) { // Vecteur . Tenseur -> Vecteur
                    for (CfemInt j = 0; j < 3; ++j) {
                        CfemInt R0 = 0 * 3 + j, R1 = 1 * 3 + j, R2 = 2 * 3 + j;
                        for (CfemInt q = 0; q < nQp; ++q) {
                            out[j * nQp + q] = spanL[0 * nQp + q] * spanR[R0 * nQp + q] +
                                               spanL[1 * nQp + q] * spanR[R1 * nQp + q] +
                                               spanL[2 * nQp + q] * spanR[R2 * nQp + q];
                        }
                    }
                }
                else if (m_dimL == 9 && m_dimR == 9) { // Tenseur . Tenseur -> Tenseur
                    for (CfemInt i = 0; i < 3; ++i) {
                        for (CfemInt j = 0; j < 3; ++j) {
                            CfemInt outIdx = i * 3 + j;
                            CfemInt L0 = i * 3 + 0, L1 = i * 3 + 1, L2 = i * 3 + 2;
                            CfemInt R0 = 0 * 3 + j, R1 = 1 * 3 + j, R2 = 2 * 3 + j;
                            for (CfemInt q = 0; q < nQp; ++q) {
                                out[outIdx * nQp + q] = spanL[L0 * nQp + q] * spanR[R0 * nQp + q] +
                                                        spanL[L1 * nQp + q] * spanR[R1 * nQp + q] +
                                                        spanL[L2 * nQp + q] * spanR[R2 * nQp + q];
                            }
                        }
                    }
                }
                else if (m_dimL == 1) { // Scalaire . (Vecteur/Tenseur)
                    for (CfemInt i = 0; i < m_dimOut; ++i) {
                        for (CfemInt q = 0; q < nQp; ++q) {
                            out[i * nQp + q] = spanL[q] * spanR[i * nQp + q];
                        }
                    }
                }
                else if (m_dimR == 1) { // (Vecteur/Tenseur) . Scalaire
                    for (CfemInt i = 0; i < m_dimOut; ++i) {
                        for (CfemInt q = 0; q < nQp; ++q) {
                            out[i * nQp + q] = spanL[i * nQp + q] * spanR[q];
                        }
                    }
                }
                break;
            }
            case MathOp::Inner: { // Double contraction sigma : epsilon -> Scalaire
                for (CfemInt q = 0; q < nQp; ++q) out[q] = 0.0;
                // Boucle externe sur les composantes pour vectorisation parfaite
                for (CfemInt i = 0; i < 9; ++i) {
                    for (CfemInt q = 0; q < nQp; ++q) {
                        out[q] += spanL[i * nQp + q] * spanR[i * nQp + q];
                    }
                }
                break;
            }
            case MathOp::Cross: { // Produit vectoriel 3D
                for (CfemInt q = 0; q < nQp; ++q) {
                    out[0 * nQp + q] = spanL[1 * nQp + q] * spanR[2 * nQp + q] - spanL[2 * nQp + q] * spanR[1 * nQp + q];
                    out[1 * nQp + q] = spanL[2 * nQp + q] * spanR[0 * nQp + q] - spanL[0 * nQp + q] * spanR[2 * nQp + q];
                    out[2 * nQp + q] = spanL[0 * nQp + q] * spanR[1 * nQp + q] - spanL[1 * nQp + q] * spanR[0 * nQp + q];
                }
                break;
            }
            case MathOp::Outer: { // Produit dyadique (Vecteur X Vecteur -> Tenseur)
                for (CfemInt i = 0; i < 3; ++i) {
                    for (CfemInt j = 0; j < 3; ++j) {
                        CfemInt outIdx = i * 3 + j;
                        for (CfemInt q = 0; q < nQp; ++q) {
                            out[outIdx * nQp + q] = spanL[i * nQp + q] * spanR[j * nQp + q];
                        }
                    }
                }
                break;
            }
            case MathOp::Pow: {
                for (CfemInt i = 0; i < m_dimOut; ++i) {
                    for (CfemInt q = 0; q < nQp; ++q) {
                        out[i * nQp + q] = std::pow(spanL[i * nQp + q], spanR[i * nQp + q]);
                    }
                }
                break;
            }
            case MathOp::Max: {
                for (CfemInt i = 0; i < m_dimOut; ++i) {
                    for (CfemInt q = 0; q < nQp; ++q) {
                        out[i * nQp + q] = std::max(spanL[i * nQp + q], spanR[i * nQp + q]);
                    }
                }
                break;
            }
            case MathOp::Min: {
                for (CfemInt i = 0; i < m_dimOut; ++i) {
                    for (CfemInt q = 0; q < nQp; ++q) {
                        out[i * nQp + q] = std::min(spanL[i * nQp + q], spanR[i * nQp + q]);
                    }
                }
                break;
            }
            case MathOp::Atan2: {
                for (CfemInt i = 0; i < m_dimOut; ++i) {
                    for (CfemInt q = 0; q < nQp; ++q) {
                        out[i * nQp + q] = std::atan2(spanL[i * nQp + q], spanR[i * nQp + q]);
                    }
                }
                break;
            }
            default:
                break;
        }

        // Libération de la mémoire de l'Arena
        ctx.releaseToBookmark(bookmark);
    }

    inline CfemInt getPrecedence(MathOp op)
    {
        switch (op)
        {
            case MathOp::Add:
            case MathOp::Sub:
                return 1;

            case MathOp::Mul:
            case MathOp::Div:
                return 2;

            case MathOp::Dot:
            case MathOp::Inner:
            case MathOp::Cross:
            case MathOp::Pow:
            case MathOp::Max:
            case MathOp::Min:
            case MathOp::Atan2:
                return 3;

            default:
                return 4;
        }
    }

    inline CfemInt BinaryExpression::precedence() const{
        return getPrecedence(m_op);
    }

    inline std::string maybeParen( const std::shared_ptr<const CfemExpression>& expr, CfemInt parentPrec)
    {
        if (expr->precedence() < parentPrec)
            return "(" + expr->toString() + ")";
        else
            return expr->toString();
    }

    inline std::string BinaryExpression::toString() const
    {
        CfemInt prec = precedence();

        std::string L = maybeParen(m_lhs, prec);
        std::string R = maybeParen(m_rhs, prec + 1);

        switch(m_op)
        {
            case MathOp::Add: return L + " + " + R;
            case MathOp::Sub: return L + " - " + R;
            case MathOp::Mul: return L + " * " + R;
            case MathOp::Div: return L + " / " + R;
            case MathOp::Pow: return L + "**" + R;

            default:
                return opToSym(m_op) + "(" +
                    m_lhs->toString() + ", " +
                    m_rhs->toString() + ")";
        }
    }
}