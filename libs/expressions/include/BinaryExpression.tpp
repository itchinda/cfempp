#include <string>

namespace cfem::sym
{
    inline void BinaryExpression::evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    {
        auto bookmark = ctx.getBookmark();

        std::span<CfemReal> spanL = ctx.getBuffer(m_dimL);
        std::span<CfemReal> spanR = ctx.getBuffer(m_dimR);

        m_lhs->evaluate(cellId, xi, spanL, ctx);
        m_rhs->evaluate(cellId, xi, spanR, ctx);

        switch (m_op)
        {
            case MathOp::Add:{
                for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[i] + spanR[i];
                break;
            }
            case MathOp::Sub:{
                for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[i] - spanR[i];
                break;
            }
            case MathOp::Mul:{
                if (m_dimL == 1) // Scalaire * (Vecteur ou Tenseur)
                    for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[0] * spanR[i];
                else if (m_dimR == 1) // (Vecteur ou Tenseur) * Scalaire
                    for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[i] * spanR[0];
                else // Multiplication élément par élément (Hadamard)
                    for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[i] * spanR[i];
                break;
            }
            case MathOp::Div:{
                if (m_dimR == 1) // (Vecteur ou Tenseur) / Scalaire
                    for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[i] / spanR[0];
                else // Division élément par élément
                    for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[i] / spanR[i];
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
                    for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                            out[i*3 + j] = spanL[i*3 + 0] * spanR[0*3 + j] +
                                        spanL[i*3 + 1] * spanR[1*3 + j] +
                                        spanL[i*3 + 2] * spanR[2*3 + j];
                        }
                    }
                }
                break;
            }
            case MathOp::Inner: { // Double contraction sigma : epsilon -> Scalaire
                CfemReal sum = 0.0;
                for (int i = 0; i < 9; ++i) sum += spanL[i] * spanR[i];
                out[0] = sum;
                break;
            }
            case MathOp::Cross:{ // Uniquement 3D
                out[0] = spanL[1] * spanR[2] - spanL[2] * spanR[1];
                out[1] = spanL[2] * spanR[0] - spanL[0] * spanR[2];
                out[2] = spanL[0] * spanR[1] - spanL[1] * spanR[0];
                break;
            }
            case MathOp::Pow:{
                for (int i = 0; i < m_dimOut; ++i) out[i] = std::pow(spanL[i], spanR[i]);
                break;
            }
            case MathOp::Max:{
                for (int i = 0; i < m_dimOut; ++i) out[i] = std::max(spanL[i], spanR[i]);
                break;
            }
            case MathOp::Min:{
                for (int i = 0; i < m_dimOut; ++i) out[i] = std::min(spanL[i], spanR[i]);
                break;
            }
            case MathOp::Atan2:{
                for (int i = 0; i < m_dimOut; ++i) out[i] = std::atan2(spanL[i], spanR[i]);
                break;
            }
            default:
                break;
        }
        // On libère la mémoire utilisée par spanL et spanR pour les prochains nœuds
        ctx.releaseToBookmark(bookmark);
    }

    inline int getPrecedence(MathOp op)
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