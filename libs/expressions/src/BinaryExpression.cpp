#include "BinaryExpression.h"    
#include "UnaryExpression.h"
#include "ConstantExpression.h"
    
namespace cfem::sym
{    
    BinaryExpression::BinaryExpression(std::shared_ptr<const CfemExpression> lhs,
                                   std::shared_ptr<const CfemExpression> rhs,
                                   MathOp op)
    : m_lhs(std::move(lhs)), m_rhs(std::move(rhs)), m_op(op)
    {
        CFEM_ASSERT(m_lhs && m_rhs, "Operands cannot be null.");
        m_dimL = m_lhs->getDimension();
        m_dimR = m_rhs->getDimension();

        // --- Dimension Logic & Validation ---
        switch (m_op)
        {
        case MathOp::Add:
        case MathOp::Sub:
            // Rigid Physics: Only allow same-dimension operations
            if (m_dimL != m_dimR)
            {
                CFEM_ERROR << "Addition/Subtraction requires identical dimensions (L:" << std::to_string(m_dimL) << " R:" << std::to_string(m_dimR) << ").";
            }
            m_dimOut = m_dimL;
            break;

        case MathOp::Mul:
        case MathOp::Div:
        case MathOp::Pow:
        case MathOp::Max:
        case MathOp::Min:
        case MathOp::Atan2:
            // Broadcasting: Scalar with Vector/Tensor, or element-wise matching
            if (m_dimL != m_dimR && m_dimL != 1 && m_dimR != 1)
            {
                CFEM_ERROR << "Operation requires matching dimensions or a scalar operand.";
            }
            m_dimOut = std::max(m_dimL, m_dimR);
            break;

        case MathOp::Dot:
            // Supports: 
            // Vec(3) . Vec(3) -> Scal(1)
            // Ten(9) . Vec(3) -> Vec(3)
            // Vec(3) . Ten(9) -> Vec(3)
            // Ten(9) . Ten(9) -> Ten(9) (Matrix Multiplication)
            if (m_dimL == 3 && m_dimR == 3)      m_dimOut = 1;
            else if (m_dimL == 9 && m_dimR == 3) m_dimOut = 3;
            else if (m_dimL == 3 && m_dimR == 9) m_dimOut = 3;
            else if (m_dimL == 9 && m_dimR == 9) m_dimOut = 9;
            else {
                CFEM_ERROR << "Dot product expects combinations of Vector(3) and Tensor(9).";
            }
            break;

        case MathOp::Cross:
            if (m_dimL != 3 || m_dimR != 3)
            {
                CFEM_ERROR << "Cross product is only defined for 3D vectors.";
            }
            m_dimOut = 3;
            break;

        case MathOp::Inner:
            // Double contraction (Frobenius inner product). 
            // Dimensions must match to perform component-wise sum.
            if (m_dimL != m_dimR)
            {
                CFEM_ERROR << "Inner product requires matching dimensions.";
            }
            m_dimOut = 1; 
            break;

        default:
            CFEM_ERROR << "Unsupported operation for BinaryExpression.";
        }
    }

    // void BinaryExpression::evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
    // {
    //     auto bookmark = ctx.getBookmark();

    //     std::span<CfemReal> spanL = ctx.getBuffer(m_dimL);
    //     std::span<CfemReal> spanR = ctx.getBuffer(m_dimR);

    //     m_lhs->evaluate(cellId, xi, spanL, ctx);
    //     m_rhs->evaluate(cellId, xi, spanR, ctx);

    //     switch (m_op)
    //     {
    //     case MathOp::Add:
    //         for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[i] + spanR[i];
    //         break;

    //     case MathOp::Sub:
    //         for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[i] - spanR[i];
    //         break;

    //     case MathOp::Mul:
    //         if (m_dimL == 1) // Scalaire * (Vecteur ou Tenseur)
    //             for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[0] * spanR[i];
    //         else if (m_dimR == 1) // (Vecteur ou Tenseur) * Scalaire
    //             for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[i] * spanR[0];
    //         else // Multiplication élément par élément (Hadamard)
    //             for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[i] * spanR[i];
    //         break;

    //     case MathOp::Div:
    //         if (m_dimR == 1) // (Vecteur ou Tenseur) / Scalaire
    //             for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[i] / spanR[0];
    //         else // Division élément par élément
    //             for (int i = 0; i < m_dimOut; ++i) out[i] = spanL[i] / spanR[i];
    //         break;

    //     case MathOp::Dot: 
    //         if (m_dimL == 3 && m_dimR == 3) { // Vecteur . Vecteur -> Scalaire
    //             out[0] = spanL[0]*spanR[0] + spanL[1]*spanR[1] + spanL[2]*spanR[2];
    //         }
    //         else if (m_dimL == 9 && m_dimR == 3) { // Tenseur . Vecteur -> Vecteur
    //             out[0] = spanL[0]*spanR[0] + spanL[1]*spanR[1] + spanL[2]*spanR[2];
    //             out[1] = spanL[3]*spanR[0] + spanL[4]*spanR[1] + spanL[5]*spanR[2];
    //             out[2] = spanL[6]*spanR[0] + spanL[7]*spanR[1] + spanL[8]*spanR[2];
    //         }
    //         else if (m_dimL == 3 && m_dimR == 9) { // Vecteur . Tenseur -> Vecteur
    //             out[0] = spanL[0]*spanR[0] + spanL[1]*spanR[3] + spanL[2]*spanR[6];
    //             out[1] = spanL[0]*spanR[1] + spanL[1]*spanR[4] + spanL[2]*spanR[7];
    //             out[2] = spanL[0]*spanR[2] + spanL[1]*spanR[5] + spanL[2]*spanR[8];
    //         }
    //         else if (m_dimL == 9 && m_dimR == 9) { // Tenseur . Tenseur -> Tenseur (MatMul)
    //             for (int i = 0; i < 3; ++i) {
    //                 for (int j = 0; j < 3; ++j) {
    //                     out[i*3 + j] = spanL[i*3 + 0] * spanR[0*3 + j] +
    //                                 spanL[i*3 + 1] * spanR[1*3 + j] +
    //                                 spanL[i*3 + 2] * spanR[2*3 + j];
    //                 }
    //             }
    //         }
    //         break;

    //     case MathOp::Inner: { // Double contraction sigma : epsilon -> Scalaire
    //         CfemReal sum = 0.0;
    //         for (int i = 0; i < 9; ++i) sum += spanL[i] * spanR[i];
    //         out[0] = sum;
    //         break;
    //     }

    //     case MathOp::Cross: // Uniquement 3D
    //         out[0] = spanL[1] * spanR[2] - spanL[2] * spanR[1];
    //         out[1] = spanL[2] * spanR[0] - spanL[0] * spanR[2];
    //         out[2] = spanL[0] * spanR[1] - spanL[1] * spanR[0];
    //         break;

    //     case MathOp::Pow:
    //         for (int i = 0; i < m_dimOut; ++i) out[i] = std::pow(spanL[i], spanR[i]);
    //         break;

    //     case MathOp::Max:
    //         for (int i = 0; i < m_dimOut; ++i) out[i] = std::max(spanL[i], spanR[i]);
    //         break;

    //     case MathOp::Min:
    //         for (int i = 0; i < m_dimOut; ++i) out[i] = std::min(spanL[i], spanR[i]);
    //         break;

    //     case MathOp::Atan2:
    //         for (int i = 0; i < m_dimOut; ++i) out[i] = std::atan2(spanL[i], spanR[i]);
    //         break;

    //     default:
    //         break;
    //     }

    //     // On libère la mémoire utilisée par spanL et spanR pour les prochains nœuds
    //     ctx.releaseToBookmark(bookmark);
    // }

    std::shared_ptr<const CfemExpression> BinaryExpression::grad() const {
        auto L = m_lhs;
        auto R = m_rhs;
        auto gL = L->grad();
        auto gR = R->grad();

        // Protection : si une branche n'est pas dérivable, l'opération globale ne l'est pas.
        if (!gL || !gR) {
            CFEM_ERROR << "One of the operands in BinaryExpression is not differentiable.";
            return nullptr;
        }

        switch (m_op) {
            // --- Linéarité : (f ± g)' = f' ± g' ---
            case MathOp::Add:
            case MathOp::Sub:
                return std::make_shared<BinaryExpression>(gL, gR, m_op);

            // --- Règle du Produit (Leibniz) : (f*g)' = f'*g + f*g' ---
            // S'applique à Mul, Dot (produit scalaire), Inner (contraction), Cross (produit vectoriel)
            case MathOp::Mul:
            case MathOp::Dot:
            case MathOp::Inner:
            case MathOp::Cross: {
                auto term1 = std::make_shared<BinaryExpression>(gL, R, m_op);
                auto term2 = std::make_shared<BinaryExpression>(L, gR, m_op);
                return std::make_shared<BinaryExpression>(term1, term2, MathOp::Add);
            }

            // --- Règle du Quotient : (f/g)' = (f'*g - f*g') / g² ---
            case MathOp::Div: {
                auto t1 = std::make_shared<BinaryExpression>(gL, R, MathOp::Mul);
                auto t2 = std::make_shared<BinaryExpression>(L, gR, MathOp::Mul);
                auto num = std::make_shared<BinaryExpression>(t1, t2, MathOp::Sub);
                auto den = std::make_shared<BinaryExpression>(R, R, MathOp::Mul);
                return std::make_shared<BinaryExpression>(num, den, MathOp::Div);
            }

            // --- Puissance Générale : (f^g)' = f^g * (g' * ln(f) + g * f' / f) ---
            case MathOp::Pow: {
                // Cas n * f^(n-1) * f' (Terme de base)
                auto un = std::make_shared<ConstantExpression>(1.0);
                auto gm1 = std::make_shared<BinaryExpression>(R, un, MathOp::Sub);
                auto f_pow_gm1 = std::make_shared<BinaryExpression>(L, gm1, MathOp::Pow);
                auto coeff = std::make_shared<BinaryExpression>(R, f_pow_gm1, MathOp::Mul);
                auto term1 = std::make_shared<BinaryExpression>(coeff, gL, MathOp::Mul);

                // Cas ln(f) * f^g * g' (Si l'exposant R est aussi une variable)
                auto ln_f = std::make_shared<UnaryExpression>(L, MathOp::Log);
                auto f_pow_g = std::make_shared<BinaryExpression>(L, R, MathOp::Pow);
                auto ln_f_f_g = std::make_shared<BinaryExpression>(ln_f, f_pow_g, MathOp::Mul);
                auto term2 = std::make_shared<BinaryExpression>(ln_f_f_g, gR, MathOp::Mul);

                return std::make_shared<BinaryExpression>(term1, term2, MathOp::Add);
            }

            // --- Atan2 : grad(atan2(y, x)) = (x*gy - y*gx) / (x² + y²) ---
            case MathOp::Atan2: {
                auto x_gy = std::make_shared<BinaryExpression>(R, gL, MathOp::Mul);
                auto y_gx = std::make_shared<BinaryExpression>(L, gR, MathOp::Mul);
                auto num = std::make_shared<BinaryExpression>(x_gy, y_gx, MathOp::Sub);
                
                auto x2 = std::make_shared<BinaryExpression>(R, R, MathOp::Mul);
                auto y2 = std::make_shared<BinaryExpression>(L, L, MathOp::Mul);
                auto den = std::make_shared<BinaryExpression>(x2, y2, MathOp::Add);
                
                return std::make_shared<BinaryExpression>(num, den, MathOp::Div);
            }

            // --- Opérateurs Logiques / Comparaison ---
            // La dérivée est nulle presque partout (on considère ces fonctions comme des constantes par morceaux)

            // --- Max / Min : Non-lisse ---
            case MathOp::Max:
            case MathOp::Min:
                // Mathématiquement : (f > g) ? f' : g'
                // Nécessite un nœud de type "ConditionalExpression" que nous n'avons pas encore détaillé.
                CFEM_ERROR << "Gradient of Max/Min is non-smooth and requires conditional node support.";
                return nullptr;

            default:
                CFEM_ERROR << "Binary operation not supported for differentiation.";
                return nullptr;
        }
    }

    std::shared_ptr<const CfemExpression> BinaryExpression::timeDerivative() const {
        auto L = m_lhs;
        auto R = m_rhs;
        auto dL = L->timeDerivative();
        auto dR = R->timeDerivative();

        // Protection: if a branch isn't differentiable over time, the global op isn't either.
        if (!dL || !dR) {
            CFEM_ERROR << "One of the operands in BinaryExpression is not differentiable wrt time.";
            return nullptr;
        }

        switch (m_op) {
            // --- Linearity: d(f ± g)/dt = df/dt ± dg/dt ---
            case MathOp::Add:
            case MathOp::Sub:
                return std::make_shared<BinaryExpression>(dL, dR, m_op);

            // --- Product Rule: d(f*g)/dt = (df/dt)*g + f*(dg/dt) ---
            // Applies to all bilinear multiplications
            case MathOp::Mul:
            case MathOp::Dot:
            case MathOp::Inner:
            case MathOp::Cross: {
                auto term1 = std::make_shared<BinaryExpression>(dL, R, m_op);
                auto term2 = std::make_shared<BinaryExpression>(L, dR, m_op);
                return std::make_shared<BinaryExpression>(term1, term2, MathOp::Add);
            }

            // --- Quotient Rule: d(f/g)/dt = ((df/dt)*g - f*(dg/dt)) / g^2 ---
            case MathOp::Div: {
                auto t1 = std::make_shared<BinaryExpression>(dL, R, MathOp::Mul);
                auto t2 = std::make_shared<BinaryExpression>(L, dR, MathOp::Mul);
                auto num = std::make_shared<BinaryExpression>(t1, t2, MathOp::Sub);
                auto den = std::make_shared<BinaryExpression>(R, R, MathOp::Mul);
                return std::make_shared<BinaryExpression>(num, den, MathOp::Div);
            }

            // --- General Power Rule: d(f^g)/dt = f^g * ( (dg/dt)*ln(f) + g*(df/dt)/f ) ---
            case MathOp::Pow: {
                auto f_pow_g = std::make_shared<BinaryExpression>(L, R, MathOp::Pow);
                
                // Term A: g * f^(g-1) * (df/dt)
                auto un = std::make_shared<ConstantExpression>(1.0, 1);
                auto gm1 = std::make_shared<BinaryExpression>(R, un, MathOp::Sub);
                auto f_pow_gm1 = std::make_shared<BinaryExpression>(L, gm1, MathOp::Pow);
                auto term1 = std::make_shared<BinaryExpression>(
                    std::make_shared<BinaryExpression>(R, f_pow_gm1, MathOp::Mul), 
                    dL, MathOp::Mul);

                // Term B: ln(f) * f^g * (dg/dt)
                auto ln_f = std::make_shared<UnaryExpression>(L, MathOp::Log);
                auto term2 = std::make_shared<BinaryExpression>(
                    std::make_shared<BinaryExpression>(ln_f, f_pow_g, MathOp::Mul), 
                    dR, MathOp::Mul);

                return std::make_shared<BinaryExpression>(term1, term2, MathOp::Add);
            }

            // --- Atan2: d(atan2(y, x))/dt = (x*(dy/dt) - y*(dx/dt)) / (x^2 + y^2) ---
            case MathOp::Atan2: {
                auto x_dy = std::make_shared<BinaryExpression>(R, dL, MathOp::Mul);
                auto y_dx = std::make_shared<BinaryExpression>(L, dR, MathOp::Mul);
                auto num = std::make_shared<BinaryExpression>(x_dy, y_dx, MathOp::Sub);
                
                auto x2 = std::make_shared<BinaryExpression>(R, R, MathOp::Mul);
                auto y2 = std::make_shared<BinaryExpression>(L, L, MathOp::Mul);
                auto den = std::make_shared<BinaryExpression>(x2, y2, MathOp::Add);
                
                return std::make_shared<BinaryExpression>(num, den, MathOp::Div);
            }

            // --- Non-Smooth Operations ---
            case MathOp::Max:
            case MathOp::Min:
                CFEM_ERROR << "Time derivative of Max/Min requires conditional node support.";
                return nullptr;

            default:
                CFEM_ERROR << "Binary operation not supported for time differentiation.";
                return nullptr;
        }
    }

}