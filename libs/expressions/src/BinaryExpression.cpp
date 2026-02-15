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
                CFEM_ERROR("Addition/Subtraction requires identical dimensions.");
            }
            m_dimOut = m_dimL;
            break;

        case MathOp::Mul:
        case MathOp::Div:
            // Broadcasting: Scalar * Vector, Vector * Scalar, or Vector * Vector (Hadamard)
            if (m_dimL != m_dimR && m_dimL != 1 && m_dimR != 1)
            {
                CFEM_ERROR("Multiplication/Division requires matching dimensions or a scalar operand.");
            }
            m_dimOut = std::max(m_dimL, m_dimR);
            break;

        case MathOp::Dot:
            if (m_dimL != m_dimR || m_dimL == 1)
            {
                CFEM_ERROR("Dot product requires matching vector dimensions (dim > 1).");
            }
            m_dimOut = 1; // Reduction to scalar
            break;

        case MathOp::Cross:
            if (m_dimL != 3 || m_dimR != 3)
            {
                CFEM_ERROR("Cross product is only defined for 3D vectors.");
            }
            m_dimOut = 3;
            break;

        case MathOp::Inner:
            m_dimOut = 1; break;

        default : CFEM_ERROR("Unsupported operation for BinaryExpression.");
        }
        
    }

    std::shared_ptr<const CfemExpression> BinaryExpression::grad() const {
        auto L = m_lhs;
        auto R = m_rhs;
        auto gL = L->grad();
        auto gR = R->grad();

        // Protection : si une branche n'est pas dérivable, l'opération globale ne l'est pas.
        if (!gL || !gR) {
            CFEM_ERROR("One of the operands in BinaryExpression is not differentiable.");
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
                CFEM_ERROR("Gradient of Max/Min is non-smooth and requires conditional node support.");
                return nullptr;

            default:
                CFEM_ERROR("Binary operation not supported for differentiation.");
                return nullptr;
        }
    }

}