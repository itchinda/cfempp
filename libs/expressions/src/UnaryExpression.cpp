#include "MathOperators.h"
#include "ZeroExpression.h"
#include "UnaryExpression.h"
#include "BinaryExpression.h"
#include "ConstantExpression.h"

namespace cfem::sym
{

    /**
     * @brief Internal helper to compute f'(u). 
     * This handles the "outer" part of the chain rule for all supported ops.
     */
    std::shared_ptr<const CfemExpression> UnaryExpression::getDerivativeWrtInner() const {
        auto u = m_expr;
        auto one = std::make_shared<ConstantExpression>(1.0, m_dimIn);
        auto zero = std::make_shared<ConstantExpression>(0.0, m_dimIn);

        switch (m_op) {
            // Basic Math
            case MathOp::Negate: return std::make_shared<ConstantExpression>(-1.0, m_dimIn);
            case MathOp::Abs:    return std::make_shared<UnaryExpression>(u, MathOp::Sign);
            case MathOp::Sqrt:   return std::make_shared<ConstantExpression>(0.5, m_dimIn) / std::make_shared<UnaryExpression>(u, MathOp::Sqrt);
            case MathOp::Cbrt:   {
                auto cbrt_u = std::make_shared<UnaryExpression>(u, MathOp::Cbrt);
                return (std::make_shared<ConstantExpression>(1.0/3.0, m_dimIn) * cbrt_u * cbrt_u);
            }

            // Step Functions (Derivatives are zero almost everywhere)
            case MathOp::Sign:      
            case MathOp::Heaviside: return zero; 

            // Trigonometry
            case MathOp::Sin: return std::make_shared<UnaryExpression>(u, MathOp::Cos);
            case MathOp::Cos: return std::make_shared<ConstantExpression>(-1.0, m_dimIn) * std::make_shared<UnaryExpression>(u, MathOp::Sin);
            case MathOp::Tan: {
                auto sec_u = std::make_shared<UnaryExpression>(u, MathOp::Sec);
                return sec_u * sec_u;
            }
            case MathOp::Cot: {
                auto csc_u = std::make_shared<UnaryExpression>(u, MathOp::Csc);
                return std::make_shared<ConstantExpression>(-1.0, m_dimIn) * csc_u * csc_u;
            }
            case MathOp::Sec: return std::make_shared<UnaryExpression>(u, MathOp::Sec) * std::make_shared<UnaryExpression>(u, MathOp::Tan);
            case MathOp::Csc: return std::make_shared<ConstantExpression>(-1.0, m_dimIn) * std::make_shared<UnaryExpression>(u, MathOp::Csc) * std::make_shared<UnaryExpression>(u, MathOp::Cot);

            // Inverse Trigonometry
            case MathOp::Asin: return one / std::make_shared<UnaryExpression>(one - u*u, MathOp::Sqrt);
            case MathOp::Acos: return std::make_shared<ConstantExpression>(-1.0, m_dimIn) / std::make_shared<UnaryExpression>(one - u*u, MathOp::Sqrt);
            case MathOp::Atan: return one / (one + u*u);
            case MathOp::Acot: return std::make_shared<ConstantExpression>(-1.0, m_dimIn) / (one + u*u);
            case MathOp::Asec: return one / (std::make_shared<UnaryExpression>(u, MathOp::Abs) * std::make_shared<UnaryExpression>(u*u - one, MathOp::Sqrt));
            case MathOp::Acsc: return std::make_shared<ConstantExpression>(-1.0, m_dimIn) / (std::make_shared<UnaryExpression>(u, MathOp::Abs) * std::make_shared<UnaryExpression>(u*u - one, MathOp::Sqrt));

            // Hyperbolic
            case MathOp::Sinh:  return std::make_shared<UnaryExpression>(u, MathOp::Cosh);
            case MathOp::Cosh:  return std::make_shared<UnaryExpression>(u, MathOp::Sinh);
            case MathOp::Tanh:  {
                auto tanh_u = std::make_shared<UnaryExpression>(u, MathOp::Tanh);
                return one - tanh_u * tanh_u;
            }
            case MathOp::Asinh: return one / std::make_shared<UnaryExpression>(u*u + one, MathOp::Sqrt);
            case MathOp::Acosh: return one / std::make_shared<UnaryExpression>(u*u - one, MathOp::Sqrt);
            case MathOp::Atanh: return one / (one - u*u);

            // Exponential & Special Functions
            case MathOp::Exp: return std::make_shared<UnaryExpression>(u, MathOp::Exp);
            case MathOp::Log: return one / u;
            case MathOp::Erf: {
                auto coeff = std::make_shared<ConstantExpression>(2.0 / std::sqrt(M_PI), m_dimIn);
                return coeff * std::make_shared<UnaryExpression>(std::make_shared<ConstantExpression>(-1.0, m_dimIn) * u * u, MathOp::Exp);
            }
            case MathOp::Erfc: {
                auto coeff = std::make_shared<ConstantExpression>(-2.0 / std::sqrt(M_PI), m_dimIn);
                return coeff * std::make_shared<UnaryExpression>(std::make_shared<ConstantExpression>(-1.0, m_dimIn) * u * u, MathOp::Exp);
            }
            
            // Note: Gamma and Zeta derivatives usually require special functions (Digamma). 
            // If you don't have Digamma in your MathOp list, you may return Zero or throw an error.
            case MathOp::Gamma:
            case MathOp::Zeta:
                CFEM_ERROR << "Symbolic derivative for Gamma/Zeta not yet implemented.";
                return zero;

            default:
                CFEM_ERROR << "Unsupported derivative for math operation.";
                return zero;
        }
    }

    std::shared_ptr<const CfemExpression> UnaryExpression::diff(CfemInt var_idx) const {
        auto d_inner = m_expr->diff(var_idx);
        if (!d_inner) return nullptr;

        if (d_inner->isZero()) return std::make_shared<ZeroExpression>(m_dimOut);

        // Get outer derivative ( f'(u) )
        auto d_op = getDerivativeWrtInner();
        if (!d_op) return nullptr;

        // Chain rule
        return d_op * d_inner; 
    }

    // --- grad() Implementation ---
    inline std::shared_ptr<const CfemExpression> UnaryExpression::grad() const {

        // Standard Chain Rule: f'(u) * grad(u)
        auto du = getDerivativeWrtInner();

        if (m_op == MathOp::Transpose) {
            // Gradient of Transpose(A) is Transpose(Grad(A)) in the first two indices.
            // This requires a specialized Tensor operation.
            return transpose(du);
        }
                
        return du * m_expr->grad(); 
    }

    // --- timeDerivative() Implementation ---
    inline std::shared_ptr<const CfemExpression> UnaryExpression::timeDerivative() const {
        if (m_op == MathOp::Transpose) {
            return std::make_shared<UnaryExpression>(m_expr->timeDerivative(), MathOp::Transpose);
        }

        // Standard Chain Rule: f'(u) * u_dot
        auto du = getDerivativeWrtInner();
        return du * m_expr->timeDerivative();
    }

    
}