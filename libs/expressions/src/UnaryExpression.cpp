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

#include "BinaryExpression.h"
#include "ConstantExpression.h"
#include "MathOperators.h"
#include "ZeroExpression.h"
#include "UnaryExpression.h"
#include "VectorExpression.h"

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
    std::shared_ptr<const CfemExpression> UnaryExpression::createGradient() const {

        // Standard Chain Rule: f'(u) * grad(u)
        auto du = getDerivativeWrtInner();

        if (m_op == MathOp::Transpose) {
            return transpose(du);
        }
                
        return du * m_expr->grad(); 
    }

    std::shared_ptr<const CfemExpression> UnaryExpression::createTimeDerivative() const {
        if (m_op == MathOp::Transpose) {
            return std::make_shared<UnaryExpression>(m_expr->timeDerivative(), MathOp::Transpose);
        }

        // Standard Chain Rule: f'(u) * u_dot
        auto du = getDerivativeWrtInner();
        return du * m_expr->timeDerivative();
    }

    std::shared_ptr<const CfemExpression> UnaryExpression::div() const {
        using namespace utils;

        if (m_dimOut == 1) {
            CFEM_ERROR << "Divergence undefined for scalar result of: " << opToSym(m_op);
        }

        // Case A: Negation (Linear)
        if (m_op == MathOp::Negate) return -m_expr->div();

        // Case B: Transpose (Indices: Row sum becomes Column sum)
        if (m_op == MathOp::Transpose && m_dimIn == 9) {
            auto dx = m_expr->diff(0); 
            auto dy = m_expr->diff(1); 
            auto dz = m_expr->diff(2); 

            return std::make_shared<const VectorExpression>(
                extractComponent(dx, 0) + extractComponent(dy, 3) + extractComponent(dz, 6), // Div Row 1
                extractComponent(dx, 1) + extractComponent(dy, 4) + extractComponent(dz, 7), // Div Row 2
                extractComponent(dx, 2) + extractComponent(dy, 5) + extractComponent(dz, 8)  // Div Row 3
            );
        }

        // Case C: General Pointwise Chain Rule (Sin, Exp, etc.)
        auto f_prime = getDerivativeWrtInner();
        auto dx = m_expr->diff(0);
        auto dy = m_expr->diff(1);
        auto dz = m_expr->diff(2);

        // Lambda helper to safely multiply and drop zeros
        auto smart_mul = [](const auto& a, const auto& b) -> std::shared_ptr<const CfemExpression> {
            if (a->isZero() || b->isZero()) {
                return std::make_shared<const ConstantExpression>(0.0, 1); // Returns scalar 0
            }
            if (a->isConstant()) {
                auto v = std::static_pointer_cast<const ConstantExpression>(a)->getValues();
                bool isOnes = std::all_of(v.begin(), v.end(), [](CfemReal x) { return NumericPolicy::areEqual(x, 1.0); });
                if (isOnes)
                    return b; // 1 * b = b
            }
            if (b->isConstant()) {
                auto v = std::static_pointer_cast<const ConstantExpression>(b)->getValues();
                bool isOnes = std::all_of(v.begin(), v.end(), [](CfemReal x) { return NumericPolicy::areEqual(x, 1.0); });
                if (isOnes)
                    return a; // 1 * b = b
            }
            return a * b; // Fallback to normal multiplication tree node
        };

        // Lambda helper to safely add and drop zeros
        auto smart_add = [](const auto& a, const auto& b) -> std::shared_ptr<const CfemExpression> {
            if (a->isZero()) return b;
            if (b->isZero()) return a;
            return a + b;
        };

        if (m_dimIn == 3) {
            // Result is a Scalar
            auto term_x = smart_mul(extractComponent(f_prime, 0), extractComponent(dx, 0));
            auto term_y = smart_mul(extractComponent(f_prime, 1), extractComponent(dy, 1));
            auto term_z = smart_mul(extractComponent(f_prime, 2), extractComponent(dz, 2));

            auto sum_xy = smart_add(term_x, term_y);
            return smart_add(sum_xy, term_z);
        } 
        
        if (m_dimIn == 9) {
            // Result is a Vector (Row-wise Divergence)
            auto row_sum = [&](int i, int j, int k) {
                auto tx = smart_mul(extractComponent(f_prime, i), extractComponent(dx, i));
                auto ty = smart_mul(extractComponent(f_prime, j), extractComponent(dy, j));
                auto tz = smart_mul(extractComponent(f_prime, k), extractComponent(dz, k));
                
                return smart_add(smart_add(tx, ty), tz);
            };

            auto row1 = row_sum(0, 1, 2);
            auto row2 = row_sum(3, 4, 5);
            auto row3 = row_sum(6, 7, 8);

            return std::make_shared<const VectorExpression>(row1, row2, row3);
        }

        CFEM_ERROR << "Divergence not implemented for operator: " << opToSym(m_op);
        return nullptr;
    }

    
}