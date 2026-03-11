#pragma once

#include "UnaryExpression.h"
#include "ConstantExpression.h"

namespace cfem::sym
{
    inline UnaryExpression::UnaryExpression(std::shared_ptr<const CfemExpression> expr, MathOp op)
            : m_expr(std::move(expr)), m_op(op) 
    {
        m_dimIn = m_expr->getDimension();

        // --- Dimension Logic & Validation ---
        if (m_op == MathOp::Transpose) {
            if (m_dimIn != 9) {
                CFEM_ERROR("Transpose is only defined for dim 9 tensors.");
            }
        }

        // All other operations in our Unary list are component-wise (Hadamard).
        // This means a Vector in gives a Vector out, Tensor in gives Tensor out.
        m_dimOut = m_dimIn;
    }

    inline void UnaryExpression::evaluate(CfemInt cellIdx, const Vector3D& point, 
                  std::span<CfemReal> out, EvalContext& ctx) const 
    {
        // Get pre-allocated thread-local scratchpad to prevent heap allocations
        auto bookmark = ctx.getBookmark(); 

        int dim = m_expr->getDimension();
        auto in_buf = ctx.getBuffer(dim);

        
        // Evaluate the child node
        m_expr->evaluate(cellIdx, point, in_buf, ctx);

        switch (m_op) {
            // Matrix Operations
            case MathOp::Transpose:
                out[0] = in_buf[0]; out[1] = in_buf[3]; out[2] = in_buf[6];
                out[3] = in_buf[1]; out[4] = in_buf[4]; out[5] = in_buf[7];
                out[6] = in_buf[2]; out[7] = in_buf[5]; out[8] = in_buf[8];
                break;
            // Basic Math
            case MathOp::Negate: for(int i=0; i<m_dimOut; ++i) out[i] = -in_buf[i]; break;
            case MathOp::Abs:    for(int i=0; i<m_dimOut; ++i) out[i] = std::abs(in_buf[i]); break;
            case MathOp::Sqrt:   for(int i=0; i<m_dimOut; ++i) out[i] = (in_buf[i] > 0) ? std::sqrt(in_buf[i]) : 0.0; break;
            case MathOp::Cbrt:   for(int i=0; i<m_dimOut; ++i) out[i] = std::cbrt(in_buf[i]); break;
            // Logical / Step Functions
            case MathOp::Sign:      for(int i=0; i<m_dimOut; ++i) out[i] = (in_buf[i] > 0) ? 1.0 : ((in_buf[i] < 0) ? -1.0 : 0.0); break;
            case MathOp::Heaviside: for(int i=0; i<m_dimOut; ++i) out[i] = (in_buf[i] >= 0) ? 1.0 : 0.0; break;
            // Trigonometry
            case MathOp::Sin: for(int i=0; i<m_dimOut; ++i) out[i] = std::sin(in_buf[i]); break;
            case MathOp::Cos: for(int i=0; i<m_dimOut; ++i) out[i] = std::cos(in_buf[i]); break;
            case MathOp::Tan: for(int i=0; i<m_dimOut; ++i) out[i] = std::tan(in_buf[i]); break;
            case MathOp::Cot: for(int i=0; i<m_dimOut; ++i) out[i] = 1.0 / std::tan(in_buf[i]); break; // Relies on IEEE 754 infinity for div by 0
            case MathOp::Sec: for(int i=0; i<m_dimOut; ++i) out[i] = 1.0 / std::cos(in_buf[i]); break;
            case MathOp::Csc: for(int i=0; i<m_dimOut; ++i) out[i] = 1.0 / std::sin(in_buf[i]); break;
            // Inverse Trigonometry
            case MathOp::Asin: for(int i=0; i<m_dimOut; ++i) out[i] = std::asin(in_buf[i]); break;
            case MathOp::Acos: for(int i=0; i<m_dimOut; ++i) out[i] = std::acos(in_buf[i]); break;
            case MathOp::Atan: for(int i=0; i<m_dimOut; ++i) out[i] = std::atan(in_buf[i]); break;
            case MathOp::Acot: for(int i=0; i<m_dimOut; ++i) out[i] = M_PI_2 - std::atan(in_buf[i]); break;
            case MathOp::Asec: for(int i=0; i<m_dimOut; ++i) out[i] = std::acos(1.0 / in_buf[i]); break;
            case MathOp::Acsc: for(int i=0; i<m_dimOut; ++i) out[i] = std::asin(1.0 / in_buf[i]); break;
            // Hyperbolic
            case MathOp::Sinh:  for(int i=0; i<m_dimOut; ++i) out[i] = std::sinh(in_buf[i]); break;
            case MathOp::Cosh:  for(int i=0; i<m_dimOut; ++i) out[i] = std::cosh(in_buf[i]); break;
            case MathOp::Tanh:  for(int i=0; i<m_dimOut; ++i) out[i] = std::tanh(in_buf[i]); break;
            case MathOp::Asinh: for(int i=0; i<m_dimOut; ++i) out[i] = std::asinh(in_buf[i]); break;
            case MathOp::Acosh: for(int i=0; i<m_dimOut; ++i) out[i] = std::acosh(in_buf[i]); break;
            case MathOp::Atanh: for(int i=0; i<m_dimOut; ++i) out[i] = std::atanh(in_buf[i]); break;
            // Exponential & Special Functions
            case MathOp::Exp:   for(int i=0; i<m_dimOut; ++i) out[i] = std::exp(in_buf[i]); break;
            case MathOp::Log:   for(int i=0; i<m_dimOut; ++i) out[i] = (in_buf[i] > 0) ? std::log(in_buf[i]) : -1e20; break;
            case MathOp::Gamma: for(int i=0; i<m_dimOut; ++i) out[i] = std::tgamma(in_buf[i]); break;
            case MathOp::Erf:   for(int i=0; i<m_dimOut; ++i) out[i] = std::erf(in_buf[i]); break;
            case MathOp::Erfc:  for(int i=0; i<m_dimOut; ++i) out[i] = std::erfc(in_buf[i]); break;
            case MathOp::Zeta:  
                #if __cplusplus >= 201703L
                    for(int i=0; i<m_dimOut; ++i) out[i] = std::riemann_zeta(in_buf[i]);
                #else
                    CFEM_ERROR("Zeta function requires C++17 or higher.");
                #endif
                break;

            default:
                CFEM_ERROR("Unsupported math operation.");
                break;

        ctx.releaseToBookmark(bookmark); 
        }
    }
}