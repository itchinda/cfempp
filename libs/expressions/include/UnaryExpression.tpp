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

#include "UnaryExpression.h"
#include "ConstantExpression.h"

namespace cfem::sym
{
    inline UnaryExpression::UnaryExpression(std::shared_ptr<const CfemExpression> expr, MathOp op)
            : m_expr(std::move(expr)), m_op(op) 
    {
        m_dimIn = m_expr->getNumComponents();

        // --- Dimension Logic & Validation ---
        if (m_op == MathOp::Transpose || m_op == MathOp::Trace ) {
            if (m_dimIn != 9) {
                CFEM_THROW("Transpose/Trace are only defined for dim 9 tensors.");
            }
        }

        // All other operations in our Unary list are component-wise (Hadamard).
        // This means a Vector in gives a Vector out, Tensor in gives Tensor out.
        m_dimOut = m_dimIn;
    }

    inline void UnaryExpression::evaluate(CfemInt cellIdx, const Vector3D& point, 
                  std::span<CfemReal> out, EvalContext& ctx) const 
    {

        CfemReal stack[9];
        std::span<CfemReal> in_buf;
        
        bool needsBookmark = (m_dimIn > 9);
        EvalContext::Bookmark bookmark;

        if (needsBookmark) {
            bookmark = ctx.getBookmark();
            in_buf = ctx.getBuffer(m_dimIn);
        } else {
            in_buf = {stack, static_cast<size_t>(m_dimIn)};
        }
        
        // Evaluate the child node
        m_expr->evaluate(cellIdx, point, in_buf, ctx);

        switch (m_op) {
            // Matrix Operations
            case MathOp::Transpose:
                out[0] = in_buf[0]; out[1] = in_buf[3]; out[2] = in_buf[6];
                out[3] = in_buf[1]; out[4] = in_buf[4]; out[5] = in_buf[7];
                out[6] = in_buf[2]; out[7] = in_buf[5]; out[8] = in_buf[8];
                break;
            case MathOp::Trace:
                out[0] = in_buf[0] + in_buf[4] + in_buf[8];
                break;
            case MathOp::Det: {
                // Determinant using Laplace expansion (Sarrus Rule)
                out[0] = in_buf[0] * (in_buf[4] * in_buf[8] - in_buf[5] * in_buf[7]) -
                         in_buf[1] * (in_buf[3] * in_buf[8] - in_buf[5] * in_buf[6]) +
                         in_buf[2] * (in_buf[3] * in_buf[7] - in_buf[4] * in_buf[6]);
                break;
            }

            case MathOp::Inv: {
                CfemReal det = in_buf[0] * (in_buf[4] * in_buf[8] - in_buf[5] * in_buf[7]) -
                               in_buf[1] * (in_buf[3] * in_buf[8] - in_buf[5] * in_buf[6]) +
                               in_buf[2] * (in_buf[3] * in_buf[7] - in_buf[4] * in_buf[6]);

                if (NumericPolicy::isZero(std::abs(det))) {
                    for(int i=0; i<9; ++i) out[i] = 0.0; 
                    break;
                }

                CfemReal invDet = 1.0 / det;

                // 3. Compute Adjugate matrix (Cofactors transposed)
                out[0] = (in_buf[4] * in_buf[8] - in_buf[5] * in_buf[7]) * invDet;
                out[1] = (in_buf[2] * in_buf[7] - in_buf[1] * in_buf[8]) * invDet;
                out[2] = (in_buf[1] * in_buf[5] - in_buf[2] * in_buf[4]) * invDet;

                out[3] = (in_buf[5] * in_buf[6] - in_buf[3] * in_buf[8]) * invDet;
                out[4] = (in_buf[0] * in_buf[8] - in_buf[2] * in_buf[6]) * invDet;
                out[5] = (in_buf[2] * in_buf[3] - in_buf[0] * in_buf[5]) * invDet;

                out[6] = (in_buf[3] * in_buf[7] - in_buf[4] * in_buf[6]) * invDet;
                out[7] = (in_buf[1] * in_buf[6] - in_buf[0] * in_buf[7]) * invDet;
                out[8] = (in_buf[0] * in_buf[4] - in_buf[1] * in_buf[3]) * invDet;
                break;
            }
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
                for(int i=0; i<m_dimOut; ++i) out[i] = std::riemann_zeta(in_buf[i]);
                break;
            default:
                CFEM_THROW("Unsupported math operation.");
                break;
        if (needsBookmark) {
            ctx.releaseToBookmark(bookmark); 
        }
        }
    }

    inline CfemInt UnaryExpression::precedence() const{
        return 5;
    }

    inline std::string UnaryExpression::toString() const {
        switch (m_op){
            case MathOp::Negate: return "(-(" + m_expr->toString() + "))";
                break;
            default:
                return opToSym(m_op) + "(" + m_expr->toString() + ")";
                break;
        }
        return opToSym(m_op) + "(" + m_expr->toString() + ")";
    }

    
}