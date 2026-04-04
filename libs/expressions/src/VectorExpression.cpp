#include "VectorExpression.h"
#include "TensorExpression.h"
#include "symengine/ScalarSymbolicExpression.h"
#include "symengine/VectorSymbolicExpression.h"
#include "symengine/TensorSymbolicExpression.h"
#include "BinaryExpression.h"
#include "UnaryExpression.h"
#include "ScaledExpression.h"
#include "LambdaExpression.h"
#include "GradExpression.h"

#include "MathOperators.h"

namespace cfem::sym
{
    std::shared_ptr<const CfemExpression> VectorExpression::grad() const
    {
        // 1. Get the vector gradients of the x, y, and z scalar components.
        // Each g_i is a Vector expression (dimension 3).
        auto g_x = m_components[0]->grad(); 
        auto g_y = m_components[1]->grad(); 
        auto g_z = m_components[2]->grad(); 

        // 2. Recursive helper to extract the i-th scalar component safely.
        // We use std::function only for the structural traversal; 
        // the actual "leaf" evaluation is wrapped in an optimized template.
        std::function<std::shared_ptr<const CfemExpression>(std::shared_ptr<const CfemExpression>, int)> extractComponent = 
            [&extractComponent](std::shared_ptr<const CfemExpression> expr, int idx) -> std::shared_ptr<const CfemExpression> 
        {
            if (!expr) return nullptr;

            // --- Container Passthrough ---
            if (auto v = std::dynamic_pointer_cast<const VectorSymbolicExpression>(expr)) return v->getComponent(idx);
            if (auto v = std::dynamic_pointer_cast<const VectorExpression>(expr))         return v->getComponent(idx);
            if (auto t = std::dynamic_pointer_cast<const TensorSymbolicExpression>(expr)) return t->getComponent(idx);
            
            // --- Structural Recursion (Preserve the tree) ---
            if (auto v = std::dynamic_pointer_cast<const ScaledExpression>(expr)) {
                return std::make_shared<const ScaledExpression>(v->getCoeff(), extractComponent(v->getBaseShared(), idx));
            }
            if (auto v = std::dynamic_pointer_cast<const UnaryExpression>(expr)) {
                return std::make_shared<const UnaryExpression>(extractComponent(v->getInner(), idx), v->getOp());
            }
            if (auto b = std::dynamic_pointer_cast<const BinaryExpression>(expr)) {
                MathOp op = b->getOp();
                if (op == MathOp::Add || op == MathOp::Sub) {
                    return std::make_shared<const BinaryExpression>(
                        extractComponent(b->getLHS(), idx), 
                        extractComponent(b->getRHS(), idx), op);
                }
                if (op == MathOp::Mul) {
                    if (b->getLHS()->getDimension() == 1) 
                        return std::make_shared<const BinaryExpression>(b->getLHS(), extractComponent(b->getRHS(), idx), op);
                    if (b->getRHS()->getDimension() == 1) 
                        return std::make_shared<const BinaryExpression>(extractComponent(b->getLHS(), idx), b->getRHS(), op);
                }
            }

            // --- Opaque Fallback (Lambdas, Fields, GradExpressions) ---
            // We capture by value to ensure thread-safety and avoid lifetime issues.
            auto eval_func = [expr, idx](const Vector3D& pos, CfemReal t) -> CfemReal {
                CfemReal temp[9] = {0.0}; // Stack-based, no heap allocation
                EvalContext ctx;
                ctx.physPos = pos;
                ctx.currentTime = t;
                expr->evaluate(-1, pos, temp, ctx); 
                return temp[idx];
            };

            // Handle Time Derivative propagation
            std::shared_ptr<const CfemExpression> comp_time_deriv = nullptr;
            if (auto td = expr->timeDerivative()) {
                if (!td->isZero()) {
                    comp_time_deriv = extractComponent(td, idx);
                }
            }

            // Wrap the fallback in our optimized template class
            return LambdaExpression<decltype(eval_func)>::fromScalar(
                std::move(eval_func), nullptr, comp_time_deriv
            );
        };

        // 3. Extract components to form the Jacobian rows
        // Row 0: Grad of u_x
        auto C00 = extractComponent(g_x, 0); auto C01 = extractComponent(g_x, 1); auto C02 = extractComponent(g_x, 2);
        // Row 1: Grad of u_y
        auto C10 = extractComponent(g_y, 0); auto C11 = extractComponent(g_y, 1); auto C12 = extractComponent(g_y, 2);
        // Row 2: Grad of u_z
        auto C20 = extractComponent(g_z, 0); auto C21 = extractComponent(g_z, 1); auto C22 = extractComponent(g_z, 2);

        // 4. Build the final TensorExpression (3x3 Jacobian)
        return std::make_shared<TensorExpression>(
            C00, C01, C02,
            C10, C11, C12,
            C20, C21, C22
        );
    }

    std::shared_ptr<const CfemExpression> VectorExpression::div() const {
        auto dux_dx = m_components[0]->diff(0); 
        auto duy_dy = m_components[1]->diff(1); 
        auto duz_dz = m_components[2]->diff(2); 

        return dux_dx + duy_dy + duz_dz;
    }
}