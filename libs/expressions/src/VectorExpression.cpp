#include "VectorExpression.h"
#include "TensorExpression.h"
#include "symengine/ScalarSymbolicExpression.h"
#include "symengine/VectorSymbolicExpression.h"
#include "BinaryExpression.h"
#include "UnaryExpression.h"
#include "ScaledExpression.h"
#include "LambdaExpression.h"

namespace cfem::sym
{
    std::shared_ptr<const CfemExpression> VectorExpression::grad() const
    {
        // 1. Get the vector gradients of the x, y, and z scalar components.
        auto g_x = m_components[0]->grad(); 
        auto g_y = m_components[1]->grad(); 
        auto g_z = m_components[2]->grad(); 

        // 2. Recursive Lambda to extract the i-th scalar component safely from ANY vector expression
        std::function<std::shared_ptr<const CfemExpression>(std::shared_ptr<const CfemExpression>, int)> extractComponent = 
            [&extractComponent](std::shared_ptr<const CfemExpression> expr, int idx) -> std::shared_ptr<const CfemExpression> 
        {
            if (!expr) return nullptr;

            // --- Transparent Container Types ---
            if (auto v = std::dynamic_pointer_cast<const VectorSymbolicExpression>(expr)) {
                return v->getComponent(idx);
            }
            else if (auto v = std::dynamic_pointer_cast<const VectorExpression>(expr)) {
                return v->getComponent(idx);
            }
            
            // --- Unary and Scaled Modifiers ---
            else if (auto v = std::dynamic_pointer_cast<const ScaledExpression>(expr)) {
                return std::make_shared<ScaledExpression>(v->getCoeff(), extractComponent(v->getBaseShared(), idx));
            }
            else if (auto v = std::dynamic_pointer_cast<const UnaryExpression>(expr)) {
                return std::make_shared<UnaryExpression>(extractComponent(v->getInner(), idx), v->getOp());
            }

            // --- Binary Operations (The Missing Piece) ---
            else if (auto b = std::dynamic_pointer_cast<const BinaryExpression>(expr)) {
                MathOp op = b->getOp();
                
                // Vector +/- Vector
                if (op == MathOp::Add || op == MathOp::Sub) {
                    return std::make_shared<BinaryExpression>(
                        extractComponent(b->getLHS(), idx), 
                        extractComponent(b->getRHS(), idx), op);
                }
                // Scalar * Vector  OR  Vector * Scalar
                else if (op == MathOp::Mul) {
                    if (b->getLHS()->getDimension() == 1) {
                        return std::make_shared<BinaryExpression>(b->getLHS(), extractComponent(b->getRHS(), idx), op);
                    } else if (b->getRHS()->getDimension() == 1) {
                        return std::make_shared<BinaryExpression>(extractComponent(b->getLHS(), idx), b->getRHS(), op);
                    }
                }
                // Vector / Scalar
                else if (op == MathOp::Div) {
                    if (b->getRHS()->getDimension() == 1) {
                        return std::make_shared<BinaryExpression>(extractComponent(b->getLHS(), idx), b->getRHS(), op);
                    }
                }
                // Cross Product: A x B
                else if (op == MathOp::Cross) {
                    // C_idx = A_i1 * B_i2 - A_i2 * B_i1
                    int i1 = (idx + 1) % 3;
                    int i2 = (idx + 2) % 3;
                    auto term1 = std::make_shared<BinaryExpression>(
                        extractComponent(b->getLHS(), i1), extractComponent(b->getRHS(), i2), MathOp::Mul);
                    auto term2 = std::make_shared<BinaryExpression>(
                        extractComponent(b->getLHS(), i2), extractComponent(b->getRHS(), i1), MathOp::Mul);
                    return std::make_shared<BinaryExpression>(term1, term2, MathOp::Sub);
                }
                // If it's a Pow, Max, Min, etc. applied to vectors, it will fall through to the lambda wrapper below.
            }

            // --- Opaque Black Box Fallback (Lambdas, Variables, Unhandled Ops) ---
            // We use the static factory 'fromScalar' from your LambdaExpression class!
            auto eval_func = [expr, idx](const Vector3D& pos, CfemReal t) -> CfemReal {
                CfemReal temp[3] = {0.0, 0.0, 0.0};
                EvalContext ctx;
                ctx.physPos = pos;
                ctx.currentTime = t;
                expr->evaluate(-1, pos, temp, ctx); 
                return temp[idx];
            };

            std::shared_ptr<const CfemExpression> comp_time_deriv = nullptr;
            try {
                if (auto td = expr->timeDerivative()) {
                    comp_time_deriv = extractComponent(td, idx);
                }
            } catch (...) { /* Safe to ignore if no time derivative exists */ }

            return LambdaExpression::fromScalar(eval_func, nullptr, comp_time_deriv);
        };

        // 3. Extract components row by row
        auto Cxx = extractComponent(g_x, 0); auto Cxy = extractComponent(g_x, 1); auto Cxz = extractComponent(g_x, 2);
        auto Cyx = extractComponent(g_y, 0); auto Cyy = extractComponent(g_y, 1); auto Cyz = extractComponent(g_y, 2);
        auto Czx = extractComponent(g_z, 0); auto Czy = extractComponent(g_z, 1); auto Czz = extractComponent(g_z, 2);

        // 4. Build the final Tensor
        return std::make_shared<TensorExpression>(Cxx, Cxy, Cxz,
                                                  Cyx, Cyy, Cyz,
                                                  Czx, Czy, Czz);
    }
}