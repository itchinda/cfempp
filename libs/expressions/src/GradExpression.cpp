#include "GradExpression.h"
#include "BinaryExpression.h"
#include "BinaryExpression.h"
#include "MagnitudeSqExpression.h"
#include "FEField.h"


namespace cfem::sym
{
    void GradExpression::evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const 
    {
        // 1. Analytical path (Symbolic/Constants) -> ZERO ALLOCATION, PARFAIT.
        // if (auto analytical = m_expr->grad()) {
        //     analytical->evaluate(cellId, xi, out, ctx);
        //     return;
        // }

        if (m_expr->getType() != ExpressionType::FEField && 
                    m_expr->getType() != ExpressionType::Trial && 
                    m_expr->getType() != ExpressionType::Test) 
        {
            if (auto analytical = m_expr->grad()) {
                // Safety check: ensure we aren't just recursing on ourselves
                if (analytical->getType() != ExpressionType::Gradient) {
                    analytical->evaluate(cellId, xi, out, ctx);
                    return;
                }
            }
        }

        // 2. Recursive path for Binary (Add/Sub)
        if (m_expr->getType() == ExpressionType::Binary){
            auto bin = dynamic_cast<const BinaryExpression*>(m_expr.get());
            MathOp op = bin->getOp();
            if (op == MathOp::Add || op == MathOp::Sub) {
                // --- DEBUT DU BLOC TEMPORAIRE ---
                auto bookmark = ctx.getBookmark();

                // On demande exactement la taille requise (out.size()), fini la limite de 9 !
                std::span<CfemReal> spanL = ctx.getBuffer(out.size());
                std::span<CfemReal> spanR = ctx.getBuffer(out.size());

                GradExpression(bin->getLHS()).evaluate(cellId, xi, spanL, ctx);
                GradExpression(bin->getRHS()).evaluate(cellId, xi, spanR, ctx);

                for (size_t i = 0; i < out.size(); ++i) {
                    out[i] = (op == MathOp::Add) ? (spanL[i] + spanR[i]) : (spanL[i] - spanR[i]);
                }
                
                // --- FIN DU BLOC : on rend la mémoire au scratchpad ---
                ctx.releaseToBookmark(bookmark);
                return;
            }
        }

        // 3. Chain Rule for MagnitudeSq: grad(f^2) = 2 * f * grad(f)
        if (m_expr->getType() == ExpressionType::MagnitudeSq){
            auto magSq = dynamic_cast<const MagnitudeSqExpression*>(m_expr.get());
            auto inner = magSq->getInner(); 
            int dimIn = inner->getDimension();
            
            // --- DEBUT DU BLOC TEMPORAIRE ---
            auto bookmark = ctx.getBookmark();

            // On demande dynamiquement l'espace pour la valeur et son gradient
            std::span<CfemReal> valIn = ctx.getBuffer(dimIn);
            std::span<CfemReal> gradIn = ctx.getBuffer(dimIn * 3);

            inner->evaluate(cellId, xi, valIn, ctx);
            GradExpression(inner).evaluate(cellId, xi, gradIn, ctx);
            
            std::fill(out.begin(), out.end(), 0.0);
            for (int d = 0; d < dimIn; ++d) {
                for (int i = 0; i < 3; ++i) {
                    out[i] += 2.0 * valIn[d] * gradIn[d * 3 + i];
                }
            }
            
            // --- FIN DU BLOC ---
            ctx.releaseToBookmark(bookmark);
            return;
        }

        if (m_expr->getType() == ExpressionType::FEField) {
            const FEField *fe = dynamic_cast<const FEField *>(m_expr.get());
            fe->evaluateGradientAt(cellId, xi, out, ctx);
            return;
        }

        if (m_expr->getType() == ExpressionType::Trial) {

            const CfemInt dim = m_expr->getDimension();
            const CfemInt numNodes = ctx.trialShape->values.size();

            const CfemInt comp = ctx.trialDof / numNodes;
            const CfemInt node = ctx.trialDof % numNodes;

            const auto& g = ctx.trialMapping->physicalShapeGradients[node];

            // zero output (size = dim * 3 for gradient tensor)
            std::fill(out.begin(), out.end(), CfemReal(0));

            // ∂u_comp / ∂x_j
            out[comp*3 + 0] = g.x;
            out[comp*3 + 1] = g.y;
            out[comp*3 + 2] = g.z;

            return;
        }

        if (m_expr->getType() == ExpressionType::Test) {

            const CfemInt dim = m_expr->getDimension();
            const CfemInt numNodes = ctx.testShape->values.size();

            const CfemInt comp = ctx.testDof / numNodes;
            const CfemInt node = ctx.testDof % numNodes;

            const auto& g = ctx.testMapping->physicalShapeGradients[node];

            // zero output (size = dim * 3 for gradient tensor)
            std::fill(out.begin(), out.end(), CfemReal(0));

            // ∂u_comp / ∂x_j
            out[comp*3 + 0] = g.x;
            out[comp*3 + 1] = g.y;
            out[comp*3 + 2] = g.z;

            return;
        }

        CFEM_ERROR << "GradExpression: Failed to differentiate expression tree.";
    }

}