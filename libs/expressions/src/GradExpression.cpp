#include "GradExpression.h"
#include "BinaryExpression.h"
#include "BinaryExpression.h"
#include "MagnitudeSqExpression.h"
#include "FEField.h"


namespace cfem::sym
{

    void GradExpression::evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const 
    {
        const auto type = m_expr->getType();

        // 1. Chemin Analytique (Symbolique)
        // On ne le fait que si ce n'est pas un champ/trial/test pour éviter les boucles infinies
        if (type != ExpressionType::FEField && type != ExpressionType::Trial && type != ExpressionType::Test) 
        {
            if (auto analytical = m_expr->grad()) {
                if (analytical->getType() != ExpressionType::Gradient) {
                    analytical->evaluate(cellId, xi, out, ctx);
                    return;
                }
            }
        }

        // 2. Cas des opérations linéaires (Add/Sub) : grad(a + b) = grad(a) + grad(b)
        if (type == ExpressionType::Binary)
        {
            auto bin = static_cast<const BinaryExpression*>(m_expr.get());
            MathOp op = bin->getOp();
            
            if (op == MathOp::Add || op == MathOp::Sub) {
                const size_t sz = out.size();
                CfemReal stackL[9], stackR[9];
                std::span<CfemReal> spanL, spanR;
                
                bool useScratch = (sz > 9);
                EvalContext::Bookmark bookmark;

                if (useScratch) {
                    bookmark = ctx.getBookmark();
                    spanL = ctx.getBuffer(sz);
                    spanR = ctx.getBuffer(sz);
                } else {
                    spanL = {stackL, sz};
                    spanR = {stackR, sz};
                }

                // Appel récursif via de nouveaux GradExpression (temporaires légers)
                GradExpression(bin->getLHS()).evaluate(cellId, xi, spanL, ctx);
                GradExpression(bin->getRHS()).evaluate(cellId, xi, spanR, ctx);

                if (op == MathOp::Add) {
                    for (size_t i = 0; i < sz; ++i) out[i] = spanL[i] + spanR[i];
                } else {
                    for (size_t i = 0; i < sz; ++i) out[i] = spanL[i] - spanR[i];
                }

                if (useScratch) ctx.releaseToBookmark(bookmark);
                return;
            }
        }

        // 3. Règle de la chaîne pour MagnitudeSq : grad(f^2) = 2 * f * grad(f)
        if (type == ExpressionType::MagnitudeSq)
        {
            auto magSq = static_cast<const MagnitudeSqExpression*>(m_expr.get());
            auto inner = magSq->getInner(); 
            const int dimIn = inner->getDimension();
            const size_t szGrad = static_cast<size_t>(dimIn) * 3;

            CfemReal stackVal[3], stackGrad[9]; // Supporte jusqu'à un vecteur 3D
            std::span<CfemReal> valIn, gradIn;
            
            bool useScratch = (szGrad > 9);
            EvalContext::Bookmark bookmark;

            if (useScratch) {
                bookmark = ctx.getBookmark();
                valIn = ctx.getBuffer(dimIn);
                gradIn = ctx.getBuffer(szGrad);
            } else {
                valIn = {stackVal, static_cast<size_t>(dimIn)};
                gradIn = {stackGrad, szGrad};
            }

            inner->evaluate(cellId, xi, valIn, ctx);
            GradExpression(inner).evaluate(cellId, xi, gradIn, ctx);
            
            std::fill(out.begin(), out.end(), 0.0);
            for (int d = 0; d < dimIn; ++d) {
                const CfemReal twoVal = 2.0 * valIn[d];
                out[0] += twoVal * gradIn[d * 3 + 0];
                out[1] += twoVal * gradIn[d * 3 + 1];
                out[2] += twoVal * gradIn[d * 3 + 2];
            }
            
            if (useScratch) ctx.releaseToBookmark(bookmark);
            return;
        }

        // 4. Cas FEField (Le plus fréquent en calcul de résidus)
        if (type == ExpressionType::FEField) {
            static_cast<const FEField *>(m_expr.get())->evaluateGradientAt(cellId, xi, out, ctx);
            return;
        }

        // 5. Cas Trial / Test (Fonctions de forme)
        if (type == ExpressionType::Trial || type == ExpressionType::Test) {
            const bool isTrial = (type == ExpressionType::Trial);
            const auto* shape = isTrial ? ctx.trialShape : ctx.testShape;
            const auto* map   = isTrial ? ctx.trialMapping : ctx.testMapping;
            const int dof     = isTrial ? ctx.trialDof : ctx.testDof;

            const int numNodes = shape->values.size();
            const int comp = dof / numNodes;
            const int node = dof % numNodes;

            const auto& g = map->physicalShapeGradients[node];

            std::fill(out.begin(), out.end(), 0.0);
            out[comp*3 + 0] = g.x;
            out[comp*3 + 1] = g.y;
            out[comp*3 + 2] = g.z;
            return;
        }

        CFEM_ERROR << "GradExpression: Failed to differentiate expression tree.";
    }

}