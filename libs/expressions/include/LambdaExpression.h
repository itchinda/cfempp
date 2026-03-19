#pragma once 

#include "cfemutils.h"
#include "Vector3D.h"
#include "TensorO2.h"
#include "CfemExpression.h"
#include "EvalContext.h"
#include <vector>
#include <span>
#include <functional>
#include <memory>

namespace cfem::sym
{

    /** *************************************************************************
     * @brief Wrapper to use arbitrary C++ lambdas/functions as symbolic expressions.
     * **************************************************************************** */
    class LambdaExpression : public CfemExpression
    {
        std::function<void(const Vector3D&, CfemReal, std::span<CfemReal>)> m_func;       //! Evaluation function
        CfemInt                                                             m_dim;        //! Dimension of the output (1, 3, 9, etc.)
        std::shared_ptr<const CfemExpression>                               m_gradExpr;   //! Spatial gradient expression
        std::shared_ptr<const CfemExpression>                               m_timeDeriv;  //! Time derivative expression

    public:
        // Constructor for purely spatial lambdas
        LambdaExpression(
                    std::function<void(const Vector3D&, std::span<CfemReal>)> f,
                    CfemInt                                                   dim,
                    std::shared_ptr<const CfemExpression>                     gradExpr = nullptr);

        // Constructor for spatiotemporal lambdas
        LambdaExpression(
                    std::function<void(const Vector3D&, CfemReal, std::span<CfemReal>)> f,
                    CfemInt                                                             dim,
                    std::shared_ptr<const CfemExpression>                               gradExpr = nullptr,
                    std::shared_ptr<const CfemExpression>                               timeDeriv = nullptr);

        CfemInt getDimension() const override { return m_dim; }
        bool isPointwise() const override { return true; }

        void evaluate(CfemInt, const Vector3D&, std::span<CfemReal> out, EvalContext& ctx) const override {
            CFEM_ASSERT(out.size() >= m_dim, "LambdaExpression: output span too small.");
            m_func(ctx.physPos, ctx.currentTime, out);
        }

        void setTimeDerivative(std::shared_ptr<const CfemExpression> expr) {
            if (expr && expr->getDimension() != m_dim) {
                CFEM_ERROR << "LambdaExpression::setTimeDerivative: Dimension mismatch.";
            }
            m_timeDeriv = std::move(expr);
        }

        /**
         * @brief Defines or updates the expression used for the time derivative.
         * Useful for assigning numerical schemes like (U - U_old) / dt.
         */
        void setGradient(std::shared_ptr<const CfemExpression> expr) {
            if (expr && expr->getDimension() != m_dim * 3) {
                CFEM_ERROR << "LambdaExpression::setGradient: Dimension mismatch.";
            }
            m_gradExpr = std::move(expr);
        }

        std::shared_ptr<const CfemExpression> grad() const override {
            if (!m_gradExpr) {
                CFEM_ERROR << "Analytical gradient was not provided for this LambdaExpression. Cannot compute spatial derivatives for lambda.";
            }
            return m_gradExpr; 
        }
        
        // MISSING OVERRIDE ADDED HERE:
        std::shared_ptr<const CfemExpression> timeDerivative() const override { 
            if (!m_timeDeriv){
                CFEM_ERROR << "Analytical time derivative was not provided for this LambdaExpression. Cannot compute time derivatives for lambda.";
            }
            return m_timeDeriv; 
        }
        
        std::string toString() const override { return "LambdaExpression_dim_" + std::to_string(m_dim); }
        ExpressionType getType() const override { return ExpressionType::Lambda; }
        
        // --- Static Factories for purely spatial lambdas ---

        static std::shared_ptr<const LambdaExpression> fromScalar(
                                                    const std::function<CfemReal(const Vector3D&)>& f,
                                                    std::shared_ptr<const CfemExpression>           gradient = nullptr);

        static std::shared_ptr<const LambdaExpression> fromVector(
                                                    const std::function<Vector3D(const Vector3D&)>& f,
                                                    std::shared_ptr<const CfemExpression>           gradient = nullptr);
                
        static std::shared_ptr<const LambdaExpression> fromTensor(
                                                    const std::function<TensorO2(const Vector3D&)>& f,
                                                    std::shared_ptr<const CfemExpression>           gradient = nullptr);
        
        // --- Static Factories for spatiotemporal lambdas ---

        static std::shared_ptr<const LambdaExpression> fromScalar(
                                                    const std::function<CfemReal(const Vector3D&, CfemReal)>& f,
                                                    std::shared_ptr<const CfemExpression>                     gradient = nullptr,
                                                    std::shared_ptr<const CfemExpression>                     timeDeriv = nullptr);

        static std::shared_ptr<const LambdaExpression> fromVector(
                                                    const std::function<Vector3D(const Vector3D&, CfemReal)>& f,
                                                    std::shared_ptr<const CfemExpression>                     gradient = nullptr,
                                                    std::shared_ptr<const CfemExpression>                     timeDeriv = nullptr);
                
        static std::shared_ptr<const LambdaExpression> fromTensor(
                                                    const std::function<TensorO2(const Vector3D&, CfemReal)>& f,
                                                    std::shared_ptr<const CfemExpression>                     gradient = nullptr,
                                                    std::shared_ptr<const CfemExpression>                     timeDeriv = nullptr);
    };

}