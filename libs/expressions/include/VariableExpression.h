#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
#include "TensorO2.h"
#include "CfemExpression.h"
#include "ZeroExpression.h"
#include "EvalContext.h"
#include <span>
#include <memory>

namespace cfem::sym
{

    class VariableExpression : public CfemExpression {
    private:
        std::string m_name;
        const CfemReal* m_dataPtr;
        CfemInt m_dim;
        std::shared_ptr<const CfemExpression> m_gradExpr; 
        std::shared_ptr<const CfemExpression> m_timeDerivExpr;

    public:
        VariableExpression(const std::string& name, const CfemReal* ptr, 
                        CfemInt dim, 
                        std::shared_ptr<const CfemExpression> gradExpr = nullptr,
                        std::shared_ptr<const CfemExpression> timeDerivExpr = nullptr) 
            : m_name(name), m_dataPtr(ptr), m_dim(dim), m_gradExpr(gradExpr), m_timeDerivExpr(timeDerivExpr)
        {
            NamingPolicy::validate(name);
        }

        void evaluate(CfemInt cellId, const Vector3D& xi, std::span<CfemReal> out, EvalContext& ctx) const override {
            for (CfemInt i = 0; i < m_dim; ++i) out[i] = m_dataPtr[i];
        }

        std::shared_ptr<const CfemExpression> grad() const override {
            if (!m_gradExpr) {
                CFEM_ERROR << "Gradient was not provided. Cannot compute spatial derivatives for VariableExpression.";
            }
            return m_gradExpr;
        }

        /**
         * @brief Defines or updates the expression used for the time derivative.
         * Useful for assigning numerical schemes like (U - U_old) / dt.
         */
        void setTimeDerivative(std::shared_ptr<const CfemExpression> expr) {
            if (expr && expr->getDimension() != m_dim) {
                CFEM_ERROR << "VariableExpression::setTimeDerivative: Dimension mismatch.";
            }
            m_timeDerivExpr = std::move(expr);
        }

        /**
         * @brief Defines or updates the expression used for the time derivative.
         * Useful for assigning numerical schemes like (U - U_old) / dt.
         */
        void setGradient(std::shared_ptr<const CfemExpression> expr) {
            if (expr && expr->getDimension() != m_dim * 3) {
                CFEM_ERROR << "VariableExpression::setGradient: Dimension mismatch.";
            }
            m_gradExpr = std::move(expr);
        }

        /**
         * @brief Returns the time derivative expression.
         * @throw CFEM_ERROR if no time derivative or numerical scheme has been assigned.
         */
        std::shared_ptr<const CfemExpression> timeDerivative() const override {
            if (!m_timeDerivExpr) {
                CFEM_ERROR << "VariableExpression '" << m_name << 
                    "' has no time derivative defined. Please assign a derivative or a numerical scheme (e.g., BDF1).";
            }
            return m_timeDerivExpr;
        }

        CfemInt getDimension() const override { return m_dim; }
        CfemInt precedence() const override{ return 10; }
        std::string toString() const override { return m_name; }
        ExpressionType getType() const override { return ExpressionType::Variable; }

    };

}