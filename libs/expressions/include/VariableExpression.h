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

#include "cfemutils.h"
#include "Vector3D.h"
#include "Tensor3D.h"
#include "CfemExpression.h"
#include "ZeroExpression.h"
#include "EvalContext.h"
#include <span>
#include <memory>

namespace cfem::sym
{

    class VariableExpression : public CfemExpression {
    protected:
        std::string m_name;
        const CfemReal* m_dataPtr;
        CfemInt m_dim;
        std::shared_ptr<const CfemExpression> m_gradExpr; 
        std::shared_ptr<const CfemExpression> m_timeDerivExpr;

        std::shared_ptr<const CfemExpression> createGradient() const override {
            if (!m_gradExpr) {
                CFEM_ERROR << "Gradient was not provided. Cannot compute spatial derivatives for VariableExpression.";
            }
            return m_gradExpr;
        }

        /**
         * @brief Returns the time derivative expression.
         * @throw CFEM_ERROR if no time derivative or numerical scheme has been assigned.
         */
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override {
            if (!m_timeDerivExpr) {
                CFEM_ERROR << "VariableExpression '" << m_name << 
                    "' has no time derivative defined. Please assign a derivative or a numerical scheme (e.g., BDF1).";
            }
            return m_timeDerivExpr;
        }

        void doCollectSpaceRequirements(SpaceRequirements& reqs, 
                                        ShapeEvaluationFlags inheritedFlags) const override 
        {
            ShapeEvaluationFlags gradFlags = ShapeEvaluationFlags::None;

            if (hasFlag(inheritedFlags, ShapeEvaluationFlags::FirstDerivatives)) {
                gradFlags |= ShapeEvaluationFlags::Value;
            }

            if (hasFlag(inheritedFlags, ShapeEvaluationFlags::SecondDerivatives)) {
                gradFlags |= ShapeEvaluationFlags::FirstDerivatives;
            }

            if (m_gradExpr && gradFlags != ShapeEvaluationFlags::None) {
                m_gradExpr->collectSpaceRequirements(reqs, gradFlags);
            }
            else if (!m_gradExpr && gradFlags != ShapeEvaluationFlags::None) {
                CFEM_THROW("Error in AST: Derivatives requested from a LambdaExpression, "
                           "but no analytical gradient expression was provided during its construction.");
            }
        }

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

        // std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const override {
        //     if (var_idx == 3) {
        //         return timeDerivative();
        //     }
            
        //     // 2. Spatial derivative (x=0, y=1, z=2)
        //     auto g = grad();
        //     if (!g) {
        //         CFEM_ERROR << "Cannot compute diff() because grad() is not set for variable: " << m_name;
        //     }
            
        //     // Assuming grad() returns a Vector Expression for a scalar variable:
        //     // grad() = [du/dx, du/dy, du/dz]. We just extract the 'var_idx' component.
        //     return extractComponent(g, var_idx); 
        // }


        /**
         * @brief Defines or updates the expression used for the time derivative.
         * Useful for assigning numerical schemes like (U - U_old) / dt.
         */
        void setTimeDerivative(std::shared_ptr<const CfemExpression> expr) {
            if (expr && expr->getNumComponents() != m_dim) {
                CFEM_ERROR << "VariableExpression::setTimeDerivative: Dimension mismatch.";
            }
            m_timeDerivExpr = std::move(expr);
        }

        /**
         * @brief Defines or updates the expression used for the time derivative.
         * Useful for assigning numerical schemes like (U - U_old) / dt.
         */
        void setGradient(std::shared_ptr<const CfemExpression> expr) {
            if (expr && expr->getNumComponents() != m_dim * 3) {
                CFEM_ERROR << "VariableExpression::setGradient: Dimension mismatch.";
            }
            m_gradExpr = std::move(expr);
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            if (m_gradExpr) m_gradExpr->collectRequiredSpaces(spaces); 
            if (m_timeDerivExpr) m_timeDerivExpr->collectRequiredSpaces(spaces);       
        }

        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            if (m_gradExpr) m_gradExpr->bind(spaceToId);
            if (m_timeDerivExpr) m_timeDerivExpr->bind(spaceToId);
        }

        CfemInt getNumComponents() const override { return m_dim; }
        CfemInt precedence() const override{ return 10; }
        std::string toString() const override { return m_name; }
        ExpressionType getType() const override { return ExpressionType::Variable; }

    };

}