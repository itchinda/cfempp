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
#include "CfemExpression.h"
#include "LambdaExpression.h"
#include <span>
#include <memory>


namespace cfem::sym
{

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class GradExpression : public CfemExpression
    {
    private:
        std::shared_ptr<const CfemExpression> m_expr;

        std::shared_ptr<const CfemExpression> createGradient() const override{
            CFEM_THROW("GradExpression does not yet support grad(). Cannot compute spatial derivatives for GradExpression.");
            return nullptr; 
        }
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override{
            return std::make_shared<GradExpression>(m_expr->timeDerivative());
        }

    public:
        explicit GradExpression(std::shared_ptr<const CfemExpression> expr) : m_expr(std::move(expr)) {
            if (auto lambda = std::dynamic_pointer_cast<const LambdaExpressionBase>(expr) ){
                if (!(lambda->grad()))
                    CFEM_THROW("LambdaExpressions are currentely not supported by the automatic differentation.");
            }
        }

        bool isContinuous() const override {
            return  m_expr->getType() == ExpressionType::FEField? true : m_expr->isContinuous();
        }
        CfemInt getNumComponents() const override { return m_expr->getNumComponents() * 3; }
        void prepare() const override { m_expr->prepare(); }
        void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override;
        ExpressionType getType() const override { return ExpressionType::Gradient; }
        bool dependsOnFEField() const override { return m_expr->dependsOnFEField(); }
        bool dependsOnShape() const override { return m_expr->dependsOnShape(); }
        bool dependsOnNormalGeometricNormal() const override { 
            return m_expr->dependsOnNormalGeometricNormal(); 
        }

        ShapeEvaluationFlags requiredShapeFlags() const override {
            if (!dependsOnShape())
                return ShapeEvaluationFlags::None;

            return m_expr->requiredShapeFlags() | ShapeEvaluationFlags::FirstDerivatives;
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            if (m_expr) m_expr->collectRequiredSpaces(spaces);
        }

        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            if (m_expr) m_expr->bind(spaceToId);
        }

        std::string toString() const override{ return "grad(" + m_expr->toString() +")"; }
    
    };


}