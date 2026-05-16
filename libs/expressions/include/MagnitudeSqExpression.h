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
#include "EvalContext.h"
#include <span>
#include <memory>

namespace cfem::sym
{


    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class MagnitudeSqExpression : public CfemExpression
    {
    private:
        CfemInt m_innerDim = 0;
        std::shared_ptr<const CfemExpression> m_expr;
        std::shared_ptr<const CfemExpression> createGradient() const override;
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override;
        

    public:
        explicit MagnitudeSqExpression(std::shared_ptr<const CfemExpression> expr)
            : m_innerDim(expr->getNumComponents()), m_expr(std::move(expr))
        {
            if (!m_expr) CFEM_ERROR << "MagnitudeSqExpression requires a valid inner expression.";
        }

        CfemInt getNumComponents() const override { return 1; }

        bool isContinuous() const override {
            return  m_expr->isContinuous();
        }

        bool dependsOnFEField() const override { return m_expr->dependsOnFEField(); }
        bool dependsOnShape() const override { return m_expr->dependsOnShape(); }
        ShapeEvaluationFlags requiredShapeFlags() const override {
            return m_expr->requiredShapeFlags();
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            if (m_expr) m_expr->collectRequiredSpaces(spaces);
        }

        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            m_expr->bind(spaceToId);
        }

        void prepare() const override { m_expr->prepare(); }

        void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override;

        std::shared_ptr<const CfemExpression> getInner() const { return m_expr; }


        CfemInt precedence() const override { return 5; }
        ExpressionType getType() const override { return ExpressionType::MagnitudeSq; }
        std::string toString() const override { return "mag_sqr(" + m_expr->toString() + ")"; }
    };

}

#include "MagnitudeSqExpression.tpp"