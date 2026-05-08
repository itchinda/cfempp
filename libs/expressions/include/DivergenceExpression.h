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
#include "EvalContext.h"
#include <span>
#include <memory>
#include <vector>

namespace cfem::sym
{
    /** *************************************************************************
     * @brief Dedicated node to compute the divergence of a Vector or Tensor.
     * Evaluated lazily at runtime using partial spatial derivatives.
     **************************************************************************** */
    class DivergenceExpression : public CfemExpression
    {
    private:
        std::shared_ptr<const CfemExpression> m_expr;
        CfemInt m_dimIn;
        CfemInt m_dimOut;

        // Cache the partial derivative expressions so we only build the tree once!
        std::shared_ptr<const CfemExpression> m_dx;
        std::shared_ptr<const CfemExpression> m_dy;
        std::shared_ptr<const CfemExpression> m_dz;

        std::shared_ptr<const CfemExpression> createTimeDerivative() const override {
            // d(div(u)) / dt == div( du / dt )
            return m_expr->timeDerivative()->div();
        }

        std::shared_ptr<const CfemExpression> createGradient() const override {
            CFEM_THROW("The gradient of a divergence node requires a dedicated GradExpression fallback.");
            return nullptr;
        }

    public:
        DivergenceExpression(std::shared_ptr<const CfemExpression> expr)
            : m_expr(std::move(expr))
        {
            m_dimIn = m_expr->getDimension();
            
            if (m_dimIn == 1) {
                CFEM_THROW("Divergence cannot be applied to a scalar field: " << m_expr->toString());
            } else if (m_dimIn == 3) {
                m_dimOut = 1; // Vector -> Scalar
            } else if (m_dimIn == 9) {
                m_dimOut = 3; // Tensor -> Vector
            } else {
                CFEM_THROW("Unsupported dimension for divergence: " << m_dimIn);
            }

            m_dx = m_expr->diff(0);
            m_dy = m_expr->diff(1);
            m_dz = m_expr->diff(2);
        }

        CfemInt getDimension() const override { return m_dimOut; }
        bool isPointwise() const override { return m_expr->isPointwise(); }
        bool dependsOnFEField() const override { return m_expr->dependsOnFEField(); }
        bool dependsOnShape() const override { return m_expr->dependsOnShape(); }
        
        ShapeEvaluationFlags requiredShapeFlags() const override {
            if (!dependsOnShape())
                return ShapeEvaluationFlags::None;

            return m_expr->requiredShapeFlags() | ShapeEvaluationFlags::Gradient; 
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            if (m_expr) m_expr->collectRequiredSpaces(spaces);
        }

        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            if (m_expr) m_expr->bind(spaceToId);
        }

        ExpressionType getType() const override { return ExpressionType::Divergence; }

        std::string toString() const override {
            return "div(" + m_expr->toString() + ")";
        }

        void evaluate(CfemInt cellIdx, const Vector3D& point, std::span<CfemReal> out, EvalContext& ctx) const override 
        {
            CFEM_ASSERT(out.size() >= m_dimOut, "DivergenceExpression: 'out' span is too small!");

            CfemReal buf_dx[9] = {0};
            CfemReal buf_dy[9] = {0};
            CfemReal buf_dz[9] = {0};

            m_dx->evaluate(cellIdx, point, std::span{buf_dx, (size_t)m_dimIn}, ctx);
            m_dy->evaluate(cellIdx, point, std::span{buf_dy, (size_t)m_dimIn}, ctx);
            m_dz->evaluate(cellIdx, point, std::span{buf_dz, (size_t)m_dimIn}, ctx);

            if (m_dimIn == 3) {
                // div(u) = dux/dx + duy/dy + duz/dz
                out[0] = buf_dx[0] + buf_dy[1] + buf_dz[2];
            } 
            else if (m_dimIn == 9) {
                out[0] = buf_dx[0] + buf_dy[1] + buf_dz[2]; // Row 1
                out[1] = buf_dx[3] + buf_dy[4] + buf_dz[5]; // Row 2
                out[2] = buf_dx[6] + buf_dy[7] + buf_dz[8]; // Row 3
            }
        }

        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const override {
            // By Clairaut's theorem (symmetry of second derivatives):
            // d(div(u)) / dx_i == div( du / dx_i )
            return m_expr->diff(var_idx)->div(); 
        }
    };
}