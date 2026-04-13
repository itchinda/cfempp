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
#include <sstream>

namespace cfem::sym
{

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class TensorExpression : public CfemExpression
    {
    private:
        std::shared_ptr<const CfemExpression> m_components[9];
        std::shared_ptr<const CfemExpression> createGradient() const override { 
            CFEM_ERROR << "Gradient is not yet supported for tensor expressions. It would required a 3rd order tensor which is not implemented.";
            return nullptr; 
        }
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override;

    public:
        TensorExpression (std::shared_ptr<const CfemExpression> expr_xx, std::shared_ptr<const CfemExpression> expr_xy, std::shared_ptr<const CfemExpression> expr_xz,
                          std::shared_ptr<const CfemExpression> expr_yx, std::shared_ptr<const CfemExpression> expr_yy, std::shared_ptr<const CfemExpression> expr_yz,
                          std::shared_ptr<const CfemExpression> expr_zx, std::shared_ptr<const CfemExpression> expr_zy, std::shared_ptr<const CfemExpression> expr_zz)
        { 
            m_components[0] = std::move(expr_xx);
            m_components[1] = std::move(expr_xy);
            m_components[2] = std::move(expr_xz);
            m_components[3] = std::move(expr_yx);
            m_components[4] = std::move(expr_yy);
            m_components[5] = std::move(expr_yz);
            m_components[6] = std::move(expr_zx);
            m_components[7] = std::move(expr_zy);
            m_components[8] = std::move(expr_zz);

            for (int i = 0; i < 9; ++i) {
                if (!m_components[i] || m_components[i]->getDimension() != 1) {
                    CFEM_ERROR << "TensorExpression component " << i << " must be a valid scalar expression.";
                }
            }

        }

        CfemInt getDimension() const override { return 9; }

        bool isPointwise() const override{
            return m_components[0]->isPointwise() && m_components[1]->isPointwise() && m_components[2]->isPointwise() &&
                   m_components[3]->isPointwise() && m_components[4]->isPointwise() && m_components[5]->isPointwise() &&
                   m_components[6]->isPointwise() && m_components[7]->isPointwise() && m_components[8]->isPointwise();
        }

        bool isZero() const override{
            return m_components[0]->isZero() && m_components[1]->isZero() && m_components[2]->isZero() &&
                   m_components[3]->isZero() && m_components[4]->isZero() && m_components[5]->isZero() &&
                   m_components[6]->isZero() && m_components[7]->isZero() && m_components[8]->isZero();
        }

        bool isConstant() const override{
            return m_components[0]->isConstant() && m_components[1]->isConstant() && m_components[2]->isConstant() &&
                   m_components[3]->isConstant() && m_components[4]->isConstant() && m_components[5]->isConstant() &&
                   m_components[6]->isConstant() && m_components[7]->isConstant() && m_components[8]->isConstant();
        };

        inline void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override
        {
            CFEM_ASSERT(out.size() >= 9, "Expected a buffer of size >=9.");

            for (CfemInt i=0; i<9; ++i){
                m_components[i]->evaluate(cellId, xi, out.subspan(i, 1), ctx);
            }
        }

        inline std::shared_ptr<const CfemExpression> getComponent(int i) const
        {
            CFEM_ASSERT(i < 9 && i >= 0, "Index out of bounds.");
            return m_components[i];
        }

        bool dependsOnFEField() const override { 
            return m_components[0]->dependsOnFEField() || 
                   m_components[1]->dependsOnFEField() ||
                   m_components[2]->dependsOnFEField() ||
                   m_components[3]->dependsOnFEField() || 
                   m_components[4]->dependsOnFEField() ||
                   m_components[5]->dependsOnFEField() ||
                   m_components[6]->dependsOnFEField() || 
                   m_components[7]->dependsOnFEField() ||
                   m_components[8]->dependsOnFEField(); 
        }

        bool dependsOnShape() const override {
            return m_components[0]->dependsOnShape() || 
                   m_components[1]->dependsOnShape() ||
                   m_components[2]->dependsOnShape() ||
                   m_components[3]->dependsOnShape() || 
                   m_components[4]->dependsOnShape() ||
                   m_components[5]->dependsOnShape() ||
                   m_components[6]->dependsOnShape() || 
                   m_components[7]->dependsOnShape() ||
                   m_components[8]->dependsOnShape(); 
        
        }

        ShapeEvaluationFlags requiredShapeFlags() const override {
            return m_components[0]->requiredShapeFlags() | 
                   m_components[1]->requiredShapeFlags() | 
                   m_components[2]->requiredShapeFlags() | 
                   m_components[3]->requiredShapeFlags() | 
                   m_components[4]->requiredShapeFlags() | 
                   m_components[5]->requiredShapeFlags() | 
                   m_components[6]->requiredShapeFlags() | 
                   m_components[7]->requiredShapeFlags() | 
                   m_components[8]->requiredShapeFlags();
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            m_components[0]->collectRequiredSpaces(spaces); 
            m_components[1]->collectRequiredSpaces(spaces); 
            m_components[2]->collectRequiredSpaces(spaces); 
            m_components[3]->collectRequiredSpaces(spaces); 
            m_components[4]->collectRequiredSpaces(spaces); 
            m_components[5]->collectRequiredSpaces(spaces); 
            m_components[6]->collectRequiredSpaces(spaces); 
            m_components[7]->collectRequiredSpaces(spaces); 
            m_components[8]->collectRequiredSpaces(spaces);        
        }

        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            m_components[0]->bind(spaceToId);
            m_components[1]->bind(spaceToId);
            m_components[2]->bind(spaceToId); 
            m_components[3]->bind(spaceToId); 
            m_components[4]->bind(spaceToId); 
            m_components[5]->bind(spaceToId); 
            m_components[6]->bind(spaceToId); 
            m_components[7]->bind(spaceToId); 
            m_components[8]->bind(spaceToId);
        }

        bool isContinuous() const override {
            return  m_components[0]->isContinuous() &&
                    m_components[1]->isContinuous() &&
                    m_components[2]->isContinuous() && 
                    m_components[3]->isContinuous() && 
                    m_components[4]->isContinuous() && 
                    m_components[5]->isContinuous() && 
                    m_components[6]->isContinuous() && 
                    m_components[7]->isContinuous() && 
                    m_components[8]->isContinuous();
        }

        ExpressionType getType() const override { return ExpressionType::TensorExpression; }

        std::string toString() const override {
            std::ostringstream os;
            os << "[[" << m_components[0]->toString() << ", " << m_components[1]->toString() << ", " << m_components[2]->toString() << "], "
               <<  "[" << m_components[3]->toString() << ", " << m_components[4]->toString() << ", " << m_components[5]->toString() << "], "
               <<  "[" << m_components[6]->toString() << ", " << m_components[7]->toString() << ", " << m_components[8]->toString() << "]]";
            return os.str();
        }

        std::shared_ptr<const CfemExpression> div() const override;

    };

}