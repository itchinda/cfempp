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

#include <memory>
#include <span>
#include <sstream>

#include "CfemExpression.h"
#include "cfemutils.h"
#include "EvalContext.h"
#include "Vector3D.h"

namespace cfem::sym
{

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class VectorExpression : public CfemExpression
    {
    protected:
        std::shared_ptr<const CfemExpression> m_components[3];
        std::shared_ptr<const CfemExpression> createGradient() const override;
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override;

        void doCollectSpaceRequirements(SpaceRequirements& reqs, 
                                        ShapeEvaluationFlags inheritedFlags) const override 
        {
            m_components[0]->collectSpaceRequirements(reqs, inheritedFlags); 
            m_components[1]->collectSpaceRequirements(reqs, inheritedFlags); 
            m_components[2]->collectSpaceRequirements(reqs, inheritedFlags);  
        }

    public:
        VectorExpression (std::shared_ptr<const CfemExpression> expr_x,
                          std::shared_ptr<const CfemExpression> expr_y,
                          std::shared_ptr<const CfemExpression> expr_z)
        { 
            m_components[0] = std::move(expr_x);
            m_components[1] = std::move(expr_y);
            m_components[2] = std::move(expr_z);

            for (int i = 0; i < 3; ++i) {
                if (!m_components[i] || m_components[i]->getNumComponents() != 1) {
                    CFEM_THROW("VectorExpression component " << i << " must be a valid scalar expression. Got dim = " << m_components[i]->getNumComponents() );
                }
            }

        }

        CfemInt getNumComponents() const override { return 3; }

        bool isPointwise() const override{
            return m_components[0]->isPointwise() && m_components[1]->isPointwise() && m_components[2]->isPointwise();
        }

        bool isZero() const override{
            return m_components[0]->isZero() && m_components[1]->isZero() && m_components[2]->isZero();
        }

        bool isConstant() const override{
            return m_components[0]->isConstant() && m_components[1]->isConstant() && m_components[2]->isConstant();
        };

        inline void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override
        {
            CFEM_ASSERT(out.size() >= 3, "VectorExpression: Expected a buffer of size >= 3.");

            for (CfemInt i = 0; i < 3; ++i){
                m_components[i]->evaluate(cellId, xi, out.subspan(i, 1), ctx);
            }
        }

        inline void evaluateBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const override
        {
            CFEM_ASSERT(outValues.size() >= 3 * nQp, "Batch evaluation output buffer too small for 3D Vector");

            // =========================================================================
            // FORMAT SoA (Structure of Arrays) - Optimisé pour le cache et le SIMD
            // [x_0, x_1, ..., x_n,  y_0, y_1, ..., y_n,  z_0, z_1, ..., z_n]
            // =========================================================================
            
            m_components[0]->evaluateBatch(cellId, nQp, outValues.subspan(0 * nQp, nQp), ctx);
            m_components[1]->evaluateBatch(cellId, nQp, outValues.subspan(1 * nQp, nQp), ctx);
            m_components[2]->evaluateBatch(cellId, nQp, outValues.subspan(2 * nQp, nQp), ctx);
        }

        inline std::shared_ptr<const CfemExpression> getComponent(int i) const
        {
            CFEM_ASSERT(i < 3 && i >= 0, "Index out of bounds.");
            return m_components[i];
        }

        bool dependsOnFEField() const override { 
            return m_components[0]->dependsOnFEField() || 
                   m_components[1]->dependsOnFEField() ||
                   m_components[2]->dependsOnFEField(); 
        }

        bool dependsOnShape() const override { 
            return m_components[0]->dependsOnShape() || 
                   m_components[1]->dependsOnShape() ||
                   m_components[2]->dependsOnShape(); 
        }

        bool dependsOnGeometricNormal() const override { 
            return m_components[0]->dependsOnGeometricNormal() || 
                   m_components[1]->dependsOnGeometricNormal() ||
                   m_components[2]->dependsOnGeometricNormal(); 
        }

        ShapeEvaluationFlags requiredShapeFlags() const override {
            return m_components[0]->requiredShapeFlags() | 
                   m_components[1]->requiredShapeFlags() | 
                   m_components[2]->requiredShapeFlags();
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            m_components[0]->collectRequiredSpaces(spaces); 
            m_components[1]->collectRequiredSpaces(spaces); 
            m_components[2]->collectRequiredSpaces(spaces);   
        }

        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            m_components[0]->bind(spaceToId);
            m_components[1]->bind(spaceToId);
            m_components[2]->bind(spaceToId);
        }

        bool isContinuous() const override{ 
            return m_components[0]->isContinuous()&&
                   m_components[1]->isContinuous()&&
                   m_components[2]->isContinuous();
        }

        ExpressionType getType() const override { return ExpressionType::VectorExpression; }

        std::string toString() const override {
            std::ostringstream os;
            os << "[" << m_components[0]->toString() << ", " << m_components[1]->toString() << ", " << m_components[2]->toString() << "]";
            return os.str();
        }

        std::shared_ptr<const CfemExpression> div() const override;
        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const {
            return std::make_shared<const VectorExpression>(m_components[0]->diff(var_idx), m_components[1]->diff(var_idx), m_components[2]->diff(var_idx));
        }

         
    };

}