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
#include "CfemExpression.h"
#include "BinaryExpression.h"
#include <memory>

namespace cfem::sym {

    class ScaledExpression : public CfemExpression {
    public:
        ScaledExpression(CfemReal coeff, std::shared_ptr<const CfemExpression> expr)
            : m_coeff(coeff), m_base(std::move(expr))
        {
            CFEM_ASSERT(m_base, "ScaledExpression requires a valid base CfemExpression.");
        }

        virtual void prepare() const override { m_base->prepare(); }

        virtual void evaluate(CfemInt cellId, const Vector3D &xi, 
                              std::span<CfemReal> out, EvalContext &ctx) const override {

            m_base->evaluate(cellId, xi, out, ctx);
            CfemInt dim = m_base->getNumComponents();
            for (CfemInt i = 0; i < dim; ++i) {
                out[i] *= m_coeff; 
            }
        }

        inline void evaluateBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const override
        {

            CfemInt dim = nQp * m_base->getNumComponents();
            // Un tenseur 3x3 évalué sur nQp points nécessite 9 * nQp espaces mémoire
            CFEM_ASSERT(outValues.size() >= dim, "Batch evaluation output buffer too small for 3D Tensor");

            // =========================================================================
            // FORMAT SoA (Structure of Arrays) - Optimisé pour le cache et le SIMD
            // Le buffer de sortie sera organisé ainsi :
            // [xx_0..xx_n, xy_0..xy_n, xz_0..xz_n, yx_0..yx_n, yy_0..yy_n, yz_0..yz_n, zx_0..zx_n, zy_0..zy_n, zz_0..zz_n]
            // =========================================================================

            m_base->evaluateBatch(cellId, nQp, outValues, ctx);
            for (CfemInt i = 0; i < dim; ++i) {
                outValues[i] *= m_coeff; 
            }
        }

        CfemInt getNumComponents() const override { return m_base->getNumComponents(); }

        bool isContinuous() const override {
            return  m_base->isContinuous();
        }
        bool dependsOnFEField() const override { return m_base->dependsOnFEField(); }
        bool dependsOnShape() const override { return m_base->dependsOnShape(); }
        bool dependsOnGeometricNormal() const override {  return m_base->dependsOnGeometricNormal(); }
        
        ShapeEvaluationFlags requiredShapeFlags() const override {
            return m_base->requiredShapeFlags();
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            if (m_base) m_base->collectRequiredSpaces(spaces);
        }

        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            m_base->bind(spaceToId);
        }
        
        
        // Helper accessors for the "Smart Builder"
        CfemReal getCoeff() const { return m_coeff; }
        const CfemExpression* getBasePtr() const { return m_base.get(); }
        std::shared_ptr<const CfemExpression> getBaseShared() const { return m_base; }

        std::shared_ptr<const CfemExpression> div() const override{
            return std::make_shared<ScaledExpression>(m_coeff, m_base->div());
        }

        // std::string toString() const override {
        //     if (m_coeff == -1.0) return "-(" + m_base->toString() + ")";
        //     return "(" + std::to_string(m_coeff) + " * " + m_base->toString() + ")";
        // }

        ExpressionType getType() const override { return ExpressionType::Scaled; }

        std::string toString() const override{
            std::string baseStr = m_base->toString();
            
            if (m_base->getType() == ExpressionType::Binary) {
                auto bin = std::static_pointer_cast<const BinaryExpression>(m_base);
                if (bin->getOp() == MathOp::Add || bin->getOp() == MathOp::Sub) {
                    baseStr = "(" + baseStr + ")";
                }
            }

            if (NumericPolicy::areEqual(m_coeff, 1.0)){
                return baseStr;
            }
            else if (NumericPolicy::areEqual(m_coeff, -1.0)){
                return "-" + baseStr;
            }
            
            std::string coeffStr = formatCoefficient(m_coeff);
            return coeffStr + " * " + baseStr;
        }

    protected:
        CfemReal m_coeff;
        std::shared_ptr<const CfemExpression> m_base;

        std::shared_ptr<const CfemExpression> createGradient() const override {
            CFEM_ASSERT(m_base, "ScaledExpression::grad: base expression must be valid");

            auto baseGrad = m_base->grad();
            CFEM_ASSERT(baseGrad, "ScaledExpression:grad: base gradient is invalid!" );
            if (!baseGrad) return nullptr; // Or return a ZeroExpression of correct dim
            
            return std::make_shared<ScaledExpression>(m_coeff, baseGrad);
        }
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override{
            CFEM_ASSERT(m_base, "ScaledExpression::createTimeFDerivative: base expression must be valid");
            return std::make_shared<ScaledExpression>(m_coeff, m_base->timeDerivative());
        }

        void doCollectSpaceRequirements(SpaceRequirements& reqs, 
                                        ShapeEvaluationFlags inheritedFlags) const override 
        {
            if (m_base) m_base->collectSpaceRequirements(reqs, inheritedFlags);
        }
    };
}