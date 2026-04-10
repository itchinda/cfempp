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
            CfemInt dim = m_base->getDimension();
            for (CfemInt i = 0; i < dim; ++i) {
                out[i] *= m_coeff; 
            }
        }

        CfemInt getDimension() const override { return m_base->getDimension(); }

        bool isContinuous() const override {
            return  m_base->isContinuous();
        }
        bool dependsOnFEField() const override { return m_base->dependsOnFEField(); }
        bool dependsOnShape() const override { return m_base->dependsOnShape(); }
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

        std::shared_ptr<const CfemExpression> grad() const override {
            CFEM_ASSERT(m_base, "ScaledExpression::grad: base expression must be valid");

            auto baseGrad = m_base->grad();
            CFEM_ASSERT(baseGrad, "ScaledExpression:grad: base gradient is invalid!" );
            if (!baseGrad) return nullptr; // Or return a ZeroExpression of correct dim
            
            return std::make_shared<ScaledExpression>(m_coeff, baseGrad);
        }
        
        
        std::shared_ptr<const CfemExpression> div() const override{
            return std::make_shared<ScaledExpression>(m_coeff, m_base->div());
        }
        
        std::shared_ptr<const CfemExpression> timeDerivative() const override{
            return std::make_shared<ScaledExpression>(m_coeff, m_base->timeDerivative());
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

    private:
        CfemReal m_coeff;
        std::shared_ptr<const CfemExpression> m_base;
    };
}