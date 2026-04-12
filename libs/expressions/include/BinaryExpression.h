#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
#include "TensorO2.h"
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
    class BinaryExpression : public CfemExpression
    {
    private:
        std::shared_ptr<const CfemExpression> m_lhs;
        std::shared_ptr<const CfemExpression> m_rhs;
        MathOp                                m_op;
        CfemInt                               m_dimL;
        CfemInt                               m_dimR;
        CfemInt                               m_dimOut;
    
        std::shared_ptr<const CfemExpression> createGradient() const override;
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override;

    public:
        BinaryExpression(std::shared_ptr<const CfemExpression> lhs,
                         std::shared_ptr<const CfemExpression> rhs,
                         MathOp op);

        CfemInt getDimension() const override { return m_dimOut; }

        void prepare() const override
        {
            m_lhs->prepare();
            m_rhs->prepare();
        }

        bool isPointwise() const override { return m_lhs->isPointwise() && m_rhs->isPointwise(); }
        bool isBinary() const override { return true; }

        bool isContinuous() const override {
            return  m_lhs->isContinuous() && m_rhs->isContinuous();
        }
        bool dependsOnFEField() const override {
            return m_lhs->dependsOnFEField() || m_rhs->dependsOnFEField();
        }

        bool dependsOnShape() const override {
            return m_lhs->dependsOnShape() || m_rhs->dependsOnShape();
        }

        ShapeEvaluationFlags requiredShapeFlags() const override {
            return m_lhs->requiredShapeFlags() | m_rhs->requiredShapeFlags();
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            if (m_lhs) m_lhs->collectRequiredSpaces(spaces);
            if (m_rhs) m_rhs->collectRequiredSpaces(spaces);
            
        }

        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            if (m_lhs) m_lhs->bind(spaceToId);
            if (m_rhs) m_rhs->bind(spaceToId);
        }

        void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override;
        
        std::shared_ptr<const CfemExpression> getLHS() const { return m_lhs; }
        std::shared_ptr<const CfemExpression> getRHS() const { return m_rhs; }

        MathOp getOp() const { return m_op; }

        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const override;
        std::shared_ptr<const CfemExpression> div() const override;

        CfemInt precedence() const override;
        std::string toString() const override;
        ExpressionType getType() const override { return ExpressionType::Binary; }

    };

}

#include "BinaryExpression.tpp"