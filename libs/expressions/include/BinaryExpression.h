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
        std::shared_ptr<const CfemExpression> m_lhs;
        std::shared_ptr<const CfemExpression> m_rhs;
        MathOp                                m_op;
        CfemInt                               m_dimL;
        CfemInt                               m_dimR;
        CfemInt                               m_dimOut;

    public:
        BinaryExpression(std::shared_ptr<const CfemExpression> lhs,
                         std::shared_ptr<const CfemExpression> rhs,
                         MathOp                                op);

        CfemInt getDimension() const override { return m_dimOut; }
        void prepare() const override{
            m_lhs->prepare();
            m_rhs->prepare();
        }

        bool isPointwise() const override { return m_lhs->isPointwise() && m_rhs->isPointwise(); }
        bool isBinary() const override { return true; }

        void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override;

        std::shared_ptr<const CfemExpression> getLHS() const { return m_lhs; }
        std::shared_ptr<const CfemExpression> getRHS() const { return m_rhs; }
        MathOp getOp() const { return m_op; }

        std::shared_ptr<const CfemExpression> grad() const override;
    };

}
#include "BinaryExpression.tpp"