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

        CfemInt getDimension() const override { return 3; }

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

        ExpressionType getType() const override { return ExpressionType::TensorExpression; }

        std::string toString() const override {
            std::ostringstream os;
            os << "[[" << m_components[0]->toString() << ", " << m_components[1]->toString() << ", " << m_components[2]->toString() << "], "
               <<  "[" << m_components[3]->toString() << ", " << m_components[4]->toString() << ", " << m_components[5]->toString() << "], "
               <<  "[" << m_components[6]->toString() << ", " << m_components[7]->toString() << ", " << m_components[8]->toString() << "]]";
            return os.str();
        }

        std::shared_ptr<const CfemExpression> timeDerivative() const override{
           return std::make_shared<const TensorExpression>( m_components[0]->timeDerivative(), m_components[1]->timeDerivative(), m_components[2]->timeDerivative(),
                                                            m_components[3]->timeDerivative(), m_components[4]->timeDerivative(), m_components[5]->timeDerivative(),
                                                            m_components[6]->timeDerivative(), m_components[7]->timeDerivative(), m_components[8]->timeDerivative() );
        };

        std::shared_ptr<const CfemExpression> grad() const{ 
            CFEM_ERROR << "The gradient of a tensor expression is not yet implemented. Would required a 3rd order tensor expression.";
            return nullptr;
        };

        /**
         * @brief Computes the divergence of the tensor expression.
         * @return A VectorExpression representing the row-wise divergence.
         */
        // std::shared_ptr<const CfemExpression> div() const override;

         
    };

}