#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
#include "CfemExpression.h"
#include "ScalarSymbolicExpression.h"
#include "EvalContext.h"
#include <span>
#include <memory>

#include <symengine/expression.h>
#include <symengine/lambda_double.h>

#include "SymbolicEvaluator.h"


namespace cfem::sym
{

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class TensorSymbolicExpression : public CfemExpression
    {
    public:
        TensorSymbolicExpression(
            std::shared_ptr<const ScalarSymbolicExpression> dxx, std::shared_ptr<const ScalarSymbolicExpression> dxy, std::shared_ptr<const ScalarSymbolicExpression> dxz,
            std::shared_ptr<const ScalarSymbolicExpression> dyx, std::shared_ptr<const ScalarSymbolicExpression> dyy, std::shared_ptr<const ScalarSymbolicExpression> dyz,
            std::shared_ptr<const ScalarSymbolicExpression> dzx, std::shared_ptr<const ScalarSymbolicExpression> dzy, std::shared_ptr<const ScalarSymbolicExpression> dzz);

        CfemInt getDimension() const override { return 9; }
        bool isPointwise() const override { return true; }

        void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
        {
            CFEM_ASSERT(out.size() >= 9, "TensorSymbolicExpression: 'out' span is too small!");
            for (int i = 0; i < 9; ++i) m_components[i]->evaluate(cellId, xi, out.subspan(i, 1), ctx);
        }

        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const override;

        // TODO Implement grad()

        std::shared_ptr<const ScalarSymbolicExpression> getComponent(CfemInt i) const
        {
            CFEM_ASSERT(i>=0 && i< 9, "TensorSymbolicExpression: Index out of bounds.");
            return m_components[i];
        }

        ExpressionType getType() const override { return ExpressionType::TensorSymbolic; }
        std::shared_ptr<const CfemExpression> timeDerivative() const override;


        std::string toString() const override {
            std::ostringstream os;
            os << "[[" << m_components[0]->toString() << ", " << m_components[1]->toString() << ", " << m_components[2]->toString() << "], "
               <<  "[" << m_components[3]->toString() << ", " << m_components[4]->toString() << ", " << m_components[5]->toString() << "], "
               <<  "[" << m_components[6]->toString() << ", " << m_components[7]->toString() << ", " << m_components[8]->toString() << "]]";
            return os.str();
        }
    private:
        std::shared_ptr<const ScalarSymbolicExpression> m_components[9];
    };

}
