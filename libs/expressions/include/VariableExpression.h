#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
#include "TensorO2.h"
#include "CfemExpression.h"
#include "ConstantExpression.h"
#include "EvalContext.h"
#include <span>
#include <memory>

namespace cfem::sym
{

    class VariableExpression : public CfemExpression {
    private:
        const CfemReal* m_dataPtr;
        CfemInt m_dim;
        std::shared_ptr<const CfemExpression> m_gradExpr; 

    public:
        VariableExpression(const CfemReal* ptr, 
                        CfemInt dim, 
                        std::shared_ptr<const CfemExpression> gradExpr = nullptr) 
            : m_dataPtr(ptr), m_dim(dim), m_gradExpr(gradExpr) {}

        void evaluate(CfemInt cellId, const Vector3D& xi, std::span<CfemReal> out, EvalContext& ctx) const override {
            for (CfemInt i = 0; i < m_dim; ++i) out[i] = m_dataPtr[i];
        }

        std::shared_ptr<const CfemExpression> grad() const override {
            if (!m_gradExpr) {
                // Si aucun gradient n'est fourni, on suppose que c'est une constante 
                // par rapport aux coordonnées (ex: un paramètre physique fixe)
                return std::make_shared<ConstantExpression>(3); 
            }
            return m_gradExpr;
        }

        CfemInt getDimension() const override { return m_dim; }
    };

}