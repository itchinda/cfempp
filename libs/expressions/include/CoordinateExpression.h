#pragma once

#include "CfemExpression.h"
#include "ConstantExpression.h"
#include "ZeroExpression.h"

namespace cfem::sym
{

    class CoordinateExpression : public CfemExpression
    {
        CfemInt m_dimIndex; // 0 for x, 1 for y, 2 for z and 3 for t 
        std::string m_name; // e.g. x, y, z, t

    public:
        CoordinateExpression(const std::string& name, CfemInt index) : m_name(name), m_dimIndex(index) {
            CFEM_ASSERT(m_dimIndex >=0 && m_dimIndex<4, "Coordinate index should be within 0, 1, 2, and 4" );
            NamingPolicy::validate(name);
            VariableRegistry::registerVariable(name);
        }

        void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override
        {
            if (m_dimIndex == 3) 
            {
                out[0] = ctx.currentTime; return;
            }

            const Vector3D& globalX = ctx.physPos;
            out[0] = globalX[m_dimIndex];
        }

        CfemInt getDimension() const override { return 1;}

        std::shared_ptr<const CfemExpression> grad() const override
        {            
            Vector3D g{0, 0, 0};
            if (m_dimIndex == 3) return std::make_shared<ConstantExpression>(g); // grad(t) = [0,0,0]

            g[m_dimIndex] = 1.0;
            return std::make_shared<ConstantExpression>(g);
        }

        std::shared_ptr<const CfemExpression> timeDerivative() const override {
            if (m_name == "t" || m_dimIndex == 3) {
                return std::make_shared<ConstantExpression>(1.0);
            }
            return std::make_shared<ZeroExpression>(1);
        }

        std::string toString() const override { return m_name; }
        ExpressionType getType() const override { return ExpressionType::ScalarSymbolic; }
    };

}