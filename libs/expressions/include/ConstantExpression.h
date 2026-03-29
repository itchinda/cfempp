#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
#include "TensorO2.h"
#include "EvalContext.h"
#include "CfemExpression.h"
#include "ZeroExpression.h"
#include <vector>
#include <span>
#include <memory>
#include <functional>


namespace cfem::sym
{

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class ConstantExpression : public CfemExpression
    {
        std::vector<CfemReal> m_value;  //!

    public:
        ConstantExpression(std::vector<CfemReal> val)
            : m_value(std::move(val)) 
        {}

        ConstantExpression(CfemReal val)
            : m_value({val}) 
        {}
        
        ConstantExpression(const Vector3D& val)
            : m_value({val.x, val.y, val.z}) 
        {}
        
        ConstantExpression(const TensorO2& val)
            : m_value({val(0,0), val(0,1), val(0,2),
                       val(1,0), val(1,1), val(1,2),
                       val(2,0), val(1,2), val(2,2)}) 
         {}

         ConstantExpression(CfemReal val, CfemInt dim)
            : m_value(dim, val)
        {}

        int getDimension() const override { return m_value.size(); }
        bool isPointwise() const override { return true; }
        bool isConstant() const override {return true; }
        bool isZero() const override {
            for (auto val : m_value) {
                if (! NumericPolicy::isZero(std::abs(val)) ) return false;
            }
            return true;
        }

        const std::vector<CfemReal>& getValues() const { return m_value; }

        void evaluate(CfemInt, const Vector3D &, std::span<CfemReal> out, EvalContext &) const override {
            std::copy(m_value.begin(), m_value.end(), out.begin());
        }

        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const override {
            return std::make_shared<ZeroExpression>(this->getDimension());
        }

        std::shared_ptr<const CfemExpression> grad() const override;
        std::shared_ptr<const CfemExpression> div() const override;
        std::shared_ptr<const CfemExpression> timeDerivative() const override {
            return std::make_shared<ZeroExpression>(getDimension());
        }

        CfemInt precedence() const override{ return 10; }

        ExpressionType getType() const override { return ExpressionType::Constant; }
        std::string toString() const override;


    };

}