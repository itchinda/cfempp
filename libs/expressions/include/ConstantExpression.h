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
            : m_value(std::move(val)) {}

        ConstantExpression(CfemReal val)
            : m_value({val}) {}
        
        ConstantExpression(const Vector3D& val)
            : m_value({val.x, val.y, val.z}) {}
        
        ConstantExpression(const TensorO2& val)
            : m_value({val(0,0), val(0,1), val(0,2),
                       val(1,0), val(1,1), val(1,2),
                       val(2,0), val(1,2), val(2,2)}) {}

        int getDimension() const override { return m_value.size(); }
        bool isPointwise() const override { return true; }
        bool isConstant() const override {return true; }
        const std::vector<CfemReal>& getValues() const { return m_value; }

        void evaluate(CfemInt, const Vector3D &, std::span<CfemReal> out, EvalContext &) const override {
            std::copy(m_value.begin(), m_value.end(), out.begin());
        }

        std::shared_ptr<const CfemExpression> grad() const override {
            CfemInt dim = getDimension();
            if (dim == 1) return zeros::Zero3D;
            if (dim == 3) return zeros::ZeroTensor;
            return std::make_shared<ZeroExpression>(dim * 3);
        }
    };

}