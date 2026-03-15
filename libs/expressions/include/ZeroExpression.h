#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
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
    class ZeroExpression : public CfemExpression
    {
        CfemInt m_dim {0};

    public:
        ZeroExpression (CfemInt dim = 0): m_dim(dim) {
        }

        int getDimension() const override { return m_dim; }
        bool isPointwise() const override { return true; }
        bool isZero() const override { return true; };
        bool isConstant() const override {return true; }

        void evaluate(CfemInt, const Vector3D &, std::span<CfemReal> out, EvalContext &) const override{
            std::ranges::fill(out, 0.0);
        }

        std::shared_ptr<const CfemExpression> grad() const {
            return std::make_shared<ZeroExpression>(m_dim * 3);
        }
        ExpressionType getType() const override { return ExpressionType::Zero; }

        std::string toString() const override { 
            CfemInt dim = getDimension();
            if (dim==1) return "0";

            std::ostringstream out;
            CfemInt i = 0;
            out << "[";
            for (CfemInt i = 0; i< m_dim-1; ++i){
                out << "0" << ", ";
                ++i;
            }
            out << "0" << "]";
            return out.str();
         }
    };

    namespace zeros{
        inline const auto Zero1D = std::make_shared<ZeroExpression>(1);
        inline const auto Zero2D = std::make_shared<ZeroExpression>(2);
        inline const auto Zero3D = std::make_shared<ZeroExpression>(3);
        inline const auto ZeroTensor = std::make_shared<ZeroExpression>(9);

        inline const std::shared_ptr<const CfemExpression>& Zero = Zero1D;
    }

}