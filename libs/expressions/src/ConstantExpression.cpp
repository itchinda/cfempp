#include "ConstantExpression.h"

namespace cfem::sym
{

    std::shared_ptr<const CfemExpression> ConstantExpression::div() const {
        CfemInt dim = getDimension();
            
        switch (dim) {
            case 1:
                CFEM_ERROR << "Divergence operator cannot be applied to a scalar constant expression: " << toString();
            case 3:
                return zeros::Zero1D; 
            case 9:
                return zeros::Zero3D;
            default:
                CFEM_ERROR << "Divergence operator not supported for ConstantExpression of dimension " << dim;
        }
        return nullptr;
    }

    std::string ConstantExpression::toString() const { 
        CfemInt dim = getDimension();
        if (dim==1) return formatCoefficient(m_value[0]);

        std::ostringstream out;
        CfemInt i = 0;
        out << "[";
        for (auto val : m_value){
            if (i < dim-1) out << val << ", ";
            ++i;
        }
        out << m_value[dim-1] << "]";
        return out.str();
    }

}