#include "ZeroExpression.h"

namespace cfem::sym
{

    std::shared_ptr<const CfemExpression> ZeroExpression::div() const {
        switch (m_dim) {
            case 1:
                CFEM_ERROR << "Divergence operator cannot be applied to a scalar zero expression.";
            case 3:
                return zeros::Zero1D; 
            case 9:
                return zeros::Zero3D;
            default:
                CFEM_ERROR << "Divergence operator not supported for ZeroExpression of dimension " << m_dim;
        }
        return nullptr;
    }

    std::string ZeroExpression::toString() const { 
        CfemInt dim = getDimension();
        if (dim==1) return "0";
        
        std::ostringstream out;
        out << "[";
        for (CfemInt i = 0; i< m_dim-1; ++i){
            out << "0" << ", ";
        }
        out << "0" << "]";
        return out.str();
    }
}

