#include "LambdaExpressionBase.h"

namespace cfem::sym
{
    LambdaExpressionBase::LambdaExpressionBase(CfemInt dimension,
                                               std::shared_ptr<const CfemExpression> gradient,
                                               std::shared_ptr<const CfemExpression> timeDeriv): 
        m_dimension(dimension)
    {
        
        m_gradient = std::move(gradient);
        m_timeDeriv = std::move(timeDeriv);
        if (m_gradient){
            if (m_dimension == 1 && m_gradient->getDimension() != 3){
                CFEM_ERROR << "The gradient of a scalar lambda expression must be a 3D vector.";
            }
            if (m_dimension == 3 && m_gradient->getDimension() != 9){
                CFEM_ERROR << "The gradient of a vector lambda expression must be a 2nd order tensor (dimension 9).";
            }
        }

        if (m_timeDeriv){
            if (m_timeDeriv->getDimension() != m_dimension){
                CFEM_ERROR << "The time derivative must have the same dimension ("
                           << m_dimension << ") as the lambda expression.";
            }
        }
    }
}