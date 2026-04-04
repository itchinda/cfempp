#include "TensorExpression.h"
#include "VectorExpression.h"
#include "MathOperators.h"


namespace cfem::sym
{

    std::shared_ptr<const CfemExpression> TensorExpression::div() const {
        auto div_x = m_components[0]->diff(0) + 
                     m_components[1]->diff(1) + 
                     m_components[2]->diff(2);

        auto div_y = m_components[3]->diff(0) + 
                     m_components[4]->diff(1) + 
                     m_components[5]->diff(2);

        auto div_z = m_components[6]->diff(0) + 
                     m_components[7]->diff(1) + 
                     m_components[8]->diff(2);

        return std::make_shared<const VectorExpression>(div_x, div_y, div_z);
    }

    std::shared_ptr<const CfemExpression> TensorExpression::timeDerivative() const {
       return std::make_shared<const TensorExpression>( m_components[0]->timeDerivative(), m_components[1]->timeDerivative(), m_components[2]->timeDerivative(),
                                                        m_components[3]->timeDerivative(), m_components[4]->timeDerivative(), m_components[5]->timeDerivative(),
                                                        m_components[6]->timeDerivative(), m_components[7]->timeDerivative(), m_components[8]->timeDerivative() );
    }
}