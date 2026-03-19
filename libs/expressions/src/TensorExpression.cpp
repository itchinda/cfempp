#include "TensorExpression.h"
#include "VectorExpression.h"


namespace cfem::sym
{

    /**
     * @brief Computes the divergence of the tensor expression.
     * @return A VectorExpression representing the row-wise divergence.
     */
    // std::shared_ptr<const CfemExpression> TensorExpression::div() const {
    //     // We assume your framework has a VectorExpression(expr_x, expr_y, expr_z)
    //     // and that each scalar expression has a .diff(int direction) or similar method.
        
    //     // auto div_x = m_components[0]->diff(0) + m_components[1]->diff(1) + m_components[2]->diff(2);
    //     // auto div_y = m_components[3]->diff(0) + m_components[4]->diff(1) + m_components[5]->diff(2);
    //     // auto div_z = m_components[6]->diff(0) + m_components[7]->diff(1) + m_components[8]->diff(2);

    //     CFEM_ERROR << "TODO: Divergence expression not implemented.";

    //     return nullptr; //std::make_shared<const VectorExpression>(div_x, div_y, div_z);
    // }
}