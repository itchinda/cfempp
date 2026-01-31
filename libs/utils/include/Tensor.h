/**
 * @file Tensor.h
 * @brief Physical tensor definitions for CFEM.
 */

#ifndef CFEM_TENSOR_H
#define CFEM_TENSOR_H

#include "StaticMatrix.h"
#include "SymTensorO2.h" // If you have the symmetric version

/**
 * @brief General second-order tensor (9 components).
 * Used for non-symmetric quantities like the Deformation Gradient (F)
 * or the Velocity Gradient (L).
 */
using TensorO2 = StaticMatrix<CfemReal, 3, 3>;

/**
 * @brief Identity tensor constant.
 */
inline TensorO2 IdentityTensor() {
    TensorO2 I;
    I.setIdentity();
    return I;
}

#endif