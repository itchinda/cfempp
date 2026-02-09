#ifndef CFEM_SHAPE_EVALUATION_FLAGS_H
#define CFEM_SHAPE_EVALUATION_FLAGS_H

#include "cfemutils.h"

namespace cfem
{
    /**
     * @brief Bitmask flags to control which shape function quantities are computed.
     * These flags are used in @ref computeMapping to optimize performance by avoiding
     * unnecessary calculations. For example, a simple mass matrix assembly may only
     * require @ref Value, while a stiffness matrix requires @ref Gradient.
     */
    enum class ShapeEvaluationFlags : CfemIntFlag
    {
        /// No fields requested. Useful for clearing state. */
        None = 0,

        /// Compute shape function values \f$ N_i \f$. Required for physical position and value interpolation.
        Value = 1 << 0,

        /// Compute shape function gradients \f$ \nabla N_i \f$. Required for Jacobian and stiffness matrices.
        Gradient = 1 << 1,

        /// Compute shape function Hessians \f$ \nabla^2 N_i \f$. Required for high-order geometric analysis or specific physics.
        Hessian = 1 << 2,

        /// Compute all available fields (Values, Gradients, and Hessians).
         All = Value | Gradient | Hessian
    };

    
    /// Enable bitwise OR for ShapeEvaluationFlags
    inline ShapeEvaluationFlags operator|( ShapeEvaluationFlags a, ShapeEvaluationFlags b)
    {
        return static_cast<ShapeEvaluationFlags>( static_cast<CfemIntFlag>(a) | static_cast<CfemIntFlag>(b));
    }

    /// Enable bitwise AND for ShapeEvaluationFlags
    inline ShapeEvaluationFlags operator&( ShapeEvaluationFlags a, ShapeEvaluationFlags b)
    {
        return static_cast<ShapeEvaluationFlags>( static_cast<CfemIntFlag>(a) & static_cast<CfemIntFlag>(b) );
    }

    /// Test if a flag is set
    inline bool hasFlag(ShapeEvaluationFlags flags, ShapeEvaluationFlags test)
    {
        return (flags & test) != ShapeEvaluationFlags::None;
    }


} //namespace
#endif