//************************************************************************
// --- CFEM++ Library
// ---
// --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.
// --- ALL RIGHTS RESERVED
// ---
// --- This software is protected by international copyright laws.
// --- Use, distribution, or modification of this software in any form,
// --- source or binary, for personal, academic, or non-commercial 
// --- purposes is permitted free of charge, provided that this 
// --- copyright notice and this permission notice appear in all 
// --- copies and supporting documentation.
// ---
// --- The authors and contributors provide this code "as is" without 
// --- any express or implied warranty. In no event shall they be 
// --- held liable for any damages arising from the use of this software.
//************************************************************************

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
        None             = 0,

        /// Compute shape function values \f$ N_i \f$. Required for physical position and value interpolation.
        Value            = 1 << 0,
        
        // Compute the physical position
        PhysicalPosition = (1 << 1) | Value,
        
        /// Compute the Jacobian of the mapping from the reference element to the physical element. 
        Jacobian         = 1 << 2,
        
        /// Compute shape function gradients \f$ \nabla N_i \f$. Required for Jacobian.
        Gradient         = (1 << 3) | Jacobian,
        
        /// Compute shape function Hessians \f$ \nabla^2 N_i \f$. Required Jacobian and Shape Values
        Hessian          = (1 << 4) | Gradient,

        /// Compute all available fields (Values, Gradients, and Hessians).
        All              = Value | PhysicalPosition | Jacobian | Gradient | Hessian
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

    inline ShapeEvaluationFlags& operator|=(ShapeEvaluationFlags& a, ShapeEvaluationFlags b)
    {
        a = a | b;
        return a;
    }

    inline ShapeEvaluationFlags& operator&=(ShapeEvaluationFlags& a, ShapeEvaluationFlags b)
    {
        a = a & b;
        return a;
    }
    
    inline bool hasFlag(ShapeEvaluationFlags flags, ShapeEvaluationFlags test)
    {
        return static_cast<CfemIntFlag>(flags & test) != 0;
    }

} //namespace
#endif