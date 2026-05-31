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

#pragma once

#include <type_traits>

#include "cfemutils.h"

namespace cfem
{
    /**
     * @brief Bitmask flags to control the evaluation of shape functions on the reference element.
     * 
     * These flags dictate which quantities (values, gradients, etc.) are computed 
     * purely in the reference coordinate system (xi, eta, zeta).
     */
    enum class ShapeEvaluationFlags : CfemIntFlag
    {
        None              = 0,       ///< No fields requested.
        Value             = 1 << 0,  ///< Evaluate shape function values (N_i).
        FirstDerivatives  = 1 << 1,  ///< Evaluate first derivatives w.r.t reference coordinates (dN/dxi, dN/deta, dN/dzeta).
        SecondDerivatives = 1 << 2,  ///< Evaluate second derivatives w.r.t reference coordinates (Hessians).
        
        All = Value | FirstDerivatives | SecondDerivatives ///< Request all available reference fields.
    };

    /**
     * @brief Bitmask flags to control the computation of geometric and physical mapping fields.
     * 
     * These flags utilize a cascading bitwise structure to automatically resolve mathematical 
     * dependencies. For instance, requesting the InverseJacobian will implicitly flag the 
     * Determinant and the Jacobian for computation.
     */
    enum class MappingEvaluationFlags : CfemIntFlag
    {
        None                      = 0,                                   ///< No fields requested. Useful for clearing state. 
        PhysicalPosition          = 1 << 0,                              ///< Compute the interpolated physical coordinates (x, y, z).
        Jacobian                  = 1 << 1,                              ///< Compute the Jacobian matrix (dx/dxi).
        Determinant               = (1 << 2) | Jacobian,                 ///< Compute the determinant of the Jacobian (requires Jacobian).
        InverseJacobian           = (1 << 3) | Determinant,              ///< Compute the inverse Jacobian matrix (requires Determinant).
        PhysicalFirstDerivatives  = (1 << 4) | InverseJacobian,          ///< Compute shape gradients in the physical domain (dN/dx).
        PhysicalSecondDerivatives = (1 << 5) | PhysicalFirstDerivatives, ///< Compute shape Hessians in the physical domain.
        Normal                    = (1 << 6) | Jacobian,                 ///< Compute the outward unit normal vector (for boundaries/interfaces).

        All = PhysicalPosition | Jacobian | Determinant | InverseJacobian | PhysicalFirstDerivatives | PhysicalSecondDerivatives | Normal      ///< Request all available physical mapping fields.
    };
    
    // Bitwise Operators for Enumerations 
    // Note: SFINAE is used to restrict these operators to CFEM specific enums only,
    // preventing global namespace pollution for other standard enums.

    template<typename Enum, typename = std::enable_if_t<std::is_same_v<Enum, ShapeEvaluationFlags> || 
                                                        std::is_same_v<Enum, MappingEvaluationFlags>>>
    constexpr Enum operator|(Enum a, Enum b) 
    {
        return static_cast<Enum>( static_cast<CfemIntFlag>(a) | static_cast<CfemIntFlag>(b) );
    }

    template<typename Enum, typename = std::enable_if_t<std::is_same_v<Enum, ShapeEvaluationFlags> || 
                                                        std::is_same_v<Enum, MappingEvaluationFlags>>>
    constexpr Enum operator&(Enum a, Enum b) 
    {
        return static_cast<Enum>( static_cast<CfemIntFlag>(a) & static_cast<CfemIntFlag>(b) );
    }

    template<typename Enum, typename = std::enable_if_t<std::is_same_v<Enum, ShapeEvaluationFlags> || 
                                                        std::is_same_v<Enum, MappingEvaluationFlags>>>
    constexpr Enum& operator|=(Enum& a, Enum b) 
    {
        a = a | b;
        return a;
    }

    template<typename Enum, typename = std::enable_if_t<std::is_same_v<Enum, ShapeEvaluationFlags> || 
                                                        std::is_same_v<Enum, MappingEvaluationFlags>>>
    constexpr Enum& operator&=(Enum& a, Enum b)
    {
        a = a & b;
        return a;
    }

    /**
     * @brief Checks if a specific flag (or combination of flags) is fully set within a bitmask.
     * 
     * @tparam Enum The enumeration type.
     * @param flags The current state bitmask.
     * @param test The specific flag(s) to check for.
     * @return true if ALL bits of 'test' are present in 'flags', false otherwise.
     */
    template<typename Enum, typename = std::enable_if_t<std::is_same_v<Enum, ShapeEvaluationFlags> || 
                                                        std::is_same_v<Enum, MappingEvaluationFlags>>>
    [[nodiscard]] constexpr bool hasFlag(Enum flags, Enum test)
    {
        // STRICT CHECK: Ensure all bits of the 'test' flag are present in 'flags'.
        // This is crucial for cascading flags where testing a child should not 
        // mistakenly return true if only the parent dependency is present.
        return (static_cast<CfemIntFlag>(flags) & static_cast<CfemIntFlag>(test)) == static_cast<CfemIntFlag>(test);
    }


    /**
     * @brief Explicitly translates reference element evaluation prerequisites into their corresponding
     * physical geometric mapping requirements.
     * @param shapeFlags Active bitmask representing required reference element quantities.
     * @return MappingEvaluationFlags The minimized resolved physical mapping bitmask.
     */
    [[nodiscard]] constexpr MappingEvaluationFlags transformShapeToMappingFlags(ShapeEvaluationFlags shapeFlags)
    {
        MappingEvaluationFlags result = MappingEvaluationFlags::None;

        // Physical positions (x,y,z) are only computed if explicitly requested by a spatial symbolic node
        // if (hasFlag(shapeFlags, ShapeEvaluationFlags::Value))
        // {
        // }

        if (hasFlag(shapeFlags, ShapeEvaluationFlags::FirstDerivatives))
        {
            result |= MappingEvaluationFlags::PhysicalFirstDerivatives;
        }
        
        if (hasFlag(shapeFlags, ShapeEvaluationFlags::SecondDerivatives))
        {
            result |= MappingEvaluationFlags::PhysicalSecondDerivatives;
        }

        return result;
    }

    [[nodiscard]] constexpr ShapeEvaluationFlags transformMappingToShapeFlags(MappingEvaluationFlags mappingFlags)
    {
        ShapeEvaluationFlags result = ShapeEvaluationFlags::None;

        // Physical positions (x,y,z) are only computed if explicitly requested by a spatial symbolic node
        if (hasFlag(mappingFlags, MappingEvaluationFlags::PhysicalPosition))
        {
            result |= ShapeEvaluationFlags::Value;
        }

        if (hasFlag(mappingFlags, MappingEvaluationFlags::Jacobian) ||
            hasFlag(mappingFlags, MappingEvaluationFlags::Determinant) ||
            hasFlag(mappingFlags, MappingEvaluationFlags::InverseJacobian) ||
            hasFlag(mappingFlags, MappingEvaluationFlags::PhysicalFirstDerivatives) ||
            hasFlag(mappingFlags, MappingEvaluationFlags::Normal) )
        {
            result |= ShapeEvaluationFlags::FirstDerivatives;
        }
        
        if (hasFlag(mappingFlags, MappingEvaluationFlags::PhysicalSecondDerivatives))
        {
            result |= ShapeEvaluationFlags::FirstDerivatives | ShapeEvaluationFlags::SecondDerivatives;
        }

        return result;
    }

} // namespace cfem
