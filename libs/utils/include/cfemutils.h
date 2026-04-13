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

#ifndef CFEMUTILS_
#define CFEMUTILS_

#include <petscsys.h>
#include <vector>
#include <type_traits>
#include "ContractException.h"
#include "CfemLogger.h"
#include "ChronoManager.h"
#include <sstream>
#include <iomanip>
//#include "numeric_policy.h"

// CFEM TYPES
using CfemInt = int;
using CfemIntFlag = uint8_t;
using CfemReal = double;
    CFEM_STATIC_ASSERT((std::is_same_v<CfemReal,double>), "SymbolicEvaluator requires CfemReal = double.");

/**
 * @brief Namspace for CFEM++ types and objects
 *
 */
namespace cfem
{
    /**
     * @brief Threshold for numerical zero comparisons.
     * * Used to handle floating-point noise in geometric checks, such as
     * determining if a point lies on an element face or if a Jacobian determinant
     * is effectively zero.
     */
    constexpr CfemReal ZERO_TOLERANCE = 1e-14;
    inline constexpr CfemReal pi = static_cast<CfemReal>(3.141592653589793238462643383279502884L);

    /**
     * @brief The smallest representable difference for the @ref CfemReal type.
     * * Provides the theoretical lower bound of precision for the current
     * floating-point precision (float or double). Useful for scaling tolerances.
     */
    constexpr CfemReal MACHINE_EPSILON = std::numeric_limits<CfemReal>::epsilon();

    /**
     * @brief Internal identifier representing the global simulation universe.
     * * This ID is used to tag entities or regions that belong to the overall 
     * domain background, typically used when an entity is not associated with 
     * a specific subdomain or physical material.
     */
    constexpr CfemInt UNIVERSE_ID = -1;

    /**
     * @brief Internal identifier reserved for the mesh boundary entity.
     * This ID distinguishes entities that lie on the geometric frontier of the 
     * mesh. It is used during assembly and boundary condition mapping to 
     * identify nodes or faces that belong to the external mesh surface.
     */
    constexpr CfemInt SYS_BOUNDARY_ID = -2;


    struct NumericPolicy
    {
        static constexpr CfemReal relative_tol =
            1e3 * std::numeric_limits<CfemReal>::epsilon();

        static constexpr CfemReal absolute_tol =
            1e2 * std::numeric_limits<CfemReal>::epsilon();

        static bool isZero(CfemReal value, CfemReal scale = 1.0)
        {
            return std::abs(value) <=
                   relative_tol * scale + absolute_tol;
        }

        static bool areEqual(CfemReal a, CfemReal b, CfemReal scale = 1.0)
        {
            return std::abs(a - b) <=
                   relative_tol * scale + absolute_tol;
        }
    };



    inline std::string formatCoefficient(CfemReal val) {
        std::ostringstream out;
        // 'noshowpoint' removes the decimal if it's not needed (e.g., 7.0 -> 7)
        // 'precision' controls how many significant digits to keep
        out << std::noshowpoint << std::setprecision(15) << val;
        return out.str();
    }

} // namespace

#endif