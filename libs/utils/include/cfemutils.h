#ifndef CFEMUTILS_
#define CFEMUTILS_

#include <petscsys.h>
#include <vector>
#include <type_traits>
#include "ContractException.h"
#include "ChronoManager.h"
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
     * @brief Sentinel value indicating an unassigned or invalid ID.
     * * Commonly used to initialize indices before they are mapped to the mesh
     * or to represent a "null" relationship in connectivity arrays.
     */
    constexpr CfemInt UNIVERSE_ID = -1;

    /**
     * @brief Sentinel value reserved for system-level boundary conditions.
     * * Used to distinguish internal mesh interfaces from external boundaries
     * where constraints (Dirichlet/Neumann) are applied.
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

} // namespace

#endif