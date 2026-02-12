#ifndef CFEMUTILS_
#define CFEMUTILS_

#include <petscsys.h>
#include <vector>
#include "ContractException.h"

// CFEM TYPES
using CfemInt = PetscInt;
using CfemIntFlag = uint8_t;
using CfemReal = PetscReal;

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


    // class Environment {
    // public:
    //     // The user calls this once in main
    //     static void initialize(int &argc, char **&argv, const char* help = nullptr) {
    //         instance(); // ensures life
    //         PetscInitialize(&argc, &argv, nullptr, help);
    //     }

    //     static void ensureInitialized() {
    //     PetscBool initialized = PETSC_FALSE;
    //     PetscInitialized(&initialized);
    //     if (!initialized) {
    //         throw std::runtime_error(
    //             "[CFEM ERROR] PETSc not initialized.\n"
    //             "You must call cfem::Environment::initialize(...) "
    //             "before using CFEM objects."
    //         );
    //     }
    // }

    // private:
    //     Environment() {} 
    //     ~Environment() { PetscFinalize(); }
        
    //     static Environment& instance() {
    //         static Environment env; // Guaranteed to be destroyed at program end
    //         return env;
    //     }
    // };

} // namespace

#endif