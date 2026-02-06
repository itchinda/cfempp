#ifndef CFEMUTILS_
#define CFEMUTILS_

#include <petscsys.h>

#include "ContractException.h"


using CfemInt = PetscInt;
using CfemIntFlag = uint8_t;
using CfemReal = PetscReal;

        
#ifdef NDEBUG
    #define CFEM_ERROR(msg) \
        throw std::runtime_error(std::string("[CFEM ERROR] ") + msg)

    #define CFEM_WARNING(msg) \
        std::cerr << "[CFEM WARNING] " << msg << std::endl
#else
    #define CFEM_ERROR(msg) \
        throw std::runtime_error(std::string("[CFEM ERROR] ") + msg + \
        " File: " + __FILE__ + ":" + std::to_string(__LINE__))

    #define CFEM_WARNING(msg) \
        std::cerr << "[CFEM WARNING] " << msg << \
        " (File: " << __FILE__ << ":" << __LINE__ << ")" << std::endl
#endif

// cfemutils.h
namespace cfem {
    constexpr CfemReal ZERO_TOLERANCE = 1e-14;
    constexpr CfemReal MACHINE_EPSILON = std::numeric_limits<CfemReal>::epsilon();
    constexpr CfemInt UNIVERSE_ID = -1;
    constexpr CfemInt SYS_BOUNDARY_ID = -2;

    inline bool areNearlyEqual(CfemReal a, CfemReal b) {
        // Relative comparison: scales the tolerance based on the magnitude of the numbers
        return std::abs(a - b) <= ( (std::abs(a) < std::abs(b) ? std::abs(b) : std::abs(a)) * MACHINE_EPSILON);
    }

    /** @brief Énumération des champs activables pour les fonctions de forme. */
    enum ShapeActiveField : CfemIntFlag {
        None        = 0,
        Value       = 1 << 0,
        Gradient    = 1 << 1,
        Hessian     = 1 << 2,
        All         = 0xFF
    };
    
}

#endif