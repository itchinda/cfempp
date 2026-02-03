#ifndef CFEMUTILS_
#define CFEMUTILS_

#include <petscsys.h>

#include "ContractException.h"


using CfemInt = PetscInt;
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
    
}

#endif