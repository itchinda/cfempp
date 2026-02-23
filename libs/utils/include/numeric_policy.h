#ifndef CFEM_NUMERIC_POLICY_H
#define CFEM_NUMERIC_POLICY_H

#include <limits>
#include "cfemutils.h"

namespace cfem
{
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
}

#endif