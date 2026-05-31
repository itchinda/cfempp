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

#include "cfemutils.h"

namespace cfem
{
    namespace kernel_utils
    {
        // ========================================================================
        // Harware parameters
        // ========================================================================
        
        static constexpr size_t CACHELINE_SIZE = 64; ///< Required alignment for SIMD operations (cache-line aware).

        /**
         * @brief Calcule un "stride" (pas) garantissant l'alignement sur une ligne de cache.
         * Essentiel pour l'auto-vectorisation AVX2/AVX-512/SVE et pour éviter
         * que des variables indépendantes ne partagent la même cacheline.
         */
        inline constexpr CfemInt computePaddedStride(CfemInt nNodes) noexcept
        {
            constexpr CfemInt alignment = CACHELINE_SIZE / sizeof(CfemReal);
            return ((nNodes + alignment - 1) / alignment) * alignment;
        }

    }
}