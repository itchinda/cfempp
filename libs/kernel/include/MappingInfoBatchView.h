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
#include "kernel_utils.h"
#include "EvaluationFlags.h"
#include <span>

namespace cfem 
{
    /**
     * @struct MappingInfoBatchView
     * @brief A zero-copy, non-owning view of metric data for a batch of quadrature points.
     * 
     * This structure is designed for High-Performance Computing (HPC). It does not own any 
     * memory. Instead, it holds lightweight `std::span` objects pointing to the global geometric 
     * transformation data (owned by a master MappingInfoBatch) and raw pointers to the temporary 
     * memory arena (`EvalContext`) for space-specific physical derivatives.
     * 
     * By enforcing a Structure of Arrays (SoA) layout and SIMD-aligned memory padding, this view 
     * ensures optimal L1 cache utilization and compiler auto-vectorization during finite element integration.
     */
    struct MappingInfoBatchView 
    {
        // ====================================================================
        // METADATA (POD)
        // ====================================================================
        CfemInt numQuadraturePoints = 0;                                  ///< Number of quadrature points in the current batch.
        CfemInt numNodes = 0;                                             ///< Number of basis functions (nodes) for the current FE space.
        MappingEvaluationFlags validFlags = MappingEvaluationFlags::None; ///< Bitmask indicating which metric fields are currently computed and valid.
        bool isCellAffine = false;                                        ///< True if the cell transformation is affine (constant Jacobian).

        // ====================================================================
        // SHARED GEOMETRIC VIEWS (Zero-Copy from MappingInfoBatch)
        // ====================================================================
        std::span<const Vector3D> physicalPositions;            ///< Physical coordinates of the quadrature points.
        std::span<const Tensor3D> invJacobians;                 ///< Inverse of the geometric Jacobian matrix per quadrature point.
        std::span<const CfemReal> detJs;                        ///< Determinant of the geometric Jacobian per quadrature point.
        std::span<const Vector3D> normals;                      ///< Physical normal vectors (for boundary/surface integration).

        // ====================================================================
        // SPACE-SPECIFIC DERIVATIVES (Pointers to EvalContext Arena)
        // ====================================================================
        
        // Pointers are mutable here to allow initialization, but exposed as const to the user via getters.
        CfemReal* dN_dx_data = nullptr;                         ///< Base pointer to physical x-derivatives.
        CfemReal* dN_dy_data = nullptr;                         ///< Base pointer to physical y-derivatives.
        CfemReal* dN_dz_data = nullptr;                         ///< Base pointer to physical z-derivatives.

        CfemReal* d2N_dxx_data = nullptr;                       ///< Base pointer to physical xx-second derivatives.
        CfemReal* d2N_dyy_data = nullptr;                       ///< Base pointer to physical yy-second derivatives.
        CfemReal* d2N_dzz_data = nullptr;                       ///< Base pointer to physical zz-second derivatives.
        CfemReal* d2N_dxy_data = nullptr;                       ///< Base pointer to physical xy-second derivatives.
        CfemReal* d2N_dyz_data = nullptr;                       ///< Base pointer to physical yz-second derivatives.
        CfemReal* d2N_dxz_data = nullptr;                       ///< Base pointer to physical xz-second derivatives.

        CfemInt paddedNodeStride = 0;                           ///< Memory stride between quadrature points, aligned for SIMD vectorization.
        
        // ====================================================================
        // ACCESSORS
        // ====================================================================
        
        /** 
         * @name Read-Only Accessors
         * @brief Safe getters used by FEFields and integration loops to read mathematical derivatives.
         * Memory layout: `data[q * paddedNodeStride + i]`.
         * @param q The index of the local quadrature point.
         * @return A contiguous, SIMD-aligned constant pointer to the derivatives for all nodes at point q.
         */
        ///@{
        [[nodiscard]] inline const CfemReal* dN_dx(CfemInt q) const { return dN_dx_data + q * paddedNodeStride; }
        [[nodiscard]] inline const CfemReal* dN_dy(CfemInt q) const { return dN_dy_data + q * paddedNodeStride; }
        [[nodiscard]] inline const CfemReal* dN_dz(CfemInt q) const { return dN_dz_data + q * paddedNodeStride; }
        
        [[nodiscard]] inline const CfemReal* d2N_dxx(CfemInt q) const { return d2N_dxx_data + q * paddedNodeStride; }
        [[nodiscard]] inline const CfemReal* d2N_dyy(CfemInt q) const { return d2N_dyy_data + q * paddedNodeStride; }
        [[nodiscard]] inline const CfemReal* d2N_dzz(CfemInt q) const { return d2N_dzz_data + q * paddedNodeStride; }
        [[nodiscard]] inline const CfemReal* d2N_dxy(CfemInt q) const { return d2N_dxy_data + q * paddedNodeStride; }
        [[nodiscard]] inline const CfemReal* d2N_dyz(CfemInt q) const { return d2N_dyz_data + q * paddedNodeStride; }
        [[nodiscard]] inline const CfemReal* d2N_dxz(CfemInt q) const { return d2N_dxz_data + q * paddedNodeStride; }
        ///@}

        /** 
         * @name Mutable Accessors
         * @brief Internal getters used EXCLUSIVELY by mapping routines (e.g., computePhysicalFirstDerivatives) 
         * to populate the memory arena.
         * @param q The index of the local quadrature point.
         * @return A contiguous, SIMD-aligned mutable pointer.
         */
        ///@{
        [[nodiscard]] inline CfemReal* dN_dx_mut(CfemInt q) { return dN_dx_data + q * paddedNodeStride; }
        [[nodiscard]] inline CfemReal* dN_dy_mut(CfemInt q) { return dN_dy_data + q * paddedNodeStride; }
        [[nodiscard]] inline CfemReal* dN_dz_mut(CfemInt q) { return dN_dz_data + q * paddedNodeStride; }
        ///@}

        // ====================================================================
        // ARENA ALLOCATION
        // ====================================================================

        /**
         * @brief Allocates contiguous, SIMD-aligned memory in the temporary EvalContext arena 
         * for space-specific physical derivatives.
         * 
         * This method handles the offset distribution to enforce a pure Structure of Arrays (SoA) layout.
         * It must ONLY be called for FE Space views, never for the pure geometric view.
         * 
         * @param scratchStorage Contiguous memory bloc allocated by the caller.
         * @param spaceNumNodes The number of basis functions in the target FE space.
         * @param spatialDim The spatial dimension of the mesh (1, 2, or 3).
         * @param nQp The number of quadrature points in the current batch.
         */
        inline void setupPhysicalDerivativesMemory(std::span<CfemReal> scratchStorage, 
                                                   CfemInt spaceNumNodes, 
                                                   CfemInt spatialDim, 
                                                   CfemInt nQp)
        {
            this->numNodes = spaceNumNodes;
            this->paddedNodeStride = kernel_utils::computePaddedStride(spaceNumNodes);
            
            CfemInt reqSize = spatialDim * nQp * this->paddedNodeStride;
            CFEM_ASSERT(scratchStorage.size() >= static_cast<size_t>(reqSize), 
                        "Insufficient memory provided for physical derivatives");

            this->dN_dx_data = scratchStorage.data();
            if (spatialDim >= 2) this->dN_dy_data = scratchStorage.data() + (1 * nQp * this->paddedNodeStride);
            if (spatialDim == 3) this->dN_dz_data = scratchStorage.data() + (2 * nQp * this->paddedNodeStride);

            this->validFlags |= MappingEvaluationFlags::PhysicalFirstDerivatives;
        }

        /**
         * @brief Quick check to verify if physical first derivatives have been computed and allocated.
         * @return True if valid, false otherwise.
         */
        [[nodiscard]] inline bool hasPhysicalFirstDerivatives() const noexcept 
        {
            return hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives);
        }
    };

} // namespace cfem