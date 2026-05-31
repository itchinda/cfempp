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
#include "Tensor3D.h"
#include "DynamicVector.h"
#include "DynamicMatrix.h"

namespace cfem
{
    /**
     * @struct ShapeInfoBatch
     * @brief Contiguous, cache-aligned storage buffer for shape function evaluations across a batch of quadrature points.
     *
     * ### HPC Design (Structure of Arrays - SoA)
     * To maximize SIMD auto-vectorization and cache hits, all evaluations (values, gradients, Hessians) 
     * for all quadrature points and all nodes are flattened into a single 1D `DynamicVector`.
     * 
     * Memory Layout for a given quadrature point `q`:
     * `[ N_nodes... | dxi_nodes... | deta_nodes... | dzeta_nodes... | (Hessians...) ]`
     * 
     * Each field block is padded to ensure it starts on a SIMD-friendly memory boundary (e.g., 64 bytes).
     * Lightweight view accessors (e.g., `dN_dxi(q)`) return raw pointers to these aligned blocks.
     */
    struct ShapeInfoBatch
    {
        /** @brief Number of fields for first-order data: {N, dN/dxi, dN/deta, dN/dzeta} */
        static constexpr CfemInt NUM_FIRST_ORDER_FIELDS  = 4; 
        
        /** @brief Number of fields for second-order data: {d2N/dxi2, d2N/deta2, d2N/dzeta2, d2N/dxideta, d2N/dxidzeta, d2N/detadzeta} */
        static constexpr CfemInt NUM_SECOND_ORDER_FIELDS = 6; 

        CfemInt numQuadraturePoints = 0;
        CfemInt numNodes = 0;
        
        bool shapeMustBeUpdated = true;
        ShapeEvaluationFlags validFlags = ShapeEvaluationFlags::None;

        // ====================================================================
        // HPC Contiguous Storage
        // ====================================================================
        DynamicVector<CfemReal> storage; ///< Centralized contiguous buffer for all QPs and all fields
        
        CfemInt paddedNodeStride = 0;    ///< Padded stride per field to ensure SIMD alignment
        CfemInt qpStride = 0;            ///< Total stride per quadrature point (all fields combined)
        bool hasSecondDerivatives = false;

        /** 
         * @brief Prepares buffer sizes and strides to prevent runtime reallocations.
         * @param nQp Number of quadrature points in the batch.
         * @param nNodes Number of interpolation nodes in the reference element.
         * @param includeSecondDerivatives If true, memory for Hessian-related quantities is allocated.
         */
        inline void setup(CfemInt nQp, CfemInt nNodes, bool includeSecondDerivatives = false)
        {
            numQuadraturePoints = nQp;
            numNodes = nNodes;
            hasSecondDerivatives = includeSecondDerivatives;

            paddedNodeStride = kernel_utils::computePaddedStride(numNodes);
            const CfemInt nFields = NUM_FIRST_ORDER_FIELDS + 
                                   (includeSecondDerivatives ? NUM_SECOND_ORDER_FIELDS : 0);
            
            qpStride = nFields * paddedNodeStride;
            const size_t totalStorageSize = static_cast<size_t>(nQp) * static_cast<size_t>(qpStride);

            if (storage.size() != totalStorageSize)
            {
                storage.resize(totalStorageSize);
            }

            validFlags = ShapeEvaluationFlags::None;
        }

        // ================================================================================================
        // LIGHTWEIGHT VIEWS (Writable) - Return raw pointers for a specific Quadrature Point 'q'
        // ================================================================================================

        [[nodiscard]] inline CfemReal* values(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 0 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* dN_dxi(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 1 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* dN_deta(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 2 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* dN_dzeta(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 3 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxi2(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 4 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_deta2(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 5 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dzeta2(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 6 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxideta(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 7 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxidzeta(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 8 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_detadzeta(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 9 * paddedNodeStride;
        }

        // ================================================================================================
        // LIGHTWEIGHT VIEWS (Read-Only) - Include validFlags safety checks
        // ================================================================================================

        [[nodiscard]] inline const CfemReal* values(CfemInt q) const {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(has(ShapeEvaluationFlags::Value), "Shape values not computed");
            return storage.data() + q * qpStride + 0 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* dN_dxi(CfemInt q) const {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(has(ShapeEvaluationFlags::FirstDerivatives), "First derivatives not computed");
            return storage.data() + q * qpStride + 1 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* dN_deta(CfemInt q) const {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(has(ShapeEvaluationFlags::FirstDerivatives), "First derivatives not computed");
            return storage.data() + q * qpStride + 2 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* dN_dzeta(CfemInt q) const {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(has(ShapeEvaluationFlags::FirstDerivatives), "First derivatives not computed");
            return storage.data() + q * qpStride + 3 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dxi2(CfemInt q) const {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasSecondDerivatives && has(ShapeEvaluationFlags::SecondDerivatives), "Second derivatives not computed");
            return storage.data() + q * qpStride + 4 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* d2N_deta2(CfemInt q) const {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasSecondDerivatives && has(ShapeEvaluationFlags::SecondDerivatives), "Second derivatives not computed");
            return storage.data() + q * qpStride + 5 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dzeta2(CfemInt q) const {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasSecondDerivatives && has(ShapeEvaluationFlags::SecondDerivatives), "Second derivatives not computed");
            return storage.data() + q * qpStride + 6 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dxideta(CfemInt q) const {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasSecondDerivatives && has(ShapeEvaluationFlags::SecondDerivatives), "Second derivatives not computed");
            return storage.data() + q * qpStride + 7 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dxidzeta(CfemInt q) const {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasSecondDerivatives && has(ShapeEvaluationFlags::SecondDerivatives), "Second derivatives not computed");
            return storage.data() + q * qpStride + 8 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* d2N_detadzeta(CfemInt q) const {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasSecondDerivatives && has(ShapeEvaluationFlags::SecondDerivatives), "Second derivatives not computed");
            return storage.data() + q * qpStride + 9 * paddedNodeStride;
        }

        // ============================================================
        // Accessors & Getters
        // ============================================================

        /** 
         * @brief Reconstructs the reference gradient vector for a specific quadrature point and node.
         * @param q Quadrature point index.
         * @param i Node index.
         * @return Vector3D containing the reference gradients.
         */
        inline Vector3D gradient(CfemInt q, CfemInt i) const
        {
            CFEM_ASSERT(i >= 0 && i < numNodes, "Invalid node index");
            return Vector3D(dN_dxi(q)[i], dN_deta(q)[i], dN_dzeta(q)[i]);
        }

        /** 
         * @brief Reconstructs the reference Hessian tensor for a specific quadrature point and node. 
         * @param q Quadrature point index.
         * @param i Node index.
         * @return SymTensor3D containing the reference second derivatives.
         */
        inline SymTensor3D hessian(CfemInt q, CfemInt i) const
        {
            CFEM_ASSERT(i >= 0 && i < numNodes, "Invalid node index");
            CFEM_ASSERT(has(ShapeEvaluationFlags::SecondDerivatives), "Shape second derivatives not computed?");
            return SymTensor3D(d2N_dxi2(q)[i], d2N_deta2(q)[i], d2N_dzeta2(q)[i], 
                               d2N_detadzeta(q)[i], d2N_dxidzeta(q)[i], d2N_dxideta(q)[i]);
        }

        /**
         * @brief Checks if a specific evaluation field is marked as valid.
         */
        inline bool has(ShapeEvaluationFlags field) const
        {
            return hasFlag(validFlags, field);
        }
    };
    
}