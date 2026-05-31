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
    // ========================================================================
    // ShapeData : POD HPC Pur (Aucun pointeur stocké, calculs à la volée)
    // ========================================================================
    
    /**
     * @struct ShapeInfo
     * @brief Contiguous, cache-aligned storage buffer for element shape function evaluations.
     *
     * This structure implements a Structure-of-Arrays (SoA) layout to store values, gradients, 
     * and Hessians for a set of quadrature points.
     *
     * ### Performance Design
     * - **Contiguous Memory:** All fields are stored in a single `DynamicVector` allocation.
     * - **Cache Locality:** Optimized for L1/L2 cache prefetching by using padded strides.
     * - **Zero Overhead:** Accessors provide lightweight views (raw pointers) to data without 
     * requiring heavy abstractions in hot loops.
     * - **Memory Safety:** Designed as a move-only type (Rule of 0) to prevent dangling 
     * pointers resulting from accidental object cloning.
     */
    struct ShapeInfo
    {
        static constexpr CfemInt NUM_FIRST_ORDER_FIELDS  = 4; ///< Number of fields for first-order derivatives including the shapes values: \f$ \left\{N, \frac{\partial N}{\partial \xi}, \frac{\partial N}{\partial \eta}, \frac{\partial N}{\partial \zeta} \right} \f$
        static constexpr CfemInt NUM_SECOND_ORDER_FIELDS = 6; ///< Number of fields for second-order derivatives: \f$ \left\{\frac{\partial^2 N}{\partial \xi^2}, \frac{\partial^2 N}{\partial \eta^2}, \frac{\partial^2 N}{\partial \zeta^2}, \frac{\partial^2 N}{\partial \xi \partial \eta}, \frac{\partial^2 N}{\partial \xi \partial \zeta}, \frac{\partial^2 N}{\partial \eta\partial \zeta}  \right} \f$

        DynamicVector<CfemReal> storage;                      ///< Centralized contiguous buffer

        CfemInt numNodes = 0;                                 ///< Number of nodes per element
        CfemInt stride = 0;                                   ///< Padded stride to ensure SIMD alignment
        bool hasSecondDerivatives = false;                   

        bool shapeMustBeUpdated = true;
        ShapeEvaluationFlags validFlags = ShapeEvaluationFlags::None;

        /**
         * @brief Allocates and configures the internal storage layout.
         *
         * @param nNodes Number of interpolation nodes in the reference element.
         * @param includeSecondDerivatives If true, memory for Hessian-related quantities is allocated.
         *
         * @note
         * Calling this function invalidates all previously acquired raw field pointers.
         */
        inline void setup(CfemInt nNodes, bool includeSecondDerivatives = false)
        {
            numNodes = nNodes;
            hasSecondDerivatives = includeSecondDerivatives;
            stride = kernel_utils::computePaddedStride(numNodes);

            const CfemInt nFields = NUM_FIRST_ORDER_FIELDS + 
                                   (includeSecondDerivatives ? NUM_SECOND_ORDER_FIELDS : 0);

            const size_t totalSize = static_cast<size_t>(nFields) * static_cast<size_t>(stride);

            if (storage.size() != totalSize)
            {
                storage.resize(totalSize);
            }

            validFlags = ShapeEvaluationFlags::None;
        }

        // ================================================================================================
        // SHAPES VALUES
        // ================================================================================================

        /** @brief Returns a writable pointer to shape function values. */
        [[nodiscard]] inline CfemReal* values()  noexcept {
            return storage.data();
        }

        // ================================================================================================
        // SHAPES FIRST DERIVATIVES
        // ================================================================================================

        /** @brief Returns a writable pointer to shape derivatives along the reference \f$\xi\f$ direction. */
        [[nodiscard]] inline CfemReal* dN_dxi()  noexcept {
            return storage.data() + 1 * stride; 
        }

        /** @brief Returns a writable pointer to shape derivatives along the reference \f$\eta\f$ direction. */
        [[nodiscard]] inline CfemReal* dN_deta() noexcept {
            return storage.data() + 2 * stride; 
        }

        /** @brief Returns a writable pointer to shape derivatives along the reference \f$\zeta\f$ direction. */
        [[nodiscard]] inline CfemReal* dN_dzeta() noexcept {
            return storage.data() + 3 * stride;
        }

        // ================================================================================================
        // SHAPES SECOND DERIVATIVES
        // ================================================================================================

        /** @brief Returns a writable pointer to shape second derivatives \f$\partial^2 N \partial \zeta^2\f$ direction. */
        [[nodiscard]] inline CfemReal* d2N_dxi2()   noexcept   {
            return storage.data() + 4 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_deta2()  noexcept   {
            return storage.data() + 5 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dzeta2()  noexcept  {
            return storage.data() + 6 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxideta()  noexcept {
            return storage.data() + 7 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxidzeta() noexcept {
            return storage.data() + 8 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_detadzeta() noexcept {
            return storage.data() + 9 * stride;
        }


        // Versions const (pour la lecture)
        [[nodiscard]] inline const CfemReal* values() const {
            CFEM_ASSERT(has(ShapeEvaluationFlags::Value), "Shape values not computed?");
            return storage.data();
        }

        [[nodiscard]] inline const CfemReal* dN_dxi() const {
            CFEM_ASSERT( has(ShapeEvaluationFlags::FirstDerivatives), "Shape first derivatives not computed?");
            return storage.data() + 1 * stride;
        }

        [[nodiscard]] inline const CfemReal* dN_deta() const {
            CFEM_ASSERT( has(ShapeEvaluationFlags::FirstDerivatives), "Shape first derivatives not computed?");
            return storage.data() + 2 * stride;
        }

        [[nodiscard]] inline const CfemReal* dN_dzeta() const {
            CFEM_ASSERT( has(ShapeEvaluationFlags::FirstDerivatives), "Shape first derivatives not computed?");
            return storage.data() + 3 * stride;
        }


        [[nodiscard]] inline const CfemReal* d2N_dxi2() const {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not computed?");  
            return storage.data() + 4 * stride; 
        }

        [[nodiscard]] inline const CfemReal* d2N_deta2() const {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not computed?");  
            return storage.data() + 5 * stride; 
        }

        [[nodiscard]] inline const CfemReal* d2N_dzeta2() const {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not computed?");  
            return storage.data() + 6 * stride; 
        }

        [[nodiscard]] inline const CfemReal* d2N_dxideta() const {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not computed?");  
            return storage.data() + 7 * stride; 
        }

        [[nodiscard]] inline const CfemReal* d2N_dxidzeta() const {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not computed?");  
            return storage.data() + 8 * stride; 
        }

        [[nodiscard]] inline const CfemReal* d2N_detadzeta() const {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not computed?");  
            return storage.data() + 9 * stride; 
        }


        // ============================================================
        // Getters
        // ============================================================

        /**
         * @brief Returns a reference-space gradient vector for a specific node index.
         * @param i Node index (0 to numNodes-1).
         */
        inline Vector3D gradient(CfemInt i) const
        {
            CFEM_ASSERT(i >= 0 && i < numNodes);
            return Vector3D(dN_dxi()[i], dN_deta()[i], dN_dzeta()[i]);
        }

        /**
         * @brief Returns the symmetric Hessian tensor for a specific node index.
         * @param i Node index (0 to numNodes-1).
         */
        inline SymTensor3D hessian(CfemInt i) const
        {
            CFEM_ASSERT(i >= 0 && i < numNodes);
            CFEM_ASSERT(hasSecondDerivatives);
            return SymTensor3D(
                d2N_dxi2()[i], d2N_deta2()[i], d2N_dzeta2()[i],
                d2N_detadzeta()[i], d2N_dxidzeta()[i], d2N_dxideta()[i]
            );
        }

        inline bool has(ShapeEvaluationFlags field) const
        {
            return hasFlag(validFlags, field);
        }

        inline size_t storageSize() const
        {
            return storage.size();
        }

        inline size_t storageBytes() const
        {
            return storage.size() * sizeof(CfemReal);
        }
    };

}