#pragma once

#include "cfemutils.h"
#include "kernel_utils.h"
#include "EvaluationFlags.h"
#include "DynamicMatrix.h"
#include "DynamicMatrix.h"
#include "Tensor3D.h"
#include "Tensor2D.h"

#include "ShapeInfo.h"


namespace cfem
{
/**
     * @struct MappingInfo
     * @brief Contiguous, cache-aligned storage buffer for physical shape function derivatives and mapping data.
     *
     * This structure implements a Structure-of-Arrays (SoA) layout to store physical gradients
     * and Hessians for a set of quadrature points.
     * * ### Performance Design
     * - **Contiguous Memory:** Derivatives are stored in a single `DynamicVector` allocation.
     * - **Cache Locality:** Optimized for L1/L2 cache prefetching by using padded strides.
     * - **Zero Overhead:** Accessors provide lightweight views (raw pointers) to data without 
     * requiring heavy abstractions in hot loops.
     */
    struct MappingInfo
    {
        // ====================================================================
        // Paramètres de taille pour le buffer SoA
        // ====================================================================
        static constexpr CfemInt NUM_FIRST_ORDER_FIELDS  = 3; ///< \f$ \left\{\frac{\partial N}{\partial x}, \frac{\partial N}{\partial y}, \frac{\partial N}{\partial z}\right\} \f$
        static constexpr CfemInt NUM_SECOND_ORDER_FIELDS = 6; ///< \f$ \left\{\frac{\partial^2 N}{\partial x^2}, \frac{\partial^2 N}{\partial y^2}, \frac{\partial^2 N}{\partial z^2}, \frac{\partial^2 N}{\partial x \partial y}, d2N/dxdz, \frac{\partial^2 N}{\partial y \partial z} \right\} \f$

        // ====================================================================
        // Données d'état et Géométrie Locale (Small Tensors / Scalars)
        // ====================================================================
        CfemInt numNodes = 0;
        MappingEvaluationFlags validFlags = MappingEvaluationFlags::None;

        Vector3D physicalPosition;
        Tensor3D jacobian;
        Tensor3D invJacobian;
        CfemReal detJ = 0.0;
        Vector3D normal;

        CfemReal metric1D = 0.0;
        SymTensor2D metric2D;

        // ====================================================================
        // Stockage Contigu HPC pour les Dérivées Physiques Nodales
        // ====================================================================
        DynamicVector<CfemReal> storage; ///< Centralized contiguous buffer
        CfemInt stride = 0;              ///< Padded stride to ensure SIMD alignment
        bool hasSecondDerivatives = false;

        /**
         * @brief Allocates and configures the internal storage layout.
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

            validFlags = MappingEvaluationFlags::None;
        }

        inline void invalidate()
        {
            validFlags = MappingEvaluationFlags::None;
        }

        // ================================================================================================
        // SHAPES FIRST DERIVATIVES (Physical Space)
        // ================================================================================================

        [[nodiscard]] inline CfemReal* dN_dx() {
            return storage.data() + 0 * stride; 
        }

        [[nodiscard]] inline CfemReal* dN_dy() {
            return storage.data() + 1 * stride; 
        }

        [[nodiscard]] inline CfemReal* dN_dz() {
            return storage.data() + 2 * stride;
        }

        // Versions const
        [[nodiscard]] inline const CfemReal* dN_dx() const {
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "Physical first derivatives not computed?"); 
            return storage.data() + 0 * stride; 
        }

        [[nodiscard]] inline const CfemReal* dN_dy() const {
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "Physical first derivatives not computed?");
            return storage.data() + 1 * stride; 
        }

        [[nodiscard]] inline const CfemReal* dN_dz() const {
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "Physical first derivatives not computed?");
            return storage.data() + 2 * stride;
        }

        // ================================================================================================
        // SHAPES SECOND DERIVATIVES (Physical Space)
        // ================================================================================================

        [[nodiscard]] inline CfemReal* d2N_dx2() noexcept {
            return storage.data() + 3 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dy2() noexcept {
            return storage.data() + 4 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dz2() noexcept {
            return storage.data() + 5 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxdy() noexcept {
            return storage.data() + 6 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxdz() noexcept {
            return storage.data() + 7 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dydz() noexcept {
            return storage.data() + 8 * stride;
        }

        // Versions const
        [[nodiscard]] inline const CfemReal* d2N_dx2() const {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated?");
            return storage.data() + 3 * stride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dy2() const {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated?");
            return storage.data() + 4 * stride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dz2() const {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated?");
            return storage.data() + 5 * stride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dxdy() const {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated?");
            return storage.data() + 6 * stride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dxdz() const {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated?");
            return storage.data() + 7 * stride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dydz() const {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated?");
            return storage.data() + 8 * stride;
        }

        // ============================================================
        // Getters
        // ============================================================

        inline Vector3D gradient(CfemInt node) const
        {
            CFEM_ASSERT(node >= 0 && node < numNodes);
            return Vector3D(dN_dx()[node], dN_dy()[node], dN_dz()[node]);
        }

        inline SymTensor3D hessian(CfemInt node) const
        {
            CFEM_ASSERT(node >= 0 && node < numNodes);
            CFEM_ASSERT(hasSecondDerivatives);
            return SymTensor3D(d2N_dx2()[node], d2N_dy2()[node], d2N_dz2()[node],
                               d2N_dydz()[node], d2N_dxdz()[node], d2N_dxdy()[node]);
        }

        /**
         * @brief Checks if a specific evaluation field is marked as valid.
         */
        inline bool has(MappingEvaluationFlags field) const
        {
            return hasFlag(validFlags, field);
        }
    };

}