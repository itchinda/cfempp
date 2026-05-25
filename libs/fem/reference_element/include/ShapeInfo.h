#pragma once

#include "cfemutils.h"

#include "EvaluationFlags.h"
#include "Tensor3D.h"
#include "DynamicVector.h"
#include "DynamicMatrix.h"

namespace cfem
{

    // ========================================================================
    // Paramètres matériels (Hardware-aware)
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
            stride = computePaddedStride(numNodes);

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
        [[nodiscard]] inline CfemReal* values()   {
            CFEM_ASSERT(has(ShapeEvaluationFlags::Value), "Shape values not computed?");
            return storage.data();
        }

        // ================================================================================================
        // SHAPES FIRST DERIVATIVES
        // ================================================================================================

        /** @brief Returns a writable pointer to shape derivatives along the reference \f$\xi\f$ direction. */
        [[nodiscard]] inline CfemReal* dN_dxi() {
            CFEM_ASSERT(has(ShapeEvaluationFlags::FirstDerivatives), "Shape first derivatives not computed?"); 
            return storage.data() + 1 * stride; 
        }

        /** @brief Returns a writable pointer to shape derivatives along the reference \f$\eta\f$ direction. */
        [[nodiscard]] inline CfemReal* dN_deta() {
            CFEM_ASSERT(has(ShapeEvaluationFlags::FirstDerivatives), "Shape first derivatives not computed?");
            return storage.data() + 2 * stride; 
        }

        /** @brief Returns a writable pointer to shape derivatives along the reference \f$\zeta\f$ direction. */
        [[nodiscard]] inline CfemReal* dN_dzeta() {
            CFEM_ASSERT(has(ShapeEvaluationFlags::FirstDerivatives), "Shape first derivatives not computed?");
            return storage.data() + 3 * stride;
        }

        // ================================================================================================
        // SHAPES SECOND DERIVATIVES
        // ================================================================================================

        /** @brief Returns a writable pointer to shape second derivatives \f$\partial^2 N \partial \zeta^2\f$ direction. */
        [[nodiscard]] inline CfemReal* d2N_dxi2()      {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivative not computed?");
            return storage.data() + 4 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_deta2()     {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivative not computed?");
            return storage.data() + 5 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dzeta2()    {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivative not computed?");
            return storage.data() + 6 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxideta()   {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivative not computed?");
            return storage.data() + 7 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxidzeta()  {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivative not computed?");
            return storage.data() + 8 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_detadzeta() {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivative not computed?");
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


    /**
     * @struct ShapeInfoBatch
     * @brief Contiguous, cache-aligned storage buffer for shape function evaluations across a batch of quadrature points.
     *
     * ### HPC Design (Structure of Arrays - SoA)
     * To maximize SIMD auto-vectorization and cache hits, all evaluations (values, gradients, Hessians) 
     * for all quadrature points and all nodes are flattened into a single 1D `DynamicVector`.
     * * Memory Layout for a given quadrature point `q`:
     * `[ N_nodes... | dxi_nodes... | deta_nodes... | dzeta_nodes... | (Hessians...) ]`
     * * Each field block is padded to ensure it starts on a SIMD-friendly memory boundary (e.g., 64 bytes).
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

        /** * @brief Prepares buffer sizes and strides to prevent runtime reallocations.
         * * @param nQp Number of quadrature points in the batch.
         * @param nNodes Number of interpolation nodes in the reference element.
         * @param includeSecondDerivatives If true, memory for Hessian-related quantities is allocated.
         */
        inline void setup(CfemInt nQp, CfemInt nNodes, bool includeSecondDerivatives = false)
        {
            numQuadraturePoints = nQp;
            numNodes = nNodes;
            hasSecondDerivatives = includeSecondDerivatives;

            paddedNodeStride = computePaddedStride(numNodes);
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
            CFEM_ASSERT(has(ShapeEvaluationFlags::Value), "Shape values not computed?");
            return storage.data() + q * qpStride + 0 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* dN_dxi(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(has(ShapeEvaluationFlags::FirstDerivatives), "Shape first derivatives not computed?");
            return storage.data() + q * qpStride + 1 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* dN_deta(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(has(ShapeEvaluationFlags::FirstDerivatives), "Shape first derivatives not computed?");
            return storage.data() + q * qpStride + 2 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* dN_dzeta(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(has(ShapeEvaluationFlags::FirstDerivatives), "Shape first derivatives not computed?");
            return storage.data() + q * qpStride + 3 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxi2(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated");
            return storage.data() + q * qpStride + 4 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_deta2(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated");
            return storage.data() + q * qpStride + 5 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dzeta2(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated");
            return storage.data() + q * qpStride + 6 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxideta(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated");
            return storage.data() + q * qpStride + 7 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxidzeta(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated");
            return storage.data() + q * qpStride + 8 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_detadzeta(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated");
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