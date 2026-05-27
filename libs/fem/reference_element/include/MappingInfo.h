#pragma once

#include "cfemutils.h"
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
            stride = computePaddedStride(numNodes);

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

    /**
     * @struct MappingInfoBatch
     * @brief Stores metric data resulting from the geometric transformation for a batch of points.
     *
     * HPC Design: Nodal derivatives for ALL quadrature points and ALL fields are flattened 
     * into a single `storage` buffer. This guarantees perfect SIMD vectorization and eliminates 
     * multiple allocations and `DynamicMatrix` overhead.
     */
    struct MappingInfoBatch {

        static constexpr CfemInt NUM_FIRST_ORDER_FIELDS  = 3; ///< {dN_dx, dN_dy, dN_dz}
        static constexpr CfemInt NUM_SECOND_ORDER_FIELDS = 6; ///< {d2N_dx2, d2N_dy2, d2N_dz2, d2N_dxdy, d2N_dxdz, d2N_dydz}

        CfemInt numQuadraturePoints = 0;
        CfemInt numNodes = 0;
        MappingEvaluationFlags validFlags = MappingEvaluationFlags::None;
        bool isCellAffine = false;

        // --- Global quantities per quadrature point [qp] ---
        // Ces éléments restent dans leurs propres DynamicVector car ils sont 
        // par point de quadrature (1 scalaire/tenseur par QP), ils sont donc déjà contigus.
        DynamicVector<Vector3D> physicalPositions;
        DynamicVector<Tensor3D> jacobians;
        DynamicVector<Tensor3D> invJacobians;
        DynamicVector<CfemReal> detJs;
        DynamicVector<Vector3D> normals;

        DynamicVector<CfemReal> metric1Ds;
        DynamicVector<SymTensor2D> metric2Ds;

        // --- Stockage Contigu HPC pour les Dérivées Physiques Nodales ---
        DynamicVector<CfemReal> storage; ///< Centralized contiguous buffer for all QPs and all fields
        
        CfemInt paddedNodeStride = 0;    ///< Stride per field (padded to SIMD boundary)
        CfemInt qpStride = 0;            ///< Stride per quadrature point (all fields combined)
        bool hasSecondDerivatives = false;

        void setup(CfemInt nQp, CfemInt nNodes, CfemInt dim, bool includeSecondDerivatives = false) {
            numQuadraturePoints = nQp;
            numNodes = nNodes;
            hasSecondDerivatives = includeSecondDerivatives;

            // Resize des quantités globales par QP
            if (physicalPositions.size() != static_cast<size_t>(nQp)) {
                physicalPositions.resize(nQp);
                jacobians.resize(nQp);
                invJacobians.resize(nQp);
                detJs.resize(nQp);
                normals.resize(nQp);

                if (dim == 1){
                    metric1Ds.resize(nQp);
                }
                else if (dim == 2) {
                    metric2Ds.resize(nQp);
                }
            }

            // Calcul des strides et de la taille totale du buffer
            paddedNodeStride = computePaddedStride(numNodes);
            const CfemInt nFields = NUM_FIRST_ORDER_FIELDS + 
                                   (includeSecondDerivatives ? NUM_SECOND_ORDER_FIELDS : 0);
            
            qpStride = nFields * paddedNodeStride;
            const size_t totalStorageSize = static_cast<size_t>(nQp) * static_cast<size_t>(qpStride);

            if (storage.size() != totalStorageSize) {
                storage.resize(totalStorageSize);
            }

            validFlags = MappingEvaluationFlags::None;
        }

        inline void invalidate() {
            validFlags = MappingEvaluationFlags::None;
        }

        // ================================================================================================
        // VUES LÉGÈRES (Lightweight Views) : Renvoient un pointeur vers le début du champ pour le QP 'q'
        // ================================================================================================

        [[nodiscard]] inline CfemReal* dN_dx(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 0 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* dN_dy(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 1 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* dN_dz(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 2 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dx2(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 3 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dy2(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 4 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dz2(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 5 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxdy(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 6 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxdz(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 7 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dydz(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            return storage.data() + q * qpStride + 8 * paddedNodeStride;
        }

        // const versions

        [[nodiscard]] inline const CfemReal* dN_dx(CfemInt q) const{
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "First derivatives not computed");
            return storage.data() + q * qpStride + 0 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* dN_dy(CfemInt q) const{
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "First derivatives not computed");
            return storage.data() + q * qpStride + 1 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* dN_dz(CfemInt q) const{
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "First derivatives not computed");
            return storage.data() + q * qpStride + 2 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dx2(CfemInt q) const{
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), "Second derivatives not allocated");
            return storage.data() + q * qpStride + 3 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dy2(CfemInt q) const{
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), "Second derivatives not allocated");
            return storage.data() + q * qpStride + 4 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dz2(CfemInt q) const{
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), "Second derivatives not allocated");
            return storage.data() + q * qpStride + 5 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dxdy(CfemInt q) const{
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), "Second derivatives not allocated");
            return storage.data() + q * qpStride + 6 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dxdz(CfemInt q) const{
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), "Second derivatives not allocated");
            return storage.data() + q * qpStride + 7 * paddedNodeStride;
        }

        [[nodiscard]] inline const CfemReal* d2N_dydz(CfemInt q) const{
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), "Second derivatives not allocated");
            return storage.data() + q * qpStride + 8 * paddedNodeStride;
        }

        // ============================================================
        // Getters
        // ============================================================

        inline Vector3D gradient(CfemInt q, CfemInt node) const {
            CFEM_ASSERT(node >= 0 && node < numNodes, "Invalid node index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "Physical first derivative not computed?");

            return Vector3D( dN_dx(q)[node], dN_dy(q)[node], dN_dz(q)[node] );
        }

        inline SymTensor3D hessian(CfemInt q, CfemInt node) const {
            CFEM_ASSERT(node >= 0 && node < numNodes, "Invalid node index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), "Physical second derivative not computed?");
            return SymTensor3D(d2N_dx2(q)[node], d2N_dy2(q)[node], d2N_dz2(q)[node],
                               d2N_dydz(q)[node], d2N_dxdz(q)[node], d2N_dxdy(q)[node]);
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