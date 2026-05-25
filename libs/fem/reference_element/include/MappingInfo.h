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
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "Physical first derivatives not computed?"); 
            return storage.data() + 0 * stride; 
        }

        [[nodiscard]] inline CfemReal* dN_dy() {
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "Physical first derivatives not computed?");
            return storage.data() + 1 * stride; 
        }

        [[nodiscard]] inline CfemReal* dN_dz() {
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "Physical first derivatives not computed?");
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

        [[nodiscard]] inline CfemReal* d2N_dx2() {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated?");
            return storage.data() + 3 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dy2() {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated?");
            return storage.data() + 4 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dz2() {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated?");
            return storage.data() + 5 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxdy() {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated?");
            return storage.data() + 6 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxdz() {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated?");
            return storage.data() + 7 * stride;
        }

        [[nodiscard]] inline CfemReal* d2N_dydz() {
            CFEM_ASSERT(hasSecondDerivatives, "Second derivatives not allocated?");
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

        // ============================================================
        // Calcul (Optimisé HPC)
        // ============================================================

        void computePhysicalFirstDerivativesOnly(const ShapeInfo &fieldShape)
        {
            const CfemInt n = fieldShape.numNodes;
            CFEM_ASSERT(n == numNodes);

            // Hoisting (mise en cache des variables invariantes de la boucle)
            const CfemReal i00 = invJacobian.xx();
            const CfemReal i01 = invJacobian.xy();
            const CfemReal i02 = invJacobian.xz();

            const CfemReal i10 = invJacobian.yx();
            const CfemReal i11 = invJacobian.yy();
            const CfemReal i12 = invJacobian.yz();

            const CfemReal i20 = invJacobian.zx();
            const CfemReal i21 = invJacobian.zy();
            const CfemReal i22 = invJacobian.zz();

            // Pointeurs en lecture (restrict = promesse de non-aliasing au compilateur)
            const CfemReal* __restrict dxi   = fieldShape.dN_dxi();
            const CfemReal* __restrict deta  = fieldShape.dN_deta();
            const CfemReal* __restrict dzeta = fieldShape.dN_dzeta();

            // L'astuce : On contourne temporairement l'ASSERT de validFlags pour l'écriture
            // (car ils ne seront valides qu'à la fin de cette fonction)
            CfemReal* __restrict dx = storage.data() + 0 * stride;
            CfemReal* __restrict dy = storage.data() + 1 * stride;
            CfemReal* __restrict dz = storage.data() + 2 * stride;

            // Boucle principale vectorisable par le compilateur
            for (CfemInt i = 0; i < n; ++i)
            {
                dx[i] = i00 * dxi[i] + i10 * deta[i] + i20 * dzeta[i];
                dy[i] = i01 * dxi[i] + i11 * deta[i] + i21 * dzeta[i];
                dz[i] = i02 * dxi[i] + i12 * deta[i] + i22 * dzeta[i];
            }

            validFlags |= MappingEvaluationFlags::PhysicalFirstDerivatives;
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
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "Derivatives not computed");
            return storage.data() + q * qpStride + 0 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* dN_dy(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "Derivatives not computed");
            return storage.data() + q * qpStride + 1 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* dN_dz(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "Derivatives not computed");
            return storage.data() + q * qpStride + 2 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dx2(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), "Second derivatives not allocated");
            return storage.data() + q * qpStride + 3 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dy2(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), "Second derivatives not allocated");
            return storage.data() + q * qpStride + 4 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dz2(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), "Second derivatives not allocated");
            return storage.data() + q * qpStride + 5 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxdy(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), "Second derivatives not allocated");
            return storage.data() + q * qpStride + 6 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dxdz(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), "Second derivatives not allocated");
            return storage.data() + q * qpStride + 7 * paddedNodeStride;
        }

        [[nodiscard]] inline CfemReal* d2N_dydz(CfemInt q) {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints, "Invalid QP index");
            CFEM_ASSERT(hasFlag(validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), "Second derivatives not allocated");
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

        // ============================================================
        // Calcul (Optimisé HPC)
        // ============================================================

        /**
         * @brief Computes physical gradients for the entire batch.
         * @note Hyper-optimized utilizing SIMD friendly contiguous memory rows.
         */
        void computePhysicalFirstDerivativesOnly(const ShapeInfoBatch& fieldShape) {
            CFEM_ASSERT(fieldShape.numQuadraturePoints == numQuadraturePoints, "Mismatch QP");
            CFEM_ASSERT(fieldShape.numNodes == numNodes, "Mismatch Nodes");

            for (CfemInt q = 0; q < numQuadraturePoints; ++q) 
            {
                // Cache invariant Jacobian terms in registers for this quadrature point
                const CfemReal i00 = invJacobians[q].data[Tensor3D::idx_xx], i01 = invJacobians[q].data[Tensor3D::idx_xy], i02 = invJacobians[q].data[Tensor3D::idx_xz];
                const CfemReal i10 = invJacobians[q].data[Tensor3D::idx_yx], i11 = invJacobians[q].data[Tensor3D::idx_yy], i12 = invJacobians[q].data[Tensor3D::idx_yz];
                const CfemReal i20 = invJacobians[q].data[Tensor3D::idx_zx], i21 = invJacobians[q].data[Tensor3D::idx_zy], i22 = invJacobians[q].data[Tensor3D::idx_zz];

                // Remarque : Si vous avez appliqué la même logique à ShapeInfoBatch, 
                // ceci devrait être `fieldShape.dN_dxi(q);`. 
                // Si ShapeInfoBatch utilise encore l'ancienne méthode (DynamicMatrix), laissez `.data()`.
                const CfemReal* __restrict dxi   = fieldShape.dN_dxi(q);
                const CfemReal* __restrict deta  = fieldShape.dN_deta(q);
                const CfemReal* __restrict dzeta = fieldShape.dN_dzeta(q);

                // Extraction des pointeurs en écriture (on contourne les ASSERT temporairement)
                CfemReal* __restrict dx = storage.data() + q * qpStride + 0 * paddedNodeStride;
                CfemReal* __restrict dy = storage.data() + q * qpStride + 1 * paddedNodeStride;
                CfemReal* __restrict dz = storage.data() + q * qpStride + 2 * paddedNodeStride;

                // Inner loop over nodes: perfectly contiguous, ripe for auto-vectorization
                for (CfemInt i = 0; i < numNodes; ++i) {
                    dx[i] = i00 * dxi[i] + i10 * deta[i] + i20 * dzeta[i];
                    dy[i] = i01 * dxi[i] + i11 * deta[i] + i21 * dzeta[i];
                    dz[i] = i02 * dxi[i] + i12 * deta[i] + i22 * dzeta[i];
                }
            }

            // Marquer comme valide une fois la boucle de calcul terminée
            validFlags |= MappingEvaluationFlags::PhysicalFirstDerivatives;
        }
    };

}