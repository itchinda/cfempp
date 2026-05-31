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
#include "DynamicMatrix.h"
#include "DynamicMatrix.h"
#include "Tensor3D.h"
#include "Tensor2D.h"

#include "MappingInfoBatchView.h"

namespace cfem
{
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
            paddedNodeStride = kernel_utils::computePaddedStride(numNodes);
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

        [[nodiscard]] MappingInfoBatchView createView() const 
        {
            MappingInfoBatchView view;
            
            view.numQuadraturePoints = this->numQuadraturePoints;
            view.numNodes = this->numNodes;
            view.validFlags = this->validFlags;
            view.isCellAffine = this->isCellAffine;
            
            // std::span prend (pointeur_brut, taille). Coût de l'opération : ~16 octets copiés.
            view.physicalPositions = std::span<const Vector3D>(physicalPositions.data(), numQuadraturePoints);
            view.invJacobians  = std::span<const Tensor3D>(invJacobians.data(), numQuadraturePoints);
            view.detJs         = std::span<const CfemReal>(detJs.data(), numQuadraturePoints);
            
            if (!normals.empty()) {
                view.normals = std::span<const Vector3D>(normals.data(), numQuadraturePoints);
            }

            view.paddedNodeStride = this->paddedNodeStride;
            
            // Note: paddedNodeStride et les dN_dx_data seront remplis 
            // par l'intégrateur via l'Arena, pas par cette méthode de géométrie globale.
            
            return view;
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