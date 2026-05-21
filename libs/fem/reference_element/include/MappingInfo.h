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
     * @brief Stores metric data resulting from the geometric transformation (Reference -> Physical).
     */
    struct MappingInfo
    {

        CfemInt numNodes = 0;
        MappingEvaluationFlags validFlags = MappingEvaluationFlags::None;

        Vector3D physicalPosition;
        Tensor3D jacobian;
        Tensor3D invJacobian;
        CfemReal detJ = 0.0;
        Vector3D normal;

        DynamicVector<CfemReal> dN_dx; ///< First order derivatives with respect to \$ x \$
        DynamicVector<CfemReal> dN_dy; ///< First order derivatives with respect to \$ y \$
        DynamicVector<CfemReal> dN_dz; ///< First order derivatives with respect to \$ z \$

        DynamicVector<CfemReal> d2N_dx2; ///< Second order derivatives with respect to \$ x \$
        DynamicVector<CfemReal> d2N_dy2; ///< Second order derivatives with respect to \$ y \$
        DynamicVector<CfemReal> d2N_dz2; ///< Second order derivatives with respect to \$ z \$

        DynamicVector<CfemReal> d2N_dxdy; ///< Second order derivatives with respect to \$ x and y \$
        DynamicVector<CfemReal> d2N_dxdz; ///< Second order derivatives with respect to \$ x and z \$
        DynamicVector<CfemReal> d2N_dydz; ///< Second order derivatives with respect to \$ y and z \$

        CfemReal metric1D = 0.0;
        SymTensor2D metric2D;

        void setup(CfemInt nNodes, bool includeSecondDerivatives = false)
        {
            numNodes = nNodes;

            if (dN_dx.size() != static_cast<size_t>(numNodes))
            {
                dN_dx.resize(numNodes);
                dN_dy.resize(numNodes);
                dN_dz.resize(numNodes);
            }

            if (includeSecondDerivatives && d2N_dx2.size() != static_cast<size_t>(numNodes))
            {
                d2N_dx2.resize(numNodes);
                d2N_dy2.resize(numNodes);
                d2N_dz2.resize(numNodes);

                d2N_dxdy.resize(numNodes);
                d2N_dxdz.resize(numNodes);
                d2N_dydz.resize(numNodes);
            }

            validFlags = MappingEvaluationFlags::None;
        }

        inline Vector3D gradient(CfemInt node) const
        {
            CFEM_ASSERT(node >= 0 && node < numNodes);
            return Vector3D(dN_dx[node], dN_dy[node], dN_dz[node]);
        }

        inline SymTensor3D hessian(CfemInt node) const
        {
            // CFEM_INFO << "Hiiiiiiiiiiiiiiiiiiiiiii";
            // CFEM_INFO << "Size = " << d2N_dx2.size();
            CFEM_ASSERT(node >= 0 && node < numNodes);

            return SymTensor3D(d2N_dx2[node], d2N_dy2[node], d2N_dz2[node],
                               d2N_dydz[node], d2N_dxdz[node], d2N_dxdy[node]);
        }

        inline void invalidate()
        {
            validFlags = MappingEvaluationFlags::None;
        }

        void computePhysicalFirstDerivativesOnly(const ShapeInfo &fieldShape)
        {
            const CfemInt n = fieldShape.numNodes;
            CFEM_ASSERT(n == numNodes);

            const CfemReal i00 = invJacobian(0, 0);
            const CfemReal i01 = invJacobian(0, 1);
            const CfemReal i02 = invJacobian(0, 2);

            const CfemReal i10 = invJacobian(1, 0);
            const CfemReal i11 = invJacobian(1, 1);
            const CfemReal i12 = invJacobian(1, 2);

            const CfemReal i20 = invJacobian(2, 0);
            const CfemReal i21 = invJacobian(2, 1);
            const CfemReal i22 = invJacobian(2, 2);

            const CfemReal *dxi = fieldShape.dN_dxi.data();
            const CfemReal *deta = fieldShape.dN_deta.data();
            const CfemReal *dzeta = fieldShape.dN_dzeta.data();
            // CFEM_INFO << "inv jacobian";
            // CFEM_INFO << invJacobian;

            // CFEM_INFO << "Printing de dN_dzeta";
            // for (CfemInt l = 0; l < fieldShape.dN_dzeta.size(); ++l)
            //     CFEM_INFO << fieldShape.dN_dzeta[l] << " ";

            CfemReal *dx = dN_dx.data();
            CfemReal *dy = dN_dy.data();
            CfemReal *dz = dN_dz.data();

            for (CfemInt i = 0; i < n; ++i)
            {
                dx[i] = i00 * dxi[i] + i10 * deta[i] + i20 * dzeta[i];
                dy[i] = i01 * dxi[i] + i11 * deta[i] + i21 * dzeta[i];
                dz[i] = i02 * dxi[i] + i12 * deta[i] + i22 * dzeta[i];
            }
        }
    };

    /**
     * @struct MappingInfoBatch
     * @brief Stores metric data resulting from the geometric transformation for a batch of points.
     */
    struct MappingInfoBatch {

        CfemInt numQuadraturePoints = 0;
        CfemInt numNodes = 0;
        MappingEvaluationFlags validFlags = MappingEvaluationFlags::None;

        // --- Global quantities per quadrature point [qp] ---
        DynamicVector<Vector3D> physicalPositions;
        DynamicVector<Tensor3D> jacobians;
        DynamicVector<Tensor3D> invJacobians;
        DynamicVector<CfemReal> detJs;
        DynamicVector<Vector3D> normals;

        DynamicVector<CfemReal> metric1Ds;
        DynamicVector<SymTensor2D> metric2Ds;

        // --- Nodal quantities (SoA) [qp][node] ---
        DynamicMatrix<CfemReal> dN_dx;
        DynamicMatrix<CfemReal> dN_dy;
        DynamicMatrix<CfemReal> dN_dz;

        DynamicMatrix<CfemReal> d2N_dx2;
        DynamicMatrix<CfemReal> d2N_dy2;
        DynamicMatrix<CfemReal> d2N_dz2;
        DynamicMatrix<CfemReal> d2N_dxdy;
        DynamicMatrix<CfemReal> d2N_dxdz;
        DynamicMatrix<CfemReal> d2N_dydz;

        

        void setup(CfemInt nQp, CfemInt nNodes, CfemInt dim, bool includeSecondDerivatives = false) {
            numQuadraturePoints = nQp;
            numNodes = nNodes;

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

            if (dN_dx.size() != nQp * nNodes) {
                dN_dx.resize(nQp, nNodes);
                dN_dy.resize(nQp, nNodes);
                dN_dz.resize(nQp, nNodes);
            }

            if (includeSecondDerivatives && d2N_dx2.size() != nQp * nNodes) {
                d2N_dx2.resize(nQp, nNodes);
                d2N_dy2.resize(nQp, nNodes);
                d2N_dz2.resize(nQp, nNodes);

                d2N_dxdy.resize(nQp, nNodes);
                d2N_dxdz.resize(nQp, nNodes);
                d2N_dydz.resize(nQp, nNodes);
            }

            validFlags = MappingEvaluationFlags::None;
        }

        inline Vector3D gradient(CfemInt q, CfemInt node) const {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints);
            CFEM_ASSERT(node >= 0 && node < numNodes);
            return Vector3D(dN_dx(q, node), dN_dy(q, node), dN_dz(q, node));
        }

        inline SymTensor3D hessian(CfemInt q, CfemInt node) const {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints);
            CFEM_ASSERT(node >= 0 && node < numNodes);
            return SymTensor3D(d2N_dx2(q, node), d2N_dy2(q, node), d2N_dz2(q, node),
                               d2N_dydz(q, node), d2N_dxdz(q, node), d2N_dxdy(q, node));
        }

        inline void invalidate() {
            validFlags = MappingEvaluationFlags::None;
        }

        /**
         * @brief Computes physical gradients for the entire batch.
         * @note Hyper-optimized utilizing SIMD friendly contiguous memory rows.
         */
        void computePhysicalFirstDerivativesOnly(const ShapeInfoBatch& fieldShape) {
            CFEM_ASSERT(fieldShape.numQuadraturePoints == numQuadraturePoints);
            CFEM_ASSERT(fieldShape.numNodes == numNodes);

            for (CfemInt q = 0; q < numQuadraturePoints; ++q) 
            {
                // Cache invariant Jacobian terms in registers for this quadrature point
                const CfemReal i00 = invJacobians[q].data[Tensor3D::idx_xx], i01 = invJacobians[q].data[Tensor3D::idx_xy], i02 = invJacobians[q].data[Tensor3D::idx_xz];
                const CfemReal i10 = invJacobians[q].data[Tensor3D::idx_yx], i11 = invJacobians[q].data[Tensor3D::idx_yy], i12 = invJacobians[q].data[Tensor3D::idx_yz];
                const CfemReal i20 = invJacobians[q].data[Tensor3D::idx_zx], i21 = invJacobians[q].data[Tensor3D::idx_zy], i22 = invJacobians[q].data[Tensor3D::idx_zz];

                // Extract contiguous raw pointers for the current row (quadrature point 'q')
                // This utilizes the std::span returned by DynamicMatrix::operator[]
                const CfemReal* dxi   = fieldShape.dN_dxi[q].data();
                const CfemReal* deta  = fieldShape.dN_deta[q].data();
                const CfemReal* dzeta = fieldShape.dN_dzeta[q].data();

                CfemReal* dx = dN_dx[q].data();
                CfemReal* dy = dN_dy[q].data();
                CfemReal* dz = dN_dz[q].data();

                // Inner loop over nodes: perfectly contiguous, ripe for auto-vectorization
                for (CfemInt i = 0; i < numNodes; ++i) {
                    dx[i] = i00 * dxi[i] + i10 * deta[i] + i20 * dzeta[i];
                    dy[i] = i01 * dxi[i] + i11 * deta[i] + i21 * dzeta[i];
                    dz[i] = i02 * dxi[i] + i12 * deta[i] + i22 * dzeta[i];
                }
            }
        }
    };

}