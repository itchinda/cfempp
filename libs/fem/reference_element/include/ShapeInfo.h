#pragma once

#include "cfemutils.h"

#include "EvaluationFlags.h"
#include "Tensor3D.h"
#include "DynamicVector.h"
#include "DynamicMatrix.h"

namespace cfem
{

    /**
     * @struct ShapeInfo
     * @brief Storage buffer for shape function evaluations at a specific quadrature point.
     * @note Re-use this object across multiple integration points to avoid frequent heap allocations.
     */
    struct ShapeInfo
    {
        CfemInt numNodes = 0;
        bool shapeMustBeUpdated = true;
        ShapeEvaluationFlags validFlags = ShapeEvaluationFlags::None;

        DynamicVector<CfemReal> values; ///< Shape values (N_i)

        DynamicVector<CfemReal> dN_dxi;        ///< First order shape derivatives: \$ \frac{dN_i}{d\xi}  \$
        DynamicVector<CfemReal> dN_deta;       ///< First order shape derivatives: \$ \frac{dN_i}{d\eta}  \$
        DynamicVector<CfemReal> dN_dzeta;      ///< First order shape derivatives: \$ \frac{dN_i}{d\zeta}  \$

        DynamicVector<CfemReal> d2N_dxi2;      ///< Second order local shape derivatives: \$ \frac{\partial^2N_i}{\partial\xi^2}  \$
        DynamicVector<CfemReal> d2N_deta2;     ///< Second order local shape derivatives: \$ \frac{\partial^2N_i}{\partial\eta^2}  \$
        DynamicVector<CfemReal> d2N_dzeta2;    ///< Second order local shape derivatives: \$ \frac{\partial^2N_i}{\partial\zeta^2}  \$

        DynamicVector<CfemReal> d2N_dxideta;   ///< Second order local shape derivatives: \$ \frac{\partial^2N_i}{\partial\xi \partial\eta}  \$
        DynamicVector<CfemReal> d2N_dxidzeta;  ///< Second order local shape derivatives: \$ \frac{\partial^2N_i}{\partial\xi \partial\zeta}  \$
        DynamicVector<CfemReal> d2N_detadzeta; ///< Second order local shape derivatives: \$ \frac{\partial^2N_i}{\partial\eta \partial\zeta}  \$

        /**
         * @brief Prepares buffer sizes to prevent runtime reallocations.
         * @param nNodes Number of nodes in the element.
         * @param includeSecondDerivatives If true, allocates memory for the second order derivatives
         */
        void setup(CfemInt nNodes, bool includeSecondDerivatives = false)
        {
            numNodes = nNodes;

            if (values.size() != static_cast<size_t>(numNodes)) {
                values.resize(numNodes);
                dN_dxi.resize(numNodes);
                dN_deta.resize(numNodes);
                dN_dzeta.resize(numNodes);
            }

            if (includeSecondDerivatives && d2N_dxi2.size() != static_cast<size_t>(numNodes)){
                d2N_dxi2.resize(numNodes);
                d2N_deta2.resize(numNodes);
                d2N_dzeta2.resize(numNodes);

                d2N_dxideta.resize(numNodes);
                d2N_dxidzeta.resize(numNodes);
                d2N_detadzeta.resize(numNodes);
            }

            validFlags = ShapeEvaluationFlags::None;
        }

        void resizeFirstDerivativeBuffers(CfemInt nNodes)
        {
            numNodes = nNodes;
            dN_dxi.resize(numNodes);
            dN_deta.resize(numNodes);
            dN_dzeta.resize(numNodes);

        }

        void resizeSecondDerivativeBuffers(CfemInt nNodes)
        {
            numNodes = nNodes;
            d2N_dxi2.resize(numNodes);
            d2N_deta2.resize(numNodes);
            d2N_dzeta2.resize(numNodes);

            d2N_dxideta.resize(numNodes);
            d2N_dxidzeta.resize(numNodes);
            d2N_detadzeta.resize(numNodes);

        }

        inline Vector3D gradient(CfemInt node) const
        {
            CFEM_ASSERT(node >= 0 && node < numNodes);
            return Vector3D(dN_dxi[node], dN_deta[node], dN_dzeta[node]);
        }

        inline SymTensor3D hessian(CfemInt node) const
        {
            CFEM_ASSERT(node >= 0 && node < numNodes);
            return SymTensor3D(d2N_dxi2[node], d2N_deta2[node], d2N_dzeta2[node], // The diagonal entries
                                d2N_detadzeta[node], d2N_dxidzeta[node], d2N_dxideta[node]);
        }

        /**
         * @brief Checks if a specific field (Value/Gradient/Hessian) is marked as valid.
         */
        bool has(ShapeEvaluationFlags field) const
        {
            return hasFlag(validFlags, field);
        }
    };

    /**
     * @struct ShapeInfoBatch
     * @brief Storage buffer for shape function evaluations at a collection of quadrature points.
     * @note Optimized for HPC using a Structure of Arrays (SoA) layout via DynamicMatrix.
     */
    struct ShapeInfoBatch
    {
        CfemInt numQuadraturePoints = 0;
        CfemInt numNodes = 0;
        bool shapeMustBeUpdated = true;
        ShapeEvaluationFlags validFlags = ShapeEvaluationFlags::None;

        DynamicMatrix<CfemReal> values;         // Shape values: values[qp][i]

        DynamicMatrix<CfemReal> dN_dxi;         // First order local shapes' derivatives at a collection of quadrature points: \$ \frac{\partial N_i}{\partial\xi}[q] \$
        DynamicMatrix<CfemReal> dN_deta;        // First order local shapes' derivatives at a collection of quadrature points: \$ \frac{\partial N_i}{\partial\eta}[q] \$
        DynamicMatrix<CfemReal> dN_dzeta;       // First order local shapes' derivatives at a collection of quadrature points: \$ \frac{\partial N_i}{\partial\xi}[q] \$
        
        DynamicMatrix<CfemReal> d2N_dxi2;       // Second order local shapes' derivatives at a collection of quadrature points \$ \frac{\partial^2N_i}{\partial\xi^2}[q] \$
        DynamicMatrix<CfemReal> d2N_deta2;      // Second order local shapes' derivatives at a collection of quadrature points \$ \frac{\partial^2N_i}{\partial\eta^2}[q] \$
        DynamicMatrix<CfemReal> d2N_dzeta2;     // Second order local shapes' derivatives at a collection of quadrature points \$ \frac{\partial^2N_i}{\partial\zeta}[q] \$

        DynamicMatrix<CfemReal> d2N_dxideta;    // Second order local shapes' derivatives at a collection of quadrature points \$ \frac{\partial^2N_i}{\partial\xi\partial\eta}[q] \$
        DynamicMatrix<CfemReal> d2N_dxidzeta;   // Second order local shapes' derivatives at a collection of quadrature points \$ \frac{\partial^2N_i}{\partial\xi\partial\zeta}[q] \$
        DynamicMatrix<CfemReal> d2N_detadzeta;  // Second order local shapes' derivatives at a collection of quadrature points \$ \frac{\partial^2N_i}{\partial\eta\partial\zeta}[q] \$

        /** 
         * @brief Prepares buffer sizes to prevent runtime reallocations.
         */
        void setup(CfemInt nQp, CfemInt nNodes, bool includeSecondDerivatives = false)
        {
            numQuadraturePoints = nQp;
            numNodes = nNodes;

            if (values.size() != nQp * nNodes)
            {
                values.resize(nQp, nNodes);
                dN_dxi.resize(nQp, nNodes);
                dN_deta.resize(nQp, nNodes);
                dN_dzeta.resize(nQp, nNodes);
            }

            if (includeSecondDerivatives && d2N_dxi2.size() != nQp * nNodes)
            {
                d2N_dxi2.resize(nQp, nNodes);
                d2N_deta2.resize(nQp, nNodes);
                d2N_dzeta2.resize(nQp, nNodes);

                d2N_dxideta.resize(nQp, nNodes);
                d2N_dxidzeta.resize(nQp, nNodes);
                d2N_detadzeta.resize(nQp, nNodes);
            }

            validFlags = ShapeEvaluationFlags::None;
        }

        /** 
         * @brief Reconstructs the reference gradient vector for a specific qp and node. 
         */
        inline Vector3D gradient(CfemInt q, CfemInt i) const
        {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints);
            CFEM_ASSERT(i >= 0 && i < numNodes);
            return Vector3D(dN_dxi(q, i), dN_deta(q, i), dN_dzeta(q, i));
        }

        /** 
         * @brief Reconstructs the reference Hessian tensor for a specific qp and node. 
         */
        inline SymTensor3D hessian(CfemInt q, CfemInt i) const
        {
            CFEM_ASSERT(q >= 0 && q < numQuadraturePoints);
            CFEM_ASSERT(i >= 0 && i < numNodes);
            return SymTensor3D(d2N_dxi2(q, i), d2N_deta2(q, i), d2N_dzeta2(q, i), 
                               d2N_detadzeta(q, i), d2N_dxidzeta(q, i), d2N_dxideta(q, i));
        }

        bool has(ShapeEvaluationFlags field) const
        {
            return hasFlag(validFlags, field);
        }
    };

}