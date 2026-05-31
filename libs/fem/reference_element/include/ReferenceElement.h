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
#include "DynamicMatrix.h"
#include "DynamicVector.h"
#include "EvaluationFlags.h"
#include "MappingInfo.h"
#include "MappingInfoBatch.h"
#include "MappingInfoBatchView.h"
#include "ShapeInfo.h"
#include "ShapeInfoBatch.h"
#include "Tensor3D.h"
#include "Vector3D.h"

namespace cfem::reference_elem
{
    static const DynamicVector<DynamicVector<int>> line_face_nodes = {
        {0}, // Face 0 (Left node)
        {1}  // Face 1 (Right node)
    };
}

namespace cfem
{
    enum class Order
    {
        P0 = 0,
        P1 = 1,
        P2 = 2,
        Mini = 101
    };

    /**
     * @class ReferenceElement
     * @brief Abstract base class for Finite Element reference definitions.
     * 
     * Provides the machinery to compute shape functions, geometric mappings, Jacobians, 
     * and physical derivatives. It supports both traditional single-point evaluation 
     * and high-performance, SIMD-friendly batch evaluation using Structure-of-Arrays (SoA).
     */
    class ReferenceElement
    {
    protected:
        DynamicVector<Vector3D> m_nodes; ///< Node coordinates in reference space (\f$ \xi \f$).

        CfemInt m_numNodes;   ///< Reference element nodes count
        CfemInt m_dimension;  ///< Topological dimension of the reference element

    public:
        virtual ~ReferenceElement() = default;

        /**
         * @brief Returns the topological dimension of the reference element.
         * @return 1 for lines, 2 for triangles/quads, 3 for tetrahedrons/hexahedrons.
         */
        CfemInt getDimension() const { return m_dimension; };

        /**
         * @brief Returns the total number of nodes defining the reference element.
         */
        CfemInt getNumNodes() const { return m_numNodes; };

        /**
         * @brief Retrieves the local node indices that belong to a specific face of the element.
         * @param iFace The local index of the face.
         * @return A view (span) of the local node indices making up the face.
         */
        virtual std::span<const CfemInt> getFaceNodes(CfemInt iFace) const = 0;

        /**
         * @brief Centralized method to evaluate shape functions and/or their reference derivatives at a single point.
         * @param xi Target point in reference space.
         * @param requestedShape Bitmask indicating which fields to compute (Values, FirstDerivatives, and/or SecondDerivatives.).
         * @param out The ShapeInfo structure to populate.
         */
        void evaluateShapes(const Vector3D &xi, ShapeEvaluationFlags requestedShape, ShapeInfo &out) const;

        /**
         * @brief Centralized method for SIMD-friendly batch evaluation of shape functions and reference derivatives.
         * @param xiPoints Array of target points in reference space.
         * @param requestedShape Bitmask indicating which fields to compute.
         * @param outBatch The ShapeInfoBatch structure (SoA layout) to populate.
         */
        void evaluateShapesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeEvaluationFlags requestedShape, ShapeInfoBatch &outBatch) const;

        /**
         * @brief Evaluates the shape function values (\f$ N_i \f$) at a given reference point.
         * @param xi The coordinates of the evaluation point in the reference element (\f$\xi, \eta, \zeta \f$).
         * @param info ShapeInfo object whose value buffer will be populated.
         */
        virtual void evaluateShapesValues(const Vector3D &xi, ShapeInfo &info) const = 0;

        /**
         * @brief Vectorized evaluation of shape function values (\f$ N_i \f$) for a batch of reference points.
         * @param xi Array of evaluation points.
         * @param shapeBatch SoA structure to store the computed values.
         */
        virtual void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xi,
                                               ShapeInfoBatch &shapeBatch) const = 0;

        /**
         * @brief Evaluates the first-order derivatives of the shape functions (local gradients) 
         * with respect to the reference coordinate system.
         * @param xi The coordinates of the evaluation point in the reference element (\f$\xi, \eta, \zeta \f$).
         * @param info ShapeInfo object whose derivative buffers will be populated.
         */
        virtual void evaluateShapesDerivatives(const Vector3D &xi, ShapeInfo &info) const = 0;

        /**
         * @brief Vectorized evaluation of first-order reference derivatives for a batch of points.
         * @param xi Array of evaluation points.
         * @param shapeBatch SoA structure to store the computed reference derivatives.
         */
        virtual void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xi,
                                                    ShapeInfoBatch &shapeBatch) const = 0;

        /**
         * @brief Evaluates the second-order derivatives of the shape functions (local Hessians) 
         * with respect to the reference coordinate system.
         * @param xi The coordinates of the evaluation point.
         * @param info The ShapeInfo object containing buffers to be filled.
         */
        virtual void evaluateShapesSecondDerivatives(const Vector3D &xi, ShapeInfo &info) const = 0;

        /**
         * @brief Vectorized evaluation of second-order reference derivatives for a batch of points.
         * @param xi Array of evaluation points.
         * @param shapeBatch SoA structure to store the computed reference Hessians.
         */
        virtual void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xi,
                                                          ShapeInfoBatch &shapeBatch) const = 0;

        /**
         * @brief Projects a reference coordinate from a lower-dimensional boundary facet into the parent cell's reference space.
         * @note Crucial for boundary integrations where the surface measure is computed on the facet, but 
         * physical derivatives require shape evaluations inside the parent volume.
         * @param localFacetId The local index of the facet within the parent element.
         * @param xiFacet The reference coordinate on the facet.
         * @return The projected coordinate \f$ (\xi, \eta, \zeta) \f$ in the parent element's reference space.
         */
        virtual Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const = 0;

        /**
         * @brief Returns the coordinates of the nodes defining the element in the reference space.
         */
        const DynamicVector<Vector3D> &getReferenceNodes() const { return m_nodes; }

        /**
         * @brief Returns the reference coordinate of a specific local node.
         * @param i Local node index.
         */
        const Vector3D &getNodeCoords(CfemInt i) const
        {
            CFEM_ASSERT(i < m_nodes.size() && i >= 0);
            return m_nodes[i];
        }

        /**
         * @brief Computes the full geometric mapping at a given point xi, internally evaluating required shape functions.
         * @param xi Target point in reference space.
         * @param physicalNodesCoords Actual physical coordinates of the element's nodes.
         * @param requestedShape Bitmask of reference shape fields required.
         * @param requestedMapping Bitmask of physical fields required (Jacobian, Gradients, Hessians).
         * @param shape Work buffer for reference shape functions.
         * @param out MappingInfo structure to be populated.
         */
        void computeMapping(const Vector3D &xi,
                            const DynamicVector<Vector3D> &physicalNodesCoords,
                            ShapeEvaluationFlags requestedShape,
                            MappingEvaluationFlags requestedMapping,
                            ShapeInfo &shape,
                            MappingInfo &out) const;

        /**
         * @brief Computes the full geometric mapping utilizing pre-evaluated reference shape functions.
         * @param xi Target point in reference space.
         * @param physicalNodesCoords Actual physical coordinates of the element's nodes.
         * @param requestedMapping Bitmask of physical fields required.
         * @param shape Pre-evaluated shape functions.
         * @param out MappingInfo structure to be populated.
         */
        void computeMapping(const Vector3D &xi,
                            const DynamicVector<Vector3D> &physicalNodesCoords,
                            MappingEvaluationFlags requestedMapping,
                            const ShapeInfo &shape,
                            MappingInfo &out) const;

        /**
         * @brief Computes the full geometric mapping for a batch of points, internally evaluating required shape functions.
         * @note Features an affine-bypass optimization: if the element is affine, Jacobians are evaluated only once.
         * @param xiPoints Array of target points in reference space.
         * @param physicalNodesCoords Actual physical coordinates of the element's nodes.
         * @param requestedShape Bitmask of reference shape fields required.
         * @param requestedMapping Bitmask of geometric fields required.
         * @param shapeBatch SoA work buffer for reference shape functions.
         * @param outBatch SoA MappingInfoBatch structure to be populated.
         */
        void computeMappingBatch(const DynamicVector<Vector3D> &xiPoints,
                                 const DynamicVector<Vector3D> &physicalNodesCoords,
                                 ShapeEvaluationFlags requestedShape,
                                 MappingEvaluationFlags requestedMapping,
                                 ShapeInfoBatch &shapeBatch,
                                 MappingInfoBatch &outBatch) const;

        /**
         * @brief Computes the geometric mapping for a batch of points utilizing pre-evaluated reference shapes.
         * @param xiPoints Array of target points in reference space.
         * @param physicalNodesCoords Actual physical coordinates of the element's nodes.
         * @param requestedMapping Bitmask of geometric fields required.
         * @param shapeBatch Pre-evaluated SoA reference shapes.
         * @param outBatch SoA MappingInfoBatch structure to be populated.
         */
        void computeMappingBatch(const DynamicVector<Vector3D> &xiPoints,
                                 const DynamicVector<Vector3D> &physicalNodesCoords,
                                 MappingEvaluationFlags requestedMapping,
                                 const ShapeInfoBatch &shapeBatch,
                                 MappingInfoBatch &outBatch) const;

        /**
         * @brief Transforms reference gradients to physical space using \f$ J^{-T} \f$.
         * @param shape Pre-evaluated reference shape functions.
         * @param out MappingInfo structure (must contain valid inverse Jacobian) to be populated with physical derivatives.
         */
        void computePhysicalFirstDerivatives(const ShapeInfo &shape,
                                             MappingInfo &out) const;

        /**
         * @brief High-performance SIMD kernel transforming reference gradients to physical space for a batch of points.
         * @note Utilizes dimensional masking to avoid unnecessary computations in 2D or 1D spaces, and safely handles 
         * embedded manifolds (e.g., computing 3D physical gradients for a 2D surface element). Memory is managed via the View.
         * @param shape Pre-evaluated reference shape batch.
         * @param out The Zero-Copy Hybrid View containing the target Arena memory layout.
         */
        void computePhysicalFirstDerivativesBatch(const ShapeInfoBatch &shape,
                                                  MappingInfoBatchView &out) const;

        /**
         * @brief Computes physical Hessians including the geometric curvature correction term.
         * @note Reference: \f$ H_{phys} = J^{-T} * [H_{ref} - \nabla_{phys}(N) * H_{geom}] * J^{-1} \f$
         * @param physicalNodesCoords Actual physical coordinates of the element's nodes.
         * @param shape Pre-evaluated reference shape functions.
         * @param out MappingInfo structure to be populated.
         */
        void computePhysicalSecondDerivatives(const DynamicVector<Vector3D> &physicalNodesCoords,
                                              const ShapeInfo &shape,
                                              MappingInfo &out) const;

        /**
         * @brief Computes physical Hessians for a batch of points, accounting for geometric curvature.
         * @param physicalNodesCoords Actual physical coordinates of the element's nodes.
         * @param shapeBatch Pre-evaluated reference shape batch.
         * @param outBatch SoA MappingInfoBatch to be populated with physical Hessians.
         */
        void computePhysicalSecondDerivativesBatch(const DynamicVector<Vector3D> &physicalNodesCoords,
                                                   const ShapeInfoBatch &shapeBatch,
                                                   MappingInfoBatch &outBatch) const;

        /**
         * @brief Interpolates a field (scalar or vector) at a point evaluated in ShapeInfo.
         * @tparam T The field type (e.g., CfemReal, Vector3D, or Tensor types).
         * @param shape Pre-evaluated shape functions at point xi.
         * @param nodalValues The values of the field at the element nodes.
         * @return The interpolated value.
         */
        template <typename T>
        T interpolate(const ShapeInfo &shape, const DynamicVector<T> &nodalValues) const;

        /**
         * @brief Interpolates the physical gradient of a scalar field at a point.
         * @param map MappingInfo containing physical first derivatives of the shape functions.
         * @param nodalValues The values of the field at the element nodes.
         * @return Vector3D containing \f$ [\frac{\partial u}{\partial x}, \frac{\partial u}{\partial y}, \frac{\partial u}{\partial z}] \f$.
         */
        Vector3D interpolateGradient(const MappingInfo &map,
                                     const DynamicVector<CfemReal> &nodalValues) const;

        /**
         * @brief Interpolates the physical Hessian of a scalar field.
         * @param map MappingInfo containing physical second derivatives of the shape functions.
         * @param nodalValues The values of the field at the element nodes.
         * @return SymTensor3D representing the symmetric 2nd derivative matrix.
         */
        SymTensor3D interpolateHessian(const MappingInfo &map,
                                       const DynamicVector<CfemReal> &nodalValues) const;

    private:
        /**
         * @brief Computes the Jacobian, its determinant (metric), and the pseudo-inverse for a single point.
         * @note This version safely supports 1D/2D elements embedded in 3D space using 
         * Moore-Penrose pseudo-inverse logic.
         * @param shape Pre-evaluated reference shape functions.
         * @param physicalNodesCoords Actual physical coordinates of the element's nodes.
         * @param out MappingInfo structure to be populated with the Jacobian.
         */
        inline void computeJacobian(const ShapeInfo &shape,
                                    const DynamicVector<Vector3D> &physicalNodesCoords,
                                    MappingInfo &out) const;

        /**
         * @brief Computes the Jacobian matrices for a batch of quadrature points.
         * @param shape SoA batch containing pre-evaluated reference shape derivatives.
         * @param physicalNodesCoords Actual physical coordinates of the element's nodes.
         * @param out MappingInfoBatch structure to store the computed Jacobians.
         */
        inline void computeJacobianBatch(const ShapeInfoBatch &shape,
                                         const DynamicVector<Vector3D> &physicalNodesCoords,
                                         MappingInfoBatch &out) const;

        /**
         * @brief Computes the determinant of the Jacobian (metric / volume or surface element) for a single point.
         * @note For embedded manifolds (e.g., 2D surface in 3D space), this evaluates the Gram determinant metric.
         * @param out MappingInfo structure containing a valid Jacobian, to be populated with its determinant.
         */
        inline void computeDeterminant(MappingInfo &out) const;

        /**
         * @brief SIMD-friendly computation of Jacobian determinants for a batch of points.
         * @param out SoA MappingInfoBatch structure containing valid Jacobians, to be populated with their determinants.
         */
        inline void computeDeterminantBatch(MappingInfoBatch &out) const;

        /**
         * @brief Computes the inverse of the Jacobian matrix for a single point.
         * @note Utilizes Moore-Penrose pseudo-inversion when the reference dimension is strictly less 
         * than the spatial dimension (e.g., computing the inverse for a surface element).
         * @param out MappingInfo structure containing a valid Jacobian, to be populated with its inverse.
         */
        inline void computeInverseJacobian(MappingInfo &out) const;

        /**
         * @brief SIMD-friendly computation of inverse Jacobian matrices for a batch of points.
         * @param out SoA MappingInfoBatch structure to be populated with the inverses.
         */
        inline void computeInverseJacobianBatch(MappingInfoBatch &out) const;

        /**
         * @brief Computes the physical geometric normal vector at a single evaluation point.
         * @note Primarily applicable for \f$ (N-1) \f$-dimensional boundary or interface elements (facets).
         * @param out MappingInfo structure containing a valid Jacobian, to be populated with the normal vector.
         */
        inline void computeGeometricNormal(MappingInfo &out) const;

        /**
         * @brief SIMD-friendly computation of physical geometric normal vectors for a batch of points.
         * @param outBatch SoA MappingInfoBatch structure to be populated with the normal vectors.
         */
        inline void computeGeometricNormalBatch(MappingInfoBatch &outBatch) const;

        /**
         * @brief Legacy batch derivative computation (Replaced by MappingInfoBatchView architecture).
         */
        void computePhysicalFirstDerivativesBatch(const ShapeInfoBatch &shape,
                                                  MappingInfoBatch &out) const;
    };

}

#include "ReferenceElement.tpp"
