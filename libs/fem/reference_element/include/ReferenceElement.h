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
#include "Vector3D.h"
#include "DynamicVector.h"
#include "DynamicMatrix.h"
#include "Tensor3D.h"
#include "EvaluationFlags.h"

#include "ShapeInfo.h"
#include "ShapeInfoBatch.h"
#include "MappingInfo.h"
#include "MappingInfoBatch.h"
#include "MappingInfoBatchView.h"


namespace cfem::reference_elem
{
    static const DynamicVector<DynamicVector<int>> line_face_nodes = {
        {0}, // Face 0 (Left node)
        {1}  // Face 1 (Right node)
    };
}

namespace cfem {
    enum class Order { P0 = 0, P1 = 1, P2 = 2 , Mini = 101};

/**
 * @class ReferenceElement
 * @brief Abstract base class for Finite Element reference definitions.
 * Provides the machinery to compute geometric mappings, Jacobians, and 
 * physical derivatives from local shape functions.
 */
class ReferenceElement {
protected:
    DynamicVector<Vector3D> m_nodes; ///< Node coordinates in reference space (xi)

    CfemInt m_numNodes;
    CfemInt m_dimension;


public:
    virtual ~ReferenceElement() = default;

    CfemInt dimension() const { return m_dimension; };
    CfemInt getNumNodes() const { return m_numNodes; };
    virtual std::span<const CfemInt> getFaceNodes(CfemInt iFace) const = 0;

    void evaluateShapes(const Vector3D &xi, ShapeEvaluationFlags requestedShape, ShapeInfo& out) const;
    void evaluateShapesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeEvaluationFlags requestedShape, ShapeInfoBatch& outBatch) const;
    /**
     * @brief Evaluates the shape function values (\f$ N_i \f$) at a given reference point.
     * @param xi The coordinates of the evaluation point in the reference element (\f$\xi, \eta, \zeta \f$).
     * @param values Raw pointer to a pre-allocated array of size [numNodes] to be populated.
     */
    virtual void evaluateShapesValues(const Vector3D& xi, ShapeInfo& info) const = 0;
    virtual void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xi, 
                                           ShapeInfoBatch& shapeBatch) const = 0;

    /**
     * @brief Evaluates the first-order derivatives of the shape functions (local gradients) 
     * with respect to the reference coordinate system.
     * @param xi The coordinates of the evaluation point in the reference element (\f$\xi, \eta, \zeta \f$).
     * @param dN_dxi Raw pointer to the reference coordinate \f$ \xi \f$ derivative array (size [numNodes]).
     * @param dN_deta Raw pointer to the reference coordinate \f$ \eta \f$ derivative array (size [numNodes]).
     * @param dN_dzeta Raw pointer to the reference coordinate \f$ \zeta \f$ derivative array (size [numNodes]).
     */
    virtual void evaluateShapesDerivatives(const Vector3D& xi, ShapeInfo& info) const = 0;
    virtual void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xi, 
                                                ShapeInfoBatch& shapeBatch) const = 0;

    /**
     * @brief Evaluates the second-order derivatives of the shape functions (local Hessians) 
     * with respect to the reference coordinate system.
     * 
     * @param info The ShapeInfo object containing buffers to be filled.
     */
    virtual void evaluateShapesSecondDerivatives(const Vector3D& xi, ShapeInfo& info) const = 0;
    virtual void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xi,
                                                      ShapeInfoBatch& shapeBatch) const = 0;
                          
    virtual Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const = 0;
    /**
     * @brief Returns the coordinates of the nodes in the reference space. 
     */
    const DynamicVector<Vector3D>& referenceNodes() const { return m_nodes; }
    const Vector3D& getNodeCoords( CfemInt i) const { CFEM_ASSERT( i < m_nodes.size() && i >= 0 ); return m_nodes[i]; }

    

    /**
     * @brief Computes the full geometric mapping at a given point xi.
     * @param xi Target poCfemInt in reference space.
     * @param physicalNodesCoords Actual physical coordinates of the element's nodes.
     * @param requestedFields Bitmask of physical fields required (Gradients, Hessians).
     * @param shape Work buffer for reference shape functions.
     * @param out MappingInfo structure to be populated.
     */
    void computeMapping(const Vector3D                &xi,
                        const DynamicVector<Vector3D> &physicalNodesCoords,
                        ShapeEvaluationFlags           requestedShape,
                        MappingEvaluationFlags         requestedMapping,
                        ShapeInfo                     &shape,
                        MappingInfo                   &out) const;

    void computeMapping(const Vector3D                &xi,
                        const DynamicVector<Vector3D> &physicalNodesCoords,
                        MappingEvaluationFlags         requestedMapping,
                        const ShapeInfo               &shape,
                        MappingInfo                   &out) const;

    void computeMappingBatch(const DynamicVector<Vector3D>   &xiPoints,
                        const DynamicVector<Vector3D> &physicalNodesCoords,
                        ShapeEvaluationFlags           requestedShape,
                        MappingEvaluationFlags         requestedMapping,
                        ShapeInfoBatch                &shapeBatch,
                        MappingInfoBatch              &outBatch) const;

    void computeMappingBatch(const DynamicVector<Vector3D>   &xiPoints,
                        const DynamicVector<Vector3D> &physicalNodesCoords,
                        MappingEvaluationFlags         requestedMapping,
                        const ShapeInfoBatch          &shapeBatch,
                        MappingInfoBatch              &outBatch) const;

    /** 
     * @brief Transforms reference gradients to physical space using J^-T. 
     * @param shape Work buffer for reference shape functions.
     * @param out MappingInfo structure to be populated.
     */
    void computePhysicalFirstDerivatives(const ShapeInfo& shape, 
                                         MappingInfo& out) const;

    void computePhysicalFirstDerivativesBatch(const ShapeInfoBatch &shape,
                                              MappingInfoBatchView &out) const;

    /**
     * @brief Computes physical Hessians including the geometric curvature term.
     * @param physicalNodesCoords Actual physical coordinates of the element's nodes.
     * @param shape Work buffer for reference shape functions.
     * @param out MappingInfo structure to be populated.
     * @note Reference: \f$ H_{phys} = J^{-T} * [H_{ref} - \nabla_{phys}(N) * H_{geom}] * J^{-1} \f$
     */
    void computePhysicalSecondDerivatives(const DynamicVector<Vector3D>& physicalNodesCoords, 
                                          const ShapeInfo& shape, 
                                          MappingInfo& out) const;

    void computePhysicalSecondDerivativesBatch(const DynamicVector<Vector3D>& physicalNodesCoords, 
                                               const ShapeInfoBatch& shapeBatch, 
                                               MappingInfoBatch& outBatch) const;

    /**
     * @brief Interpolates a field (scalar or vector) at a point evaluated in ShapeInfo.
     * @tparam T The field type (e.g., CfemReal, Vector3D, or even Tensor types).
     * @param shape Pre-evaluated shape functions at point xi.
     * @param nodalValues The values of the field at the element nodes.
     */
    template<typename T>
    T interpolate(const ShapeInfo& shape, const DynamicVector<T>& nodalValues) const;

    /**
     * @brief Interpolates the physical gradient of a scalar field at a point.
     * @return Vector3D containing [du/dx, du/dy, du/dz]
     */
    Vector3D interpolateGradient(const MappingInfo& map, 
                                 const DynamicVector<CfemReal>& nodalValues) const;
    

    /**
     * @brief Interpolates the physical Hessian of a scalar field.
     * @return SymTensor3D representing the 2nd derivative matrix.
     */
    SymTensor3D interpolateHessian(const MappingInfo& map, 
                                   const DynamicVector<CfemReal>& nodalValues) const;

private:
            
    /**
     * @brief Computes the Jacobian, its determinant (metric), and the pseudo-inverse.
     * @note This version supports 1D/2D elements embedded in 3D space using 
     * the Moore-Penrose pseudo-inverse logic.
     */
    inline void computeJacobian(const ShapeInfo& shape, 
                                const DynamicVector<Vector3D>& physicalNodesCoords, 
                                MappingInfo& out) const;
    inline void computeJacobianBatch(const ShapeInfoBatch& shape, 
                                const DynamicVector<Vector3D>& physicalNodesCoords, 
                                MappingInfoBatch& out) const;

    inline void computeDeterminant(MappingInfo &out) const;
    inline void computeDeterminantBatch(MappingInfoBatch &out) const;

    inline void computeInverseJacobian(MappingInfo &out) const;
    inline void computeInverseJacobianBatch(MappingInfoBatch &out) const;

    inline void computeGeometricNormal(MappingInfo &out) const;
    inline void computeGeometricNormalBatch(MappingInfoBatch &outBatch) const;

    void computePhysicalFirstDerivativesBatch(const ShapeInfoBatch &shape,
                                              MappingInfoBatch &out) const;
    };

}

#include "ReferenceElement.tpp"
