#ifndef CFEM_REFERENCE_ELEMENT_
#define CFEM_REFERENCE_ELEMENT_

#include "cfemutils.h"
#include "Vector3D.h"
#include "DynamicVector.h"
#include "DynamicMatrix.h"
#include "TensorO2.h"
#include "ShapeEvaluationFlags.h"

namespace cfem {

/**
 * @struct ShapeInfo
 * @brief Storage buffer for shape function evaluations at a specific quadrature point.
 * @note Re-use this object across multiple integration points to avoid frequent heap allocations.
 */
struct ShapeInfo {
    DynamicVector<CfemReal> values;      ///< Shape function values (N_i)
    DynamicVector<Vector3D> gradients;   ///< Local gradients w.r.t reference coordinates (dN_i/dxi)
    DynamicVector<SymTensorO2> hessians; ///< Local Hessians w.r.t reference coordinates (d2N_i/dxi2)
    bool shapeMustBeUpdated = true;
    ShapeEvaluationFlags  validFlags = ShapeEvaluationFlags::None;

    /** 
     * @brief Prepares buffer sizes to prevent runtime reallocations. 
     * @param numNodes Number of nodes in the element.
     * @param includeHessians If true, allocates memory for Hessian tensors.
     */
    void setup(int numNodes, bool includeHessians = false) {
        if (values.size() != static_cast<size_t>(numNodes)) {
            values.resize(numNodes);
            gradients.resize(numNodes);
        }
        if (includeHessians && hessians.size() != static_cast<size_t>(numNodes)) {
            hessians.resize(numNodes);
        }
        validFlags = ShapeEvaluationFlags::None;
    }
    
    /** @brief Checks if a specific field (Value/Gradient/Hessian) is marked as valid. */
    bool has(cfem::ShapeEvaluationFlags field) const {
        return hasFlag(validFlags, field);
    }
};

/**
 * @struct MappingInfo
 * @brief Stores metric data resulting from the geometric transformation (Reference -> Physical).
 */
struct MappingInfo {
    Vector3D physicalPosition; ///< The mapped x, y, z coordinates
    TensorO2 jacobian;         ///< Jacobian matrix (dx/dxi)
    TensorO2 invJacobian;      ///< Inverse Jacobian matrix (dxi/dx)
    CfemReal detJ;             ///< Determinant of the Jacobian

    DynamicVector<Vector3D> physicalShapeGradients;   ///< dN_i/dx (Gradients in physical space)
    DynamicVector<SymTensorO2> physicalShapeHessians; ///< d2N_i/dx2 (Hessians in physical space)

    /**
     * @brief Prepares the mapping buffers for a specific element size.
     * @param numNodes Number of nodes in the element.
     * @param includeHessians Whether to allocate the Hessian buffer.
     */
    void setup(int numNodes, bool includeHessians = false) {
        if (physicalShapeGradients.size() != static_cast<size_t>(numNodes)) {
            physicalShapeGradients.resize(numNodes);
        }
        
        if (includeHessians && physicalShapeHessians.size() != static_cast<size_t>(numNodes)) {
            physicalShapeHessians.resize(numNodes);
        }
    }

    // Quick accessors for physical gradients
    inline CfemReal dNdx(int i) const { return physicalShapeGradients[i].x; }
    inline CfemReal dNdy(int i) const { return physicalShapeGradients[i].y; }
    inline CfemReal dNdz(int i) const { return physicalShapeGradients[i].z; }

    // Quick accessors for physical Hessians (components)
    inline CfemReal d2Ndxdx(int nodeIdx) const { return physicalShapeHessians[nodeIdx].xx(); }
    inline CfemReal d2Ndydy(int nodeIdx) const { return physicalShapeHessians[nodeIdx].yy(); }
    inline CfemReal d2Ndzdz(int nodeIdx) const { return physicalShapeHessians[nodeIdx].zz(); }
    inline CfemReal d2Ndydz(int nodeIdx) const { return physicalShapeHessians[nodeIdx].yz(); }
    inline CfemReal d2Ndxdz(int nodeIdx) const { return physicalShapeHessians[nodeIdx].xz(); }
    inline CfemReal d2Ndxdy(int nodeIdx) const { return physicalShapeHessians[nodeIdx].xy(); }

    void computePhysicalGradientsOnly(const ShapeInfo& fieldShape) {
        // Le fait de charger dans des variables locales (registres) 
        // empêche le CPU de faire des allers-retours mémoire.
        const CfemReal i00 = invJacobian(0,0), i01 = invJacobian(0,1), i02 = invJacobian(0,2);
        const CfemReal i10 = invJacobian(1,0), i11 = invJacobian(1,1), i12 = invJacobian(1,2);
        const CfemReal i20 = invJacobian(2,0), i21 = invJacobian(2,1), i22 = invJacobian(2,2);

        for (size_t i = 0; i < fieldShape.gradients.size(); ++i) {
            const auto& g = fieldShape.gradients[i];
            physicalShapeGradients[i].x = i00 * g.x + i10 * g.y + i20 * g.z;
            physicalShapeGradients[i].y = i01 * g.x + i11 * g.y + i21 * g.z;
            physicalShapeGradients[i].z = i02 * g.x + i12 * g.y + i22 * g.z;
        }
    }
};

/**
 * @class ReferenceElement
 * @brief Abstract base class for Finite Element reference definitions.
 * Provides the machinery to compute geometric mappings, Jacobians, and 
 * physical derivatives from local shape functions.
 */
class ReferenceElement {
protected:
    DynamicVector<Vector3D> m_nodes; ///< Node coordinates in reference space (xi)

public:
    virtual ~ReferenceElement() = default;

    virtual int dimension() const = 0;
    virtual int getNumNodes() const = 0;

    /**
     * @brief Evaluates shape functions and their local derivatives at point xi.
     * @param xi Target point in reference space.
     * @param requested Bitmask of fields to calculate (Value, Gradient, Hessian).
     * @param out ShapeInfo buffer to be filled.
     */
    virtual void evaluate(const Vector3D& xi, 
                          ShapeEvaluationFlags requested, 
                          ShapeInfo& out) const = 0;
                          
    virtual Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const = 0;
    /**
     * @brief Returns the coordinates of the nodes in the reference space. 
     */
    const DynamicVector<Vector3D>& referenceNodes() const { return m_nodes; }
    const Vector3D& getNodeCoords( CfemInt i) const { CFEM_ASSERT( i < m_nodes.size() && i >= 0 ); return m_nodes[i]; }

    /**
     * @brief Computes the full geometric mapping at a given point xi.
     * @param xi Target point in reference space.
     * @param physicalNodesCoords Actual physical coordinates of the element's nodes.
     * @param requestedFields Bitmask of physical fields required (Gradients, Hessians).
     * @param shape Work buffer for reference shape functions.
     * @param out MappingInfo structure to be populated.
     */
    void computeMapping(const Vector3D& xi,
                        const DynamicVector<Vector3D>& physicalNodesCoords, 
                        ShapeEvaluationFlags requestedFields,
                        ShapeInfo& shape, 
                        MappingInfo& out) const;

    inline void computeMapping(const Vector3D &xi,
                                                 const DynamicVector<Vector3D> &physicalNodesCoords,
                                                 ShapeEvaluationFlags requestedFields,
                                                 const ShapeInfo &shape,
                                                 MappingInfo &out) const;

    /**
     * @brief Computes the Jacobian, its determinant (metric), and the pseudo-inverse.
     * @note This version supports 1D/2D elements embedded in 3D space using 
     * the Moore-Penrose pseudo-inverse logic.
     */
    void computeJacobian(const ShapeInfo& shape, 
                        const DynamicVector<Vector3D>& physicalNodesCoords, 
                        MappingInfo& out) const;


    /** 
     * @brief Transforms reference gradients to physical space using J^-T. 
     * @param shape Work buffer for reference shape functions.
     * @param out MappingInfo structure to be populated.
     */
    void computePhysicalGradients(const ShapeInfo& shape, 
                                  MappingInfo& out) const;

    /**
     * @brief Computes physical Hessians including the geometric curvature term.
     * @param physicalNodesCoords Actual physical coordinates of the element's nodes.
     * @param shape Work buffer for reference shape functions.
     * @param out MappingInfo structure to be populated.
     * @note Reference: H_phys = J^-T * [H_ref - grad_phys(N) * H_geom] * J^-1
     */
    void computePhysicalHessians(const DynamicVector<Vector3D>& physicalNodesCoords, 
                                const ShapeInfo& shape, 
                                MappingInfo& out) const;

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
     * @return SymTensorO2 representing the 2nd derivative matrix.
     */
    SymTensorO2 interpolateHessian(const MappingInfo& map, 
                                   const DynamicVector<CfemReal>& nodalValues) const;
};

}

#include "ReferenceElement.tpp"

#endif