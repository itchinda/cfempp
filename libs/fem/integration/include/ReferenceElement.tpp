#ifndef CFEM_REFERENCE_ELEMENT_TPP_
#define CFEM_REFERENCE_ELEMENT_TPP_

#include "ReferenceElement.h"

#endif

namespace cfem
{

    inline void ReferenceElement::computeMapping(const Vector3D &xi,
                                                 const DynamicVector<Vector3D> &physicalNodesCoords,
                                                 ShapeEvaluationFlags requestedFields,
                                                 ShapeInfo &shape,
                                                 MappingInfo &out) const
    {
        const int n = numNodes();
        const bool wantHessians = hasFlag(requestedFields, ShapeEvaluationFlags::Hessian);

        shape.setup(n, wantHessians);
        out.setup(n, wantHessians);

        // Reference gradients are strictly required to build the Jacobian
        ShapeEvaluationFlags internalFields = requestedFields | ShapeEvaluationFlags::Gradient;
        evaluate(xi, internalFields, shape);

        out.physicalPosition = {0, 0, 0};
        for (int i = 0; i < n; ++i)
        {
            out.physicalPosition += physicalNodesCoords[i] * shape.values[i];
        }

        // Compute Metrics (Jacobian, Determinant, Inverse)
        computeJacobian(shape, physicalNodesCoords, out);

#ifndef NDEBUG
        int status = checkInversion(shape, physicalNodesCoords);
        if (status == -1)
        {
            CFEM_WARNING("Inverted element detected (detJ < 0)");
        }
        else if (status == 0)
        {
            CFEM_ASSERT(false, ("Degenerated element detected: The determinant of the jacobian is near zero. detJ = " + std::to_string(out.detJ) + ".").c_str());
        }
#endif

        if (hasFlag(requestedFields,  ShapeEvaluationFlags::Gradient))
        {
            computePhysicalGradients(shape, out);
        }

        if (wantHessians)
        {
            computePhysicalHessians(physicalNodesCoords, shape, out);
        }
    }

    // Compute the jacobian and inverse and the determinant
    inline void ReferenceElement::computeJacobian(const ShapeInfo &shape,
                                                  const DynamicVector<Vector3D> &physicalNodesCoords,
                                                  MappingInfo &out) const
    {
        const int dim = dimension();
        const int n = numNodes();

        // 1. Accumulate Jacobian components
        // J_jk = dx_j / dxi_k
        out.jacobian.setZero();
        for (int i = 0; i < n; ++i)
        {
            const Vector3D &X = physicalNodesCoords[i];
            const Vector3D &gradN = shape.gradients[i];

            // Column 0: dx/dxi
            out.jacobian(0, 0) += X.x * gradN.x;
            out.jacobian(1, 0) += X.y * gradN.x;
            out.jacobian(2, 0) += X.z * gradN.x;
            // Column 1: dx/deta (for 2D/3D)
            out.jacobian(0, 1) += X.x * gradN.y;
            out.jacobian(1, 1) += X.y * gradN.y;
            out.jacobian(2, 1) += X.z * gradN.y;
            // Column 2: dx/dzeta (for 3D)
            out.jacobian(0, 2) += X.x * gradN.z;
            out.jacobian(1, 2) += X.y * gradN.z;
            out.jacobian(2, 2) += X.z * gradN.z;
        }

        if (dim == 3)
        {
            out.detJ = det(out.jacobian);
            out.invJacobian = inverse(out.jacobian, out.detJ);
        }
        else if (dim == 2)
        {
            // Vecteurs tangents (colonnes de J)
            Vector3D g1{out.jacobian(0, 0), out.jacobian(1, 0), out.jacobian(2, 0)};
            Vector3D g2{out.jacobian(0, 1), out.jacobian(1, 1), out.jacobian(2, 1)};

            // Tenseur métrique G = J^T * J (2x2)
            CfemReal g11 = dot(g1, g1);
            CfemReal g12 = dot(g1, g2);
            CfemReal g22 = dot(g2, g2);

            CfemReal detG = g11 * g22 - g12 * g12;
            out.detJ = std::sqrt(std::max(0.0, detG)); // La mesure de surface

            // Calcul de la pseudo-inverse J+ = (J^T J)^-1 * J^T
            // On la stocke "à plat" (2 lignes, 3 colonnes physiques) 
            // pour qu'elle soit le J^-1 du monde 2D.
            CfemReal invDetG = 1.0 / detG;
            for (int i = 0; i < 3; ++i)
            {
                // Ligne 0 (xi) : combinaison de g1 et g2 pour la direction i
                out.invJacobian(0, i) = (g22 * out.jacobian(i, 0) - g12 * out.jacobian(i, 1)) * invDetG;
                // Ligne 1 (eta)
                out.invJacobian(1, i) = (g11 * out.jacobian(i, 1) - g12 * out.jacobian(i, 0)) * invDetG;
                // Ligne 2 (zeta) - n'existe pas en 2D
                out.invJacobian(2, i) = 0.0;
            }
        }
        else if (dim == 1)
        {
            // Tangent vector: g1 = dx/dxi
            Vector3D g1{out.jacobian(0, 0), out.jacobian(1, 0), out.jacobian(2, 0)};
            CfemReal g11 = dot(g1, g1);

            out.detJ = std::sqrt(g11); // This is the "length" metric

            CfemReal invG11 = 1.0 / g11;
            for (int i = 0; i < 3; ++i)
            {
                out.invJacobian(i, 0) = out.jacobian(i, 0) * invG11;
                out.invJacobian(i, 1) = 0.0;
                out.invJacobian(i, 2) = 0.0;
            }
        }
    }

    // compute the Gradients of the shape functions of the physical element
    inline void ReferenceElement::computePhysicalGradients(const ShapeInfo &shape,
                                                           MappingInfo &out) const
    {
        const int n = numNodes();
        out.physicalShapeGradients.resize(n);

        const TensorO2 invJT = out.invJacobian.transpose();
        for (int a = 0; a < n; ++a)
        {
            out.physicalShapeGradients[a] = invJT * shape.gradients[a];
        }
    }


    // compute the Hessians of the shape functions of the physical element
    inline void ReferenceElement::computePhysicalHessians(const DynamicVector<Vector3D> &physicalNodesCoords,
                                                          const ShapeInfo &shape,
                                                          MappingInfo &out) const
    {

        const int n = numNodes();
        CFEM_ASSERT(out.physicalShapeGradients.size() == static_cast<size_t>(n),
                    "Physical gradients must be computed before physical Hessians.");

        out.physicalShapeHessians.resize(n);

        // Geometric Curvature (Hessian of the mapping itself)
        SymTensorO2 H_X, H_Y, H_Z;
        H_X.setZero();
        H_Y.setZero();
        H_Z.setZero();

        for (int a = 0; a < n; ++a)
        {
            H_X += shape.hessians[a] * physicalNodesCoords[a].x;
            H_Y += shape.hessians[a] * physicalNodesCoords[a].y;
            H_Z += shape.hessians[a] * physicalNodesCoords[a].z;
        }

        // Second-order chain rule transformation
        for (int a = 0; a < n; ++a)
        {
            SymTensorO2 R = shape.hessians[a];
            R -= H_X * out.dNdx(a);
            R -= H_Y * out.dNdy(a);
            R -= H_Z * out.dNdz(a);

            // Congruent transformation: J^-T * R * J^-1
            out.physicalShapeHessians[a] = R.quadform(out.invJacobian);
        }
    }

    // Checking the connectivity and quality (degeneracy status)
    inline int ReferenceElement::checkInversion(const ShapeInfo &shape,
                                                const DynamicVector<Vector3D> &physicalNodesCoords,
                                                const Vector3D &refNormal) const
    {
        const int dim = dimension();
        const int n = numNodes();

        // 1. Compute the tangent vectors (Jacobian columns)
        Vector3D g1{0, 0, 0}, g2{0, 0, 0}, g3{0, 0, 0};
        for (int i = 0; i < n; ++i)
        {
            const Vector3D &gradN = shape.gradients[i];
            const Vector3D &X = physicalNodesCoords[i];

            g1 += X * gradN.x;
            if (dim >= 2)
                g2 += X * gradN.y;
            if (dim == 3)
                g3 += X * gradN.z;
        }

        if (dim == 3)
        {
            // Scalar Triple Product: (g1 x g2) . g3
            CfemReal vol = dot(cross(g1, g2), g3);
            if (vol < -ZERO_TOLERANCE)
                return -1; // Inverted
            if (vol < ZERO_TOLERANCE)
                return 0; // Degenerate
            return 1;     // Valid
        }
        else if (dim == 2)
        {
            Vector3D localNormal = cross(g1, g2);
            CfemReal areaSq = localNormal.normSquare();

            if (areaSq < ZERO_TOLERANCE * ZERO_TOLERANCE)
                return 0;

            // Check orientation against reference
            return (dot(localNormal, refNormal) > 0) ? 1 : -1;
        }

        return 1; // 1D segments are never "inverted" in this geometric sense
    }

    // ----------------------------------------------------------------------
    // ------------------------------- interpolation ------------------------
    template <typename T>
    inline T ReferenceElement::interpolate(const ShapeInfo &shape, const DynamicVector<T> &nodalValues) const
    {
        const int n = numNodes();
        CFEM_ASSERT(shape.values.size() >= n, "ShapeInfo not evaluated for this element type");
        CFEM_ASSERT(nodalValues.size() == n, "Nodal values size mismatch");

        T result{}; // T must have a default constructor that initialize the field with 0.

        for (int i = 0; i < n; ++i)
        {
            result += nodalValues[i] * shape.values[i];
        }
        return result;
    }

    inline Vector3D ReferenceElement::interpolateGradient(const MappingInfo &map,
                                                          const DynamicVector<CfemReal> &nodalValues) const
    {
        const int n = numNodes();
        CFEM_ASSERT(map.physicalShapeGradients.size() >= n, "MappingInfo gradients not computed");
        CFEM_ASSERT(nodalValues.size() == n, "Nodal values size mismatch");

        Vector3D grad{0, 0, 0};
        for (int i = 0; i < n; ++i)
        {
            grad += map.physicalShapeGradients[i] * nodalValues[i];
        }
        return grad;
    }

    inline SymTensorO2 ReferenceElement::interpolateHessian(const MappingInfo &map,
                                                            const DynamicVector<CfemReal> &nodalValues) const
    {
        const int n = numNodes();
        CFEM_ASSERT(map.physicalShapeHessians.size() >= n, "MappingInfo hessians not computed");
        CFEM_ASSERT(nodalValues.size() == n, "Nodal values size mismatch");

        SymTensorO2 hess;
        hess.setZero();
        for (int i = 0; i < n; ++i)
        {
            hess += map.physicalShapeHessians[i] * nodalValues[i];
        }
        return hess;
    }

}  // namespace cfem