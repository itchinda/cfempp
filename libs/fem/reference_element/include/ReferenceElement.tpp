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
        const bool wantHessians = hasFlag(requestedFields, ShapeEvaluationFlags::Hessian);
        
        // 1. On prépare uniquement l'objet shape (qui est non-const ici)
        shape.setup(m_numNodes, wantHessians);

        if (shape.shapeMustBeUpdated) {
            // Reference gradients are strictly required to build the Jacobian
            // ShapeEvaluationFlags internalFields = requestedFields | ShapeEvaluationFlags::Gradient;
            evaluate(xi, requestedFields, shape);
        }

        // 2. Délégation à la version const de la fonction !
        // On cast 'shape' en référence const pour forcer le compilateur à utiliser la 2ème fonction
        computeMapping(xi, physicalNodesCoords, requestedFields, static_cast<const ShapeInfo&>(shape), out);
    }

    inline void ReferenceElement::computeMapping(const Vector3D &xi,
                                                 const DynamicVector<Vector3D> &physicalNodesCoords,
                                                 ShapeEvaluationFlags requestedFields,
                                                 const ShapeInfo &shape,
                                                 MappingInfo &out) const
    {
        const int n = m_numNodes;
        const bool wantPhysicalPosition = hasFlag(requestedFields, ShapeEvaluationFlags::PhysicalPosition);
        const bool wantJacobian = hasFlag(requestedFields, ShapeEvaluationFlags::Jacobian);
        const bool wantGradients = hasFlag(requestedFields, ShapeEvaluationFlags::Gradient);
        const bool wantHessians = hasFlag(requestedFields, ShapeEvaluationFlags::Hessian);

        out.setup(n, wantHessians);

        CFEM_ASSERT(shape.values.size() >= m_numNodes, "Shape Info is not initialized?");

        if (wantPhysicalPosition){
            out.physicalPosition = {0, 0, 0};
            for (int i = 0; i < n; ++i){
                out.physicalPosition += physicalNodesCoords[i] * shape.values[i];
            }
        }
        
        if (wantJacobian || wantGradients || wantHessians) {
            CFEM_ASSERT(n > 1, "Jacobian/Gradients requested for P0 element!");
            computeJacobian(shape, physicalNodesCoords, out);
        }

        if (wantGradients){
            computePhysicalGradients(shape, out);
        }

        if (wantHessians) {
            computePhysicalHessians(physicalNodesCoords, shape, out);
        }
    }

    inline void ReferenceElement::computeJacobian(const ShapeInfo &shape,
                                                  const DynamicVector<Vector3D> &physicalNodesCoords,
                                                  MappingInfo &out) const
    {
        const int n = static_cast<int>(shape.gradients.size()); 
        CFEM_ASSERT(n>0, "Shape Gradient Info must be ready prior to calling computeJacobian" );
        const int dim = m_dimension; // Idéalement, transformez dimension() en variable membre à l'avenir.

        // OPTIMISATION MAJEURE : Accumulation dans les registres CPU.
        // On initialise 9 doubles locaux. Le compilateur les placera dans les registres xmm/ymm.
        CfemReal j00 = 0.0, j01 = 0.0, j02 = 0.0;
        CfemReal j10 = 0.0, j11 = 0.0, j12 = 0.0;
        CfemReal j20 = 0.0, j21 = 0.0, j22 = 0.0;

        for (int i = 0; i < n; ++i)
        {
            const Vector3D &X = physicalNodesCoords[i];
            const Vector3D &gradN = shape.gradients[i];

            j00 += X.x * gradN.x; j10 += X.y * gradN.x; j20 += X.z * gradN.x;
            j01 += X.x * gradN.y; j11 += X.y * gradN.y; j21 += X.z * gradN.y;
            j02 += X.x * gradN.z; j12 += X.y * gradN.z; j22 += X.z * gradN.z;
        }

        // On écrit le résultat dans la mémoire structure une seule fois à la fin.
        out.jacobian(0, 0) = j00; out.jacobian(0, 1) = j01; out.jacobian(0, 2) = j02;
        out.jacobian(1, 0) = j10; out.jacobian(1, 1) = j11; out.jacobian(1, 2) = j12;
        out.jacobian(2, 0) = j20; out.jacobian(2, 1) = j21; out.jacobian(2, 2) = j22;

        if (dim == 3)
        {
            out.detJ = det(out.jacobian);
            out.invJacobian = inverse(out.jacobian, out.detJ);
        }
        else if (dim == 2)
        {
            // Vecteurs tangents
            CfemReal g11 = j00*j00 + j10*j10 + j20*j20;
            CfemReal g22 = j01*j01 + j11*j11 + j21*j21;
            CfemReal g12 = j00*j01 + j10*j11 + j20*j21;

            CfemReal detG = g11 * g22 - g12 * g12;
            out.detJ = std::sqrt(std::max(0.0, detG)); 

            CFEM_ASSERT(!NumericPolicy::areEqual(out.detJ, 0.0), "Degenerate 2D mapping detected");
            CfemReal invDetG = 1.0 / detG;
            for (int i = 0; i < 3; ++i)
            {
                out.invJacobian(0, i) = (g22 * out.jacobian(i, 0) - g12 * out.jacobian(i, 1)) * invDetG;
                out.invJacobian(1, i) = (g11 * out.jacobian(i, 1) - g12 * out.jacobian(i, 0)) * invDetG;
                out.invJacobian(2, i) = 0.0;
            }
        }
        else if (dim == 1)
        {
            CfemReal g11 = j00*j00 + j10*j10 + j20*j20;
            out.detJ = std::sqrt(g11); 
            CfemReal invG11 = 1.0 / g11;
            for (int i = 0; i < 3; ++i)
            {
                out.invJacobian(i, 0) = out.jacobian(i, 0) * invG11;
                out.invJacobian(i, 1) = 0.0;
                out.invJacobian(i, 2) = 0.0;
            }
        }

    }

    inline void ReferenceElement::computePhysicalGradients(const ShapeInfo &shape,
                                                           MappingInfo &out) const
    {
        const int n = static_cast<int>(shape.gradients.size());
        
        // Cache de la transposée de l'inverse du Jacobien dans des variables locales.
        // J^{-T} appliqué au vecteur local.
        const CfemReal i00 = out.invJacobian(0,0), i01 = out.invJacobian(0,1), i02 = out.invJacobian(0,2);
        const CfemReal i10 = out.invJacobian(1,0), i11 = out.invJacobian(1,1), i12 = out.invJacobian(1,2);
        const CfemReal i20 = out.invJacobian(2,0), i21 = out.invJacobian(2,1), i22 = out.invJacobian(2,2);

        for (int a = 0; a < n; ++a)
        {
            const Vector3D& gRef = shape.gradients[a];
            
            out.physicalShapeGradients[a].x = i00 * gRef.x + i10 * gRef.y + i20 * gRef.z;
            out.physicalShapeGradients[a].y = i01 * gRef.x + i11 * gRef.y + i21 * gRef.z;
            out.physicalShapeGradients[a].z = i02 * gRef.x + i12 * gRef.y + i22 * gRef.z;
        }
    }

    inline void ReferenceElement::computePhysicalHessians(const DynamicVector<Vector3D> &physicalNodesCoords,
                                                      const ShapeInfo &shape,
                                                      MappingInfo &out) const
    {
        const int n = static_cast<int>(shape.hessians.size());

        // 1. Calcul de la courbure géométrique (Hessienne de la transformation)
        // On utilise des variables locales pour accumuler les 6 composantes indépendantes (xx, yy, zz, xy, xz, yz)
        // pour chaque direction (X, Y, Z).
        CfemReal hX_xx = 0, hX_yy = 0, hX_zz = 0, hX_xy = 0, hX_xz = 0, hX_yz = 0;
        CfemReal hY_xx = 0, hY_yy = 0, hY_zz = 0, hY_xy = 0, hY_xz = 0, hY_yz = 0;
        CfemReal hZ_xx = 0, hZ_yy = 0, hZ_zz = 0, hZ_xy = 0, hZ_xz = 0, hZ_yz = 0;

        for (int a = 0; a < n; ++a)
        {
            const auto &hRef = shape.hessians[a];
            const auto &pos  = physicalNodesCoords[a];

            hX_xx += hRef.xx() * pos.x; hX_yy += hRef.yy() * pos.x; hX_zz += hRef.zz() * pos.x;
            hX_xy += hRef.xy() * pos.x; hX_xz += hRef.xz() * pos.x; hX_yz += hRef.yz() * pos.x;

            hY_xx += hRef.xx() * pos.y; hY_yy += hRef.yy() * pos.y; hY_zz += hRef.zz() * pos.y;
            hY_xy += hRef.xy() * pos.y; hY_xz += hRef.xz() * pos.y; hY_yz += hRef.yz() * pos.y;

            hZ_xx += hRef.xx() * pos.z; hZ_yy += hRef.yy() * pos.z; hZ_zz += hRef.zz() * pos.z;
            hZ_xy += hRef.xy() * pos.z; hZ_xz += hRef.xz() * pos.z; hZ_yz += hRef.yz() * pos.z;
        }

        // 2. Préparation des coefficients de l'inverse du Jacobien (Registres)
        const CfemReal i00 = out.invJacobian(0,0), i01 = out.invJacobian(0,1), i02 = out.invJacobian(0,2);
        const CfemReal i10 = out.invJacobian(1,0), i11 = out.invJacobian(1,1), i12 = out.invJacobian(1,2);
        const CfemReal i20 = out.invJacobian(2,0), i21 = out.invJacobian(2,1), i22 = out.invJacobian(2,2);

        // Helper lambda pour le calcul de chaque composante p,q
        auto computeComp = [&](int p, int q, CfemReal r_xx, CfemReal r_yy, CfemReal r_zz, 
                                             CfemReal r_xy, CfemReal r_xz, CfemReal r_yz) {
            // Puisque R est symétrique : R01=R10, R02=R20, R12=R21
            // Termes diagonaux de R : i=j
            CfemReal res = out.invJacobian(0,p)*r_xx*out.invJacobian(0,q) + 
                           out.invJacobian(1,p)*r_yy*out.invJacobian(1,q) + 
                           out.invJacobian(2,p)*r_zz*out.invJacobian(2,q);
            // Termes croisés de R : i != j (facteur 2 car R_ij = R_ji)
            res += (out.invJacobian(0,p)*r_xy*out.invJacobian(1,q) + out.invJacobian(1,p)*r_xy*out.invJacobian(0,q)) +
                   (out.invJacobian(0,p)*r_xz*out.invJacobian(2,q) + out.invJacobian(2,p)*r_xz*out.invJacobian(0,q)) +
                   (out.invJacobian(1,p)*r_yz*out.invJacobian(2,q) + out.invJacobian(2,p)*r_yz*out.invJacobian(1,q));
            return res;
        };

        // 3. Transformation vers l'espace physique
        for (int a = 0; a < n; ++a)
        {
            const auto &hRef = shape.hessians[a];
            const CfemReal dNdx = out.dNdx(a);
            const CfemReal dNdy = out.dNdy(a);
            const CfemReal dNdz = out.dNdz(a);

            // R = H_ref - sum( dN/dx_k * H_geom_k )
            // On calcule les 6 composantes de R directement
            CfemReal r_xx = hRef.xx() - (dNdx * hX_xx + dNdy * hY_xx + dNdz * hZ_xx);
            CfemReal r_yy = hRef.yy() - (dNdx * hX_yy + dNdy * hY_yy + dNdz * hZ_yy);
            CfemReal r_zz = hRef.zz() - (dNdx * hX_zz + dNdy * hY_zz + dNdz * hZ_zz);
            CfemReal r_xy = hRef.xy() - (dNdx * hX_xy + dNdy * hY_xy + dNdz * hZ_xy);
            CfemReal r_xz = hRef.xz() - (dNdx * hX_xz + dNdy * hY_xz + dNdz * hZ_xz);
            CfemReal r_yz = hRef.yz() - (dNdx * hX_yz + dNdy * hY_yz + dNdz * hZ_yz);
            
            // Transformation de congruence optimisée (équivalent à J^-T * R * J^-1) 
            // Pour une matrice symétrique R, la composante (p,q) de J^-T R J^-1 est:
            // H_phys_pq = sum_i sum_j (invJ_ip * R_ij * invJ_jq)
            auto& hPhys = out.physicalShapeHessians[a];
            
            // On ne calcule que la partie triangulaire supérieure (SymTensorO2)

            hPhys.xx() = computeComp(0, 0, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
            hPhys.yy() = computeComp(1, 1, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
            hPhys.zz() = computeComp(2, 2, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
            hPhys.xy() = computeComp(0, 1, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
            hPhys.xz() = computeComp(0, 2, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
            hPhys.yz() = computeComp(1, 2, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
        }
    }

    // ----------------------------------------------------------------------
    // ------------------------------- interpolation ------------------------
    template <typename T>
    inline T ReferenceElement::interpolate(const ShapeInfo &shape, const DynamicVector<T> &nodalValues) const
    {
        const int n = m_numNodes;
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
        const int n = m_numNodes;
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
        const int n = m_numNodes;
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