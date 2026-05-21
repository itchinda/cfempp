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

namespace cfem
{

    inline void ReferenceElement::evaluateShapes(const Vector3D &xi, ShapeEvaluationFlags requestedShape, ShapeInfo& shape) const
    {
        const bool includeSecondDerivatives = hasFlag(requestedShape, ShapeEvaluationFlags::SecondDerivatives);
        shape.setup(m_numNodes, includeSecondDerivatives);

        if (hasFlag(requestedShape, ShapeEvaluationFlags::Value)) {
            evaluateShapesValues(xi, shape.values.data());
        }
        if (hasFlag(requestedShape, ShapeEvaluationFlags::FirstDerivatives)) {
            evaluateShapesDerivatives(xi, shape.dN_dxi.data(), shape.dN_deta.data(), shape.dN_dzeta.data());
        }
        if (includeSecondDerivatives) {
            evaluateShapesSecondDerivatives(xi, 
                                      shape.d2N_dxi2.data(), shape.d2N_deta2.data(), shape.d2N_dzeta2.data(),
                                      shape.d2N_dxideta.data(), shape.d2N_dxidzeta.data(), shape.d2N_detadzeta.data());
        }
        shape.validFlags = requestedShape;

    }

    inline void ReferenceElement::evaluateShapesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeEvaluationFlags requestedShape, ShapeInfoBatch& shapeBatch) const
    {
        const bool includeSecondDerivatives = hasFlag(requestedShape, ShapeEvaluationFlags::SecondDerivatives);
        if (hasFlag(requestedShape, ShapeEvaluationFlags::Value)) {
            evaluateShapesValuesBatch(xiPoints, shapeBatch.values);
        }
        if (hasFlag(requestedShape, ShapeEvaluationFlags::FirstDerivatives)) {
            evaluateShapesDerivativesBatch(xiPoints, 
                                shapeBatch.dN_dxi, 
                                shapeBatch.dN_deta, 
                                shapeBatch.dN_dzeta);
        }
        if (includeSecondDerivatives) {
            evaluateShapesSecondDerivativesBatch(xiPoints,
                                      shapeBatch.d2N_dxi2, shapeBatch.d2N_deta2, shapeBatch.d2N_dzeta2,
                                      shapeBatch.d2N_dxideta, shapeBatch.d2N_dxidzeta, shapeBatch.d2N_detadzeta);
        }

        shapeBatch.validFlags = requestedShape;
        
    }
    // =========================================================================
    // ÉVALUATIONS PONCTUELLES (SINGLE POINT)
    // =========================================================================

    inline void ReferenceElement::computeMapping(const Vector3D &xi,
                                                 const DynamicVector<Vector3D> &physicalNodesCoords,
                                                 ShapeEvaluationFlags requestedShape,
                                                 MappingEvaluationFlags requestedMapping,
                                                 ShapeInfo &shape,
                                                 MappingInfo &out) const
    {
        // Évaluation polymorphe des fonctions de forme sur l'élément de référence
        if (shape.shapeMustBeUpdated){
            evaluateShapes(xi, requestedShape, shape);
        }
        
        // Passage direct à la variante géométrique passive via const ShapeInfo&
        computeMapping(xi, physicalNodesCoords, requestedMapping, shape, out);
    }

    inline void ReferenceElement::computeMapping(const Vector3D &xi,
                                                 const DynamicVector<Vector3D> &physicalNodesCoords,
                                                 MappingEvaluationFlags requestedMapping,
                                                 const ShapeInfo &shape,
                                                 MappingInfo &out) const
    {
        // Préparation des buffers physiques (gère l'allocation ou la réutilisation)
        out.setup(m_numNodes, hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalSecondDerivatives));

        // Position physique interpolée
        if (hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalPosition)) {
            out.physicalPosition.setZero();
            const CfemReal* N = shape.values.data();
            for (CfemInt i = 0; i < m_numNodes; ++i) {
                out.physicalPosition += physicalNodesCoords[i] * N[i];
            }
        }

        // Base locale / Jacobien
        if (hasFlag(requestedMapping, MappingEvaluationFlags::Jacobian)) {
            computeJacobian(shape, physicalNodesCoords, out);
        }

        // Déterminant (Calcule et stocke le tenseur métrique SymTensor2D/3D pour éviter le double calcul)
        if (hasFlag(requestedMapping, MappingEvaluationFlags::Determinant)) {
            computeDeterminant(shape, physicalNodesCoords, out);
        }

        // Inverse / Pseudo-inverse du Jacobien (Réutilise le tenseur métrique mis en cache)
        if (hasFlag(requestedMapping, MappingEvaluationFlags::InverseJacobian)) {
            computeInverseJacobian(shape, physicalNodesCoords, out);
        }

        // Gradients physiques SoA (dN/dx, dN/dy, dN/dz)
        if (hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalFirstDerivatives)) {
            computePhysicalFirstDerivatives(shape, out);
        }

        // Hessiens physiques translatés (avec prise en compte de la courbure géométrique)
        if (hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalSecondDerivatives)) {
            computePhysicalSecondDerivatives(physicalNodesCoords, shape, out);
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::Normal)){
            computeGeometricNormal(shape, out);
        }

        out.validFlags = requestedMapping;
    }

    // =========================================================================
    // ÉVALUATIONS PAR LOT (BATCH)
    // =========================================================================

    inline void ReferenceElement::computeMappingBatch(const DynamicVector<Vector3D> &xiPoints,
                                                      const DynamicVector<Vector3D> &physicalNodesCoords,
                                                      ShapeEvaluationFlags requestedShape,
                                                      MappingEvaluationFlags requestedMapping,
                                                      ShapeInfoBatch &shapeBatch,
                                                      MappingInfoBatch &outBatch) const
    {        
        // Évaluation par lignes contiguës dans les DynamicMatrix
        if (shapeBatch.shapeMustBeUpdated){
            evaluateShapesBatch(xiPoints, requestedShape, shapeBatch);
        }

        // Redirection stricte vers la variante const pour les transformations géométriques par lot
        computeMappingBatch(xiPoints, physicalNodesCoords, requestedMapping, shapeBatch, outBatch);
    }

    inline void ReferenceElement::computeMappingBatch(const DynamicVector<Vector3D> &xiPoints,
                                                      const DynamicVector<Vector3D> &physicalNodesCoords,
                                                      MappingEvaluationFlags requestedMapping,
                                                      const ShapeInfoBatch &shapeBatch,
                                                      MappingInfoBatch &outBatch) const
    {
        const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
        
        // Initialisation du lot physique
        outBatch.setup(nQp, m_numNodes, m_dimension, hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalSecondDerivatives));

        // Interpolation des positions physiques par point de quadrature
        if (hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalPosition)) {
            for (CfemInt q = 0; q < nQp; ++q) {
                outBatch.physicalPositions[q].setZero();
                const auto rowValues = shapeBatch.values[q];
                for (CfemInt i = 0; i < m_numNodes; ++i) {
                    outBatch.physicalPositions[q] += physicalNodesCoords[i] * rowValues[i];
                }
            }
        }

        // Évaluation des opérateurs physiques par lot (Chaque fonction traite l'ensemble du vecteur de points qp)
        if (hasFlag(requestedMapping, MappingEvaluationFlags::Jacobian)) {
            computeJacobianBatch(shapeBatch, physicalNodesCoords, outBatch);
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::Determinant)) {
            computeDeterminantBatch(shapeBatch, physicalNodesCoords, outBatch);
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::InverseJacobian)) {
            computeInverseJacobianBatch(shapeBatch, physicalNodesCoords, outBatch);
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalFirstDerivatives)) {
            computePhysicalFirstDerivativesBatch(shapeBatch, outBatch);
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalSecondDerivatives)) {
            computePhysicalSecondDerivativesBatch(physicalNodesCoords, shapeBatch, outBatch);
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::Normal)){
            computeGeometricNormalBatch(shapeBatch, outBatch);
        }

        outBatch.validFlags = requestedMapping;
    }

    inline void ReferenceElement::computeJacobian(const ShapeInfo& shape, 
                                                  const DynamicVector<Vector3D>& physicalNodesCoords, 
                                                  MappingInfo& out) const
    {
        const CfemInt n = shape.numNodes;

        CfemReal j00 = 0.0, j01 = 0.0, j02 = 0.0;
        CfemReal j10 = 0.0, j11 = 0.0, j12 = 0.0;
        CfemReal j20 = 0.0, j21 = 0.0, j22 = 0.0;

        // Pointeurs SoA pour ce point de quadrature précis
        const CfemReal* dxi   = shape.dN_dxi.data();
        const CfemReal* deta  = shape.dN_deta.data();
        const CfemReal* dzeta = shape.dN_dzeta.data();
        
        for (CfemInt i = 0; i < n; ++i)
        {
            const Vector3D &X = physicalNodesCoords[i];
            j00 += X.x * dxi[i];   j10 += X.y * dxi[i];   j20 += X.z * dxi[i];
            j01 += X.x * deta[i];  j11 += X.y * deta[i];  j21 += X.z * deta[i];
            j02 += X.x * dzeta[i]; j12 += X.y * dzeta[i]; j22 += X.z * dzeta[i];
        }

        out.jacobian.data[Tensor3D::idx_xx] = j00; out.jacobian.data[Tensor3D::idx_xy] = j01; out.jacobian.data[Tensor3D::idx_xz] = j02;
        out.jacobian.data[Tensor3D::idx_yx] = j10; out.jacobian.data[Tensor3D::idx_yy] = j11; out.jacobian.data[Tensor3D::idx_yz] = j12;
        out.jacobian.data[Tensor3D::idx_zx] = j20; out.jacobian.data[Tensor3D::idx_zy] = j21; out.jacobian.data[Tensor3D::idx_zz] = j22;
    }


    inline void ReferenceElement::computeJacobianBatch(const ShapeInfoBatch& shape, 
                                                       const DynamicVector<Vector3D>& physicalNodesCoords, 
                                                       MappingInfoBatch& out) const
    {
        const CfemInt nQp = shape.numQuadraturePoints;
        const CfemInt n = shape.numNodes;

        for (CfemInt q = 0; q < nQp; ++q)
        {
            CfemReal j00 = 0.0, j01 = 0.0, j02 = 0.0;
            CfemReal j10 = 0.0, j11 = 0.0, j12 = 0.0;
            CfemReal j20 = 0.0, j21 = 0.0, j22 = 0.0;

            // Pointeurs SoA pour ce point de quadrature précis
            const CfemReal* dxi   = shape.dN_dxi[q].data();
            const CfemReal* deta  = shape.dN_deta[q].data();
            const CfemReal* dzeta = shape.dN_dzeta[q].data();

            for (CfemInt i = 0; i < n; ++i)
            {
                const Vector3D &X = physicalNodesCoords[i];
                j00 += X.x * dxi[i];   j10 += X.y * dxi[i];   j20 += X.z * dxi[i];
                j01 += X.x * deta[i];  j11 += X.y * deta[i];  j21 += X.z * deta[i];
                j02 += X.x * dzeta[i]; j12 += X.y * dzeta[i]; j22 += X.z * dzeta[i];
            }

            out.jacobians[q].data[Tensor3D::idx_xx] = j00; out.jacobians[q].data[Tensor3D::idx_xy] = j01; out.jacobians[q].data[Tensor3D::idx_xz] = j02;
            out.jacobians[q].data[Tensor3D::idx_yx] = j10; out.jacobians[q].data[Tensor3D::idx_yy] = j11; out.jacobians[q].data[Tensor3D::idx_yz] = j12;
            out.jacobians[q].data[Tensor3D::idx_zx] = j20; out.jacobians[q].data[Tensor3D::idx_zy] = j21; out.jacobians[q].data[Tensor3D::idx_zz] = j22;
        }
    }

    inline void ReferenceElement::computeDeterminant(const ShapeInfo &shape,
                                                     const DynamicVector<Vector3D> &physicalNodesCoords,
                                                     MappingInfo &out) const
    {
        if (m_dimension == 3) {
            out.detJ =det(out.jacobian);
        } 
        else if (m_dimension == 2) {
            // Métrique de surface
            const CfemReal j00 = out.jacobian.data[Tensor3D::idx_xx], j01 = out.jacobian.data[Tensor3D::idx_xy];
            const CfemReal j10 = out.jacobian.data[Tensor3D::idx_yx], j11 = out.jacobian.data[Tensor3D::idx_yy];
            const CfemReal j20 = out.jacobian.data[Tensor3D::idx_zx], j21 = out.jacobian.data[Tensor3D::idx_zy];

            CfemReal g11 = j00*j00 + j10*j10 + j20*j20;
            CfemReal g22 = j01*j01 + j11*j11 + j21*j21;
            CfemReal g12 = j00*j01 + j10*j11 + j20*j21;

            out.metric2D = SymTensor2D(g11, g22, g12);
            out.detJ = std::sqrt(std::max(0.0, g11 * g22 - g12 * g12));
        }
        else if (m_dimension == 1) {
            const CfemReal j00 = out.jacobian.data[Tensor3D::idx_xx];
            const CfemReal j10 = out.jacobian.data[Tensor3D::idx_yx];
            const CfemReal j20 = out.jacobian.data[Tensor3D::idx_zx];
            out.metric1D = j00*j00 + j10*j10 + j20*j20;
            out.detJ = std::sqrt(out.metric1D);
        }
    }

    inline void ReferenceElement::computeDeterminantBatch(const ShapeInfoBatch &shape,
                                                     const DynamicVector<Vector3D> &physicalNodesCoords,
                                                     MappingInfoBatch &out) const
    {
        const CfemInt nQp = shape.numQuadraturePoints;

        if (m_dimension == 3) {
            for (CfemInt q = 0; q < nQp; ++q){
                out.detJs[q] = det(out.jacobians[q]);
            }
        } 
        else if (m_dimension == 2) {
            // Métrique de surface
            for (CfemInt q = 0; q < nQp; ++q){
                const CfemReal j00 = out.jacobians[q].data[Tensor3D::idx_xx], j01 = out.jacobians[q].data[Tensor3D::idx_xy];
                const CfemReal j10 = out.jacobians[q].data[Tensor3D::idx_yx], j11 = out.jacobians[q].data[Tensor3D::idx_yy];
                const CfemReal j20 = out.jacobians[q].data[Tensor3D::idx_zx], j21 = out.jacobians[q].data[Tensor3D::idx_zy];

                CfemReal g11 = j00*j00 + j10*j10 + j20*j20;
                CfemReal g22 = j01*j01 + j11*j11 + j21*j21;
                CfemReal g12 = j00*j01 + j10*j11 + j20*j21;

                out.metric2Ds[q] = SymTensor2D{g11, g22, g12};
                out.detJs[q] = std::sqrt(std::max(0.0, g11 * g22 - g12 * g12));
            }
        }
        else if (m_dimension == 1) {
            for (CfemInt q = 0; q < nQp; ++q){
                const CfemReal j00 = out.jacobians[q].data[Tensor3D::idx_xx];
                const CfemReal j10 = out.jacobians[q].data[Tensor3D::idx_yx];
                const CfemReal j20 = out.jacobians[q].data[Tensor3D::idx_zx];
                CfemReal g11 = j00*j00 + j10*j10 + j20*j20;
                
                out.metric1Ds[q] = g11;
                out.detJs[q] = std::sqrt(g11);
            }
        }

    }

    inline void ReferenceElement::computeInverseJacobian(const ShapeInfo &shape,
                                                         const DynamicVector<Vector3D> &physicalNodesCoords,
                                                         MappingInfo &out) const
    {
        CFEM_ASSERT(std::abs(out.detJ) > cfem::ZERO_TOLERANCE, "Degenerate Jacobian detected");

        if (m_dimension == 3) {
            out.invJacobian = inverse(out.jacobian, out.detJ);
        } 
        else if (m_dimension == 2) {
            const CfemReal g11 = out.metric2D.xx();
            const CfemReal g22 = out.metric2D.yy();
            const CfemReal g12 = out.metric2D.xy();

            CfemReal invDetG = 1.0 / (g11 * g22 - g12 * g12);

            for (CfemInt i = 0; i < 3; ++i) {
                out.invJacobian(0, i) = (g22 * out.jacobian(i, 0) - g12 * out.jacobian(i, 1)) * invDetG;
                out.invJacobian(1, i) = (g11 * out.jacobian(i, 1) - g12 * out.jacobian(i, 0)) * invDetG;
                out.invJacobian(2, i) = 0.0;
            }
        }
        else if (m_dimension == 1) {
            CfemReal invG11 = 1.0 / out.metric1D;
            for (CfemInt i = 0; i < 3; ++i) {
                out.invJacobian(i, 0) = out.jacobian(i, 0) * invG11;
                out.invJacobian(i, 1) = 0.0;
                out.invJacobian(i, 2) = 0.0;
            }
        }
    }

    inline void ReferenceElement::computeInverseJacobianBatch(const ShapeInfoBatch &shape,
                                                         const DynamicVector<Vector3D> &physicalNodesCoords,
                                                         MappingInfoBatch &out) const
    {
        const CfemInt nQp = shape.numQuadraturePoints;
        CFEM_ASSERT(std::abs(out.detJs[0]) > cfem::ZERO_TOLERANCE, "Degenerate Jacobian detected");


        if (m_dimension == 3) {
            for (CfemInt q = 0; q < nQp; ++q){
                out.invJacobians[q] = inverse(out.jacobians[q], out.detJs[q]);
            }
        } 
        else if (m_dimension == 2) {
            for (CfemInt q = 0; q < nQp; ++q){
                const CfemReal g11 = out.metric2Ds[q].xx();
                const CfemReal g22 = out.metric2Ds[q].yy();
                const CfemReal g12 = out.metric2Ds[q].xy();

                CfemReal invDetG = 1.0 / (g11 * g22 - g12 * g12);

                for (CfemInt i = 0; i < 3; ++i) {
                    out.invJacobians[q](0, i) = (g22 * out.jacobians[q](i, 0) - g12 * out.jacobians[q](i, 1)) * invDetG;
                    out.invJacobians[q](1, i) = (g11 * out.jacobians[q](i, 1) - g12 * out.jacobians[q](i, 0)) * invDetG;
                    out.invJacobians[q](2, i) = 0.0;
                }
            }
        }
        else if (m_dimension == 1) {
            for (CfemInt q = 0; q < nQp; ++q){
                CfemReal invG11 = 1.0 / out.metric1Ds[q];
                for (CfemInt i = 0; i < 3; ++i) {
                    out.invJacobians[q](i, 0) = out.jacobians[q](i, 0) * invG11;
                    out.invJacobians[q](i, 1) = 0.0;
                    out.invJacobians[q](i, 2) = 0.0;
                }
            }
        }
        
    }

    inline void ReferenceElement::computePhysicalFirstDerivatives(const ShapeInfo &shape,
                                                                  MappingInfo &out) const
    {
        const CfemInt n = shape.numNodes;
        
        // J^{-T} appliqué au vecteur local.
        const CfemReal i00 = out.invJacobian.xx(), i01 = out.invJacobian.xy(), i02 = out.invJacobian.xz();
        const CfemReal i10 = out.invJacobian.yx(), i11 = out.invJacobian.yy(), i12 = out.invJacobian.yz();
        const CfemReal i20 = out.invJacobian.zx(), i21 = out.invJacobian.zy(), i22 = out.invJacobian.zz();

        // Pointeurs SoA (lecture)
        const CfemReal* dxi   = shape.dN_dxi.data();
        const CfemReal* deta  = shape.dN_deta.data();
        const CfemReal* dzeta = shape.dN_dzeta.data();

        // Pointeurs SoA (écriture)
        CfemReal* dx = out.dN_dx.data();
        CfemReal* dy = out.dN_dy.data();
        CfemReal* dz = out.dN_dz.data();

        // Boucle SIMD ultra-rapide
        for (CfemInt a = 0; a < n; ++a)
        {
            dx[a] = i00 * dxi[a] + i10 * deta[a] + i20 * dzeta[a];
            dy[a] = i01 * dxi[a] + i11 * deta[a] + i21 * dzeta[a];
            dz[a] = i02 * dxi[a] + i12 * deta[a] + i22 * dzeta[a];
        }
    }

    inline void ReferenceElement::computePhysicalFirstDerivativesBatch(const ShapeInfoBatch &shape,
                                                                  MappingInfoBatch &out) const
    {
        const CfemInt nQp = shape.numQuadraturePoints;
        const CfemInt n = shape.numNodes;

        for (CfemInt q = 0; q < nQp; ++q)
        {
            // Cache des composantes de l'inverse du Jacobien dans les registres pour le point q
            const CfemReal i00 = out.invJacobians[q](0,0), i01 = out.invJacobians[q](0,1), i02 = out.invJacobians[q](0,2);
            const CfemReal i10 = out.invJacobians[q](1,0), i11 = out.invJacobians[q](1,1), i12 = out.invJacobians[q](1,2);
            const CfemReal i20 = out.invJacobians[q](2,0), i21 = out.invJacobians[q](2,1), i22 = out.invJacobians[q](2,2);

            // Pointeurs SoA vers les lignes de la matrice de référence (Lecture)
            const CfemReal* dxi   = shape.dN_dxi[q].data();
            const CfemReal* deta  = shape.dN_deta[q].data();
            const CfemReal* dzeta = shape.dN_dzeta[q].data();

            // Pointeurs SoA vers les lignes de la matrice physique (Écriture)
            CfemReal* dx = out.dN_dx[q].data();
            CfemReal* dy = out.dN_dy[q].data();
            CfemReal* dz = out.dN_dz[q].data();

            // Boucle sur les nœuds interne : vectorisation SIMD maximale
            for (CfemInt a = 0; a < n; ++a)
            {
                dx[a] = i00 * dxi[a] + i10 * deta[a] + i20 * dzeta[a];
                dy[a] = i01 * dxi[a] + i11 * deta[a] + i21 * dzeta[a];
                dz[a] = i02 * dxi[a] + i12 * deta[a] + i22 * dzeta[a];
            }
        }
    }

    inline void ReferenceElement::computePhysicalSecondDerivatives(const DynamicVector<Vector3D> &physicalNodesCoords,
                                                                   const ShapeInfo &shape,
                                                                   MappingInfo &out) const
    {
        const CfemInt n = shape.numNodes;

        // Étape 1 : Calcul de la courbure géométrique (H_geom) 
        CfemReal hX_xx = 0, hX_yy = 0, hX_zz = 0, hX_xy = 0, hX_xz = 0, hX_yz = 0;
        CfemReal hY_xx = 0, hY_yy = 0, hY_zz = 0, hY_xy = 0, hY_xz = 0, hY_yz = 0;
        CfemReal hZ_xx = 0, hZ_yy = 0, hZ_zz = 0, hZ_xy = 0, hZ_xz = 0, hZ_yz = 0;

        // Pointeurs SoA pour les Hessiens de référence
        const CfemReal* d2xi2      = shape.d2N_dxi2.data();
        const CfemReal* d2eta2     = shape.d2N_deta2.data();
        const CfemReal* d2zeta2    = shape.d2N_dzeta2.data();
        const CfemReal* d2xideta   = shape.d2N_dxideta.data();
        const CfemReal* d2xidzeta  = shape.d2N_dxidzeta.data();
        const CfemReal* d2etadzeta = shape.d2N_detadzeta.data();

        for (CfemInt a = 0; a < n; ++a)
        {
            const Vector3D &pos = physicalNodesCoords[a];
            
            // Chargement local depuis les tableaux
            CfemReal hxx = d2xi2[a];
            CfemReal hyy = d2eta2[a];
            CfemReal hzz = d2zeta2[a];
            CfemReal hxy = d2xideta[a];
            CfemReal hxz = d2xidzeta[a];
            CfemReal hyz = d2etadzeta[a];

            hX_xx += hxx * pos.x; hX_yy += hyy * pos.x; hX_zz += hzz * pos.x;
            hX_xy += hxy * pos.x; hX_xz += hxz * pos.x; hX_yz += hyz * pos.x;

            hY_xx += hxx * pos.y; hY_yy += hyy * pos.y; hY_zz += hzz * pos.y;
            hY_xy += hxy * pos.y; hY_xz += hxz * pos.y; hY_yz += hyz * pos.y;

            hZ_xx += hxx * pos.z; hZ_yy += hyy * pos.z; hZ_zz += hzz * pos.z;
            hZ_xy += hxy * pos.z; hZ_xz += hxz * pos.z; hZ_yz += hyz * pos.z;
        }

        const CfemReal i00 = out.invJacobian(0,0), i01 = out.invJacobian(0,1), i02 = out.invJacobian(0,2);
        const CfemReal i10 = out.invJacobian(1,0), i11 = out.invJacobian(1,1), i12 = out.invJacobian(1,2);
        const CfemReal i20 = out.invJacobian(2,0), i21 = out.invJacobian(2,1), i22 = out.invJacobian(2,2);

        // Helper lambda pour la transformation de congruence
        auto computeComp = [&](CfemInt p, CfemInt q, CfemReal r_xx, CfemReal r_yy, CfemReal r_zz, 
                                                     CfemReal r_xy, CfemReal r_xz, CfemReal r_yz) {
            CfemReal res = out.invJacobian(0,p)*r_xx*out.invJacobian(0,q) + 
                           out.invJacobian(1,p)*r_yy*out.invJacobian(1,q) + 
                           out.invJacobian(2,p)*r_zz*out.invJacobian(2,q);
            res += (out.invJacobian(0,p)*r_xy*out.invJacobian(1,q) + out.invJacobian(1,p)*r_xy*out.invJacobian(0,q)) +
                   (out.invJacobian(0,p)*r_xz*out.invJacobian(2,q) + out.invJacobian(2,p)*r_xz*out.invJacobian(0,q)) +
                   (out.invJacobian(1,p)*r_yz*out.invJacobian(2,q) + out.invJacobian(2,p)*r_yz*out.invJacobian(1,q));
            return res;
        };

        // Pointeurs SoA pour la lecture des gradients physiques
        const CfemReal* dNdx = out.dN_dx.data();
        const CfemReal* dNdy = out.dN_dy.data();
        const CfemReal* dNdz = out.dN_dz.data();

        // Pointeurs SoA pour l'écriture des Hessiens physiques
        CfemReal* d2x2 = out.d2N_dx2.data();
        CfemReal* d2y2 = out.d2N_dy2.data();
        CfemReal* d2z2 = out.d2N_dz2.data();
        CfemReal* dxdy = out.d2N_dxdy.data();
        CfemReal* dxdz = out.d2N_dxdz.data();
        CfemReal* dydz = out.d2N_dydz.data();

        // --- Étape 2 : Matrice R et Transformation vers l'espace physique ---
        for (CfemInt a = 0; a < n; ++a)
        {
            // R = H_ref - sum( dN/dx_k * H_geom_k )
            CfemReal r_xx = d2xi2[a]      - (dNdx[a] * hX_xx + dNdy[a] * hY_xx + dNdz[a] * hZ_xx);
            CfemReal r_yy = d2eta2[a]     - (dNdx[a] * hX_yy + dNdy[a] * hY_yy + dNdz[a] * hZ_yy);
            CfemReal r_zz = d2zeta2[a]    - (dNdx[a] * hX_zz + dNdy[a] * hY_zz + dNdz[a] * hZ_zz);
            CfemReal r_xy = d2xideta[a]   - (dNdx[a] * hX_xy + dNdy[a] * hY_xy + dNdz[a] * hZ_xy);
            CfemReal r_xz = d2xidzeta[a]  - (dNdx[a] * hX_xz + dNdy[a] * hY_xz + dNdz[a] * hZ_xz);
            CfemReal r_yz = d2etadzeta[a] - (dNdx[a] * hX_yz + dNdy[a] * hY_yz + dNdz[a] * hZ_yz);
            
            // Écriture directe et contiguë des 6 composantes (Congruence)
            d2x2[a] = computeComp(0, 0, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
            d2y2[a] = computeComp(1, 1, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
            d2z2[a] = computeComp(2, 2, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
            dxdy[a] = computeComp(0, 1, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
            dxdz[a] = computeComp(0, 2, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
            dydz[a] = computeComp(1, 2, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
        }
    }

    inline void ReferenceElement::computePhysicalSecondDerivativesBatch(const DynamicVector<Vector3D> &physicalNodesCoords,
                                                                        const ShapeInfoBatch &shape,
                                                                        MappingInfoBatch &out) const
    {
        const CfemInt nQp = shape.numQuadraturePoints;
        const CfemInt n = shape.numNodes;

        for (CfemInt q = 0; q < nQp; ++q)
        {
            // --- Étape 1 : Calcul de la courbure géométrique pour le point q ---
            CfemReal hX_xx = 0, hX_yy = 0, hX_zz = 0, hX_xy = 0, hX_xz = 0, hX_yz = 0;
            CfemReal hY_xx = 0, hY_yy = 0, hY_zz = 0, hY_xy = 0, hY_xz = 0, hY_yz = 0;
            CfemReal hZ_xx = 0, hZ_yy = 0, hZ_zz = 0, hZ_xy = 0, hZ_xz = 0, hZ_yz = 0;

            const CfemReal* d2xi2      = shape.d2N_dxi2[q].data();
            const CfemReal* d2eta2     = shape.d2N_deta2[q].data();
            const CfemReal* d2zeta2    = shape.d2N_dzeta2[q].data();
            const CfemReal* d2xideta   = shape.d2N_dxideta[q].data();
            const CfemReal* d2xidzeta  = shape.d2N_dxidzeta[q].data();
            const CfemReal* d2etadzeta = shape.d2N_detadzeta[q].data();

            for (CfemInt a = 0; a < n; ++a)
            {
                const Vector3D &pos = physicalNodesCoords[a];
                CfemReal hxx = d2xi2[a];
                CfemReal hyy = d2eta2[a];
                CfemReal hzz = d2zeta2[a];
                CfemReal hxy = d2xideta[a];
                CfemReal hxz = d2xidzeta[a];
                CfemReal hyz = d2etadzeta[a];

                hX_xx += hxx * pos.x; hX_yy += hyy * pos.x; hX_zz += hzz * pos.x;
                hX_xy += hxy * pos.x; hX_xz += hxz * pos.x; hX_yz += hyz * pos.x;

                hY_xx += hxx * pos.y; hY_yy += hyy * pos.y; hY_zz += hzz * pos.y;
                hY_xy += hxy * pos.y; hY_xz += hxz * pos.y; hY_yz += hyz * pos.y;

                hZ_xx += hxx * pos.z; hZ_yy += hyy * pos.z; hZ_zz += hzz * pos.z;
                hZ_xy += hxy * pos.z; hZ_xz += hxz * pos.z; hZ_yz += hyz * pos.z;
            }

            // Récupération de l'inverse du Jacobien pour le point q
            const auto& invJ = out.invJacobians[q];

            // Lambda helper local pour la transformation de congruence au point q
            auto computeComp = [&](CfemInt p, CfemInt qIndex, CfemReal r_xx, CfemReal r_yy, CfemReal r_zz, 
                                                              CfemReal r_xy, CfemReal r_xz, CfemReal r_yz) {
                CfemReal res = invJ(0,p)*r_xx*invJ(0,qIndex) + 
                               invJ(1,p)*r_yy*invJ(1,qIndex) + 
                               invJ(2,p)*r_zz*invJ(2,qIndex);
                res += (invJ(0,p)*r_xy*invJ(1,qIndex) + invJ(1,p)*r_xy*invJ(0,qIndex)) +
                       (invJ(0,p)*r_xz*invJ(2,qIndex) + invJ(2,p)*r_xz*invJ(0,qIndex)) +
                       (invJ(1,p)*r_yz*invJ(2,qIndex) + invJ(2,p)*r_yz*invJ(0,qIndex));
                return res;
            };

            // Recouvrement des pointeurs de lignes physiques (Lecture des gradients, Écriture des Hessiens)
            const CfemReal* dNdx = out.dN_dx[q].data();
            const CfemReal* dNdy = out.dN_dy[q].data();
            const CfemReal* dNdz = out.dN_dz[q].data();

            CfemReal* d2x2 = out.d2N_dx2[q].data();
            CfemReal* d2y2 = out.d2N_dy2[q].data();
            CfemReal* d2z2 = out.d2N_dz2[q].data();
            CfemReal* dxdy = out.d2N_dxdy[q].data();
            CfemReal* dxdz = out.d2N_dxdz[q].data();
            CfemReal* dydz = out.d2N_dydz[q].data();

            // Étape 2 : Calcul de R et application de la congruence pour chaque nœud
            for (CfemInt a = 0; a < n; ++a)
            {
                CfemReal r_xx = d2xi2[a]      - (dNdx[a] * hX_xx + dNdy[a] * hY_xx + dNdz[a] * hZ_xx);
                CfemReal r_yy = d2eta2[a]     - (dNdx[a] * hX_yy + dNdy[a] * hY_yy + dNdz[a] * hZ_yy);
                CfemReal r_zz = d2zeta2[a]    - (dNdx[a] * hX_zz + dNdy[a] * hY_zz + dNdz[a] * hZ_zz);
                CfemReal r_xy = d2xideta[a]   - (dNdx[a] * hX_xy + dNdy[a] * hY_xy + dNdz[a] * hZ_xy);
                CfemReal r_xz = d2xidzeta[a]  - (dNdx[a] * hX_xz + dNdy[a] * hY_xz + dNdz[a] * hZ_xz);
                CfemReal r_yz = d2etadzeta[a] - (dNdx[a] * hX_yz + dNdy[a] * hY_yz + dNdz[a] * hZ_yz);
                
                // Enregistrement direct dans les structures matricielles contigues du Batch
                d2x2[a] = computeComp(0, 0, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
                d2y2[a] = computeComp(1, 1, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
                d2z2[a] = computeComp(2, 2, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
                dxdy[a] = computeComp(0, 1, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
                dxdz[a] = computeComp(0, 2, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
                dydz[a] = computeComp(1, 2, r_xx, r_yy, r_zz, r_xy, r_xz, r_yz);
            }
        }
    }


    inline void ReferenceElement::computeGeometricNormal(const ShapeInfo &shape, MappingInfo &out) const
    {
        if (m_dimension == 2) 
        {
            // Variété 2D (Surface, Triangle/Quad) dans l'espace 3D
            // Les deux vecteurs tangents sont les deux premières colonnes du Jacobien
            Vector3D t1(out.jacobian(0, 0), out.jacobian(1, 0), out.jacobian(2, 0));
            Vector3D t2(out.jacobian(0, 1), out.jacobian(1, 1), out.jacobian(2, 1));
            
            out.normal = (t1.cross(t2)).normalized();

        }
        else if (m_dimension == 1) 
        {
            // Variété 1D (Ligne) dans un espace 2D (ou 3D z=0)
            CfemReal dx = out.jacobian(0, 0);
            CfemReal dy = out.jacobian(1, 0);
            CfemReal norm = std::sqrt(dx * dx + dy * dy);
            
            if (norm > cfem::ZERO_TOLERANCE) {
                // Règle de la main droite (rotation de -90 degrés du vecteur tangent)
                // L'orientation extérieure dépend de la convention de rotation de ton maillage
                out.normal = Vector3D(dy / norm, -dx / norm, 0.0); 
            } else {
                out.normal = Vector3D(0.0, 0.0, 0.0);
            }
        }
        else 
        {
            // Dimension 3 (Volume). Geometric normal undefined inside a 3D cell.
            out.normal = Vector3D(0.0, 0.0, 0.0);
        }
    }

    inline void ReferenceElement::computeGeometricNormalBatch(const ShapeInfoBatch &shapeBatch, 
                                                             MappingInfoBatch &outBatch) const
    {
        const CfemInt nQp = shapeBatch.numQuadraturePoints;
        CFEM_ASSERT(outBatch.jacobians.size() == nQp, "ReferenceElement::computeGeometricNormalBatch outBatch.jacobians not pre-computed?");
        if (m_dimension == 2) 
        {
            // Variété 2D (Surface, Triangle/Quad) dans l'espace 3D
            // Les deux vecteurs tangents sont les deux premières colonnes du Jacobien
            for (CfemInt q = 0; q < nQp; ++q){
                Vector3D t1(outBatch.jacobians[q].xx(), outBatch.jacobians[q].yx(), outBatch.jacobians[q].zx());
                Vector3D t2(outBatch.jacobians[q].xy(), outBatch.jacobians[q].yy(), outBatch.jacobians[q].zy());
                
                outBatch.normals[q] = (t1.cross(t2)).normalized();
            }

        }
        else if (m_dimension == 1) 
        {
            for (CfemInt q = 0; q < nQp; ++q){
                // Variété 1D (Ligne) dans un espace 2D (ou 3D z=0)
                CfemReal dx = outBatch.jacobians[q].xx();
                CfemReal dy = outBatch.jacobians[q].yx();
                CfemReal norm = std::sqrt(dx * dx + dy * dy);
                
                if (norm > cfem::ZERO_TOLERANCE) {
                    // Règle de la main droite (rotation de -90 degrés du vecteur tangent)
                    // L'orientation extérieure dépend de la convention de rotation de ton maillage
                    outBatch.normals[q] = Vector3D(dy / norm, -dx / norm, 0.0); 
                } else {
                    outBatch.normals[q] = Vector3D(0.0, 0.0, 0.0);
                }
            }
        }
        else 
        {
            // Dimension 3 (Volume). Geometric normal undefined inside a 3D cell.
            for (CfemInt q = 0; q < nQp; ++q){
                outBatch.normals[q] = Vector3D(0.0, 0.0, 0.0);
            }
        }
    }
    // ----------------------------------------------------------------------
    // ------------------------------- interpolation ------------------------
    template<typename T>
    inline T ReferenceElement::interpolate(const ShapeInfo& shape, 
                                           const DynamicVector<T>& nodalValues) const
    {
        const CfemInt n = shape.numNodes;
        CFEM_ASSERT(shape.has(ShapeEvaluationFlags::Value), "Shape values not computed");
        CFEM_ASSERT(nodalValues.size() == n, "Nodal values size mismatch");

        T result{};
        const CfemReal* N = shape.values.data();

        for (CfemInt i = 0; i < n; ++i){
            result += nodalValues[i] * N[i];
        }
        return result;
    }

    inline Vector3D ReferenceElement::interpolateGradient(const MappingInfo &map,
                                                          const DynamicVector<CfemReal> &nodalValues) const
    {
        const CfemInt n = map.numNodes;
        CFEM_ASSERT(hasFlag(map.validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), 
                    "Physical gradients not computed in MappingInfo");
        CFEM_ASSERT(nodalValues.size() == n, "Nodal values size mismatch");

        // Accumulateurs dans les registres CPU (ultra rapide)
        CfemReal gradX = 0.0;
        CfemReal gradY = 0.0;
        CfemReal gradZ = 0.0;

        // Pointeurs vers la mémoire SoA
        const CfemReal* dx = map.dN_dx.data();
        const CfemReal* dy = map.dN_dy.data();
        const CfemReal* dz = map.dN_dz.data();

        for (CfemInt i = 0; i < n; ++i) {
            const CfemReal u = nodalValues[i];
            gradX += dx[i] * u;
            gradY += dy[i] * u;
            gradZ += dz[i] * u;
        }

        return Vector3D(gradX, gradY, gradZ);
    }

    inline SymTensor3D ReferenceElement::interpolateHessian(const MappingInfo &map,
                                                            const DynamicVector<CfemReal> &nodalValues) const
    {
        const CfemInt n = map.numNodes;
        CFEM_ASSERT(hasFlag(map.validFlags, MappingEvaluationFlags::PhysicalSecondDerivatives), 
                    "Physical hessians not computed in MappingInfo");
        CFEM_ASSERT(nodalValues.size() == n, "Nodal values size mismatch");

        CfemReal h11 = 0.0, h22 = 0.0, h33 = 0.0;
        CfemReal h12 = 0.0, h13 = 0.0, h23 = 0.0;

        const CfemReal* d2x  = map.d2N_dx2.data();
        const CfemReal* d2y  = map.d2N_dy2.data();
        const CfemReal* d2z  = map.d2N_dz2.data();
        const CfemReal* dxdy = map.d2N_dxdy.data();
        const CfemReal* dxdz = map.d2N_dxdz.data();
        const CfemReal* dydz = map.d2N_dydz.data();

        for (CfemInt i = 0; i < n; ++i)
        {
            const CfemReal u = nodalValues[i];
            h11 += d2x[i] * u;
            h22 += d2y[i] * u;
            h33 += d2z[i] * u;
            h12 += dxdy[i] * u;
            h13 += dxdz[i] * u;
            h23 += dydz[i] * u;
        }

        return SymTensor3D(h11, h22, h33, h23, h13, h12);
    }

}  // namespace cfem