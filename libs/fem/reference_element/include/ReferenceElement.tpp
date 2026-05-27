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

    inline void ReferenceElement::evaluateShapes(const Vector3D &xi, ShapeEvaluationFlags requestedShape, ShapeInfo &shape) const
    {
        const bool includeSecondDerivatives = hasFlag(requestedShape, ShapeEvaluationFlags::SecondDerivatives);
        shape.setup(m_numNodes, includeSecondDerivatives);

        if (hasFlag(requestedShape, ShapeEvaluationFlags::Value))
        {
            evaluateShapesValues(xi, shape);
        }
        if (hasFlag(requestedShape, ShapeEvaluationFlags::FirstDerivatives))
        {
            evaluateShapesDerivatives(xi, shape);
        }
        if (includeSecondDerivatives)
        {
            evaluateShapesSecondDerivatives(xi, shape);
        }

        shape.validFlags |= requestedShape;
    }

    inline void ReferenceElement::evaluateShapesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeEvaluationFlags requestedShape, ShapeInfoBatch &shapeBatch) const
    {
        const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
        const bool includeSecondDerivatives = hasFlag(requestedShape, ShapeEvaluationFlags::SecondDerivatives);
        shapeBatch.setup(nQp, m_numNodes, includeSecondDerivatives);

        if (hasFlag(requestedShape, ShapeEvaluationFlags::Value))
        {
            evaluateShapesValuesBatch(xiPoints, shapeBatch);
        }
        if (hasFlag(requestedShape, ShapeEvaluationFlags::FirstDerivatives))
        {
            evaluateShapesDerivativesBatch(xiPoints, shapeBatch);
        }
        if (includeSecondDerivatives)
        {
            evaluateShapesSecondDerivativesBatch(xiPoints, shapeBatch);
        }

        shapeBatch.validFlags |= requestedShape;
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
        if (shape.shapeMustBeUpdated)
        {
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
        out.setup(m_numNodes, hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalSecondDerivatives));

        if (hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalPosition))
        {
            CfemReal px = 0.0, py = 0.0, pz = 0.0;
            const CfemReal *N = shape.values();
            for (CfemInt i = 0; i < m_numNodes; ++i)
            {
                const Vector3D &pos = physicalNodesCoords[i];
                const CfemReal ni = N[i];
                px += pos.x * ni;
                py += pos.y * ni;
                pz += pos.z * ni;
            }
            out.physicalPosition = Vector3D(px, py, pz);
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::Jacobian))
        {
            computeJacobian(shape, physicalNodesCoords, out);
            out.validFlags |= MappingEvaluationFlags::Jacobian;
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::Determinant))
        {
            computeDeterminant(out);
            out.validFlags |= MappingEvaluationFlags::Determinant;
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::InverseJacobian))
        {
            computeInverseJacobian(out);
            out.validFlags |= MappingEvaluationFlags::InverseJacobian;
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalFirstDerivatives))
        {
            computePhysicalFirstDerivatives(shape, out);
            out.validFlags |= MappingEvaluationFlags::PhysicalFirstDerivatives;
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalSecondDerivatives))
        {
            computePhysicalSecondDerivatives(physicalNodesCoords, shape, out);
            out.validFlags |= MappingEvaluationFlags::PhysicalSecondDerivatives;
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::Normal))
        {
            computeGeometricNormal(out);
            out.validFlags |= MappingEvaluationFlags::Normal;
        }

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
        if (shapeBatch.shapeMustBeUpdated)
        {
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

        if (hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalPosition))
        {
            for (CfemInt q = 0; q < nQp; ++q)
            {
                CfemReal px = 0.0, py = 0.0, pz = 0.0;

                const CfemReal *N = shapeBatch.values(q);

                for (CfemInt i = 0; i < m_numNodes; ++i)
                {
                    const Vector3D &pos = physicalNodesCoords[i];
                    const CfemReal ni = N[i];
                    px += pos.x * ni;
                    py += pos.y * ni;
                    pz += pos.z * ni;
                }
                outBatch.physicalPositions[q] = Vector3D(px, py, pz); // Écriture unique
            }
        }

        // Évaluation des opérateurs physiques par lot (Chaque fonction traite l'ensemble du vecteur de points qp)
        if (hasFlag(requestedMapping, MappingEvaluationFlags::Jacobian))
        {
            computeJacobianBatch(shapeBatch, physicalNodesCoords, outBatch);
            outBatch.validFlags |= MappingEvaluationFlags::Jacobian;
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::Determinant))
        {

            computeDeterminantBatch(outBatch);
            outBatch.validFlags |= MappingEvaluationFlags::Determinant;
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::InverseJacobian))
        {
            computeInverseJacobianBatch(outBatch);
            outBatch.validFlags |= MappingEvaluationFlags::InverseJacobian;
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalFirstDerivatives))
        {
            computePhysicalFirstDerivativesBatch(shapeBatch, outBatch);
            outBatch.validFlags |= MappingEvaluationFlags::PhysicalFirstDerivatives;
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::PhysicalSecondDerivatives))
        {
            computePhysicalSecondDerivativesBatch(physicalNodesCoords, shapeBatch, outBatch);
            outBatch.validFlags |= MappingEvaluationFlags::PhysicalSecondDerivatives;
        }

        if (hasFlag(requestedMapping, MappingEvaluationFlags::Normal))
        {
            computeGeometricNormalBatch(outBatch);
            outBatch.validFlags |= MappingEvaluationFlags::Normal;
        }

    }

    inline void ReferenceElement::computeJacobian(const ShapeInfo &shape,
                                                  const DynamicVector<Vector3D> &physicalNodesCoords,
                                                  MappingInfo &out) const
    {
        const CfemInt n = shape.numNodes;

        CfemReal j00 = 0.0, j01 = 0.0, j02 = 0.0;
        CfemReal j10 = 0.0, j11 = 0.0, j12 = 0.0;
        CfemReal j20 = 0.0, j21 = 0.0, j22 = 0.0;

        // Pointeurs SoA pour ce point de quadrature précis
        const CfemReal *dxi = shape.dN_dxi();
        const CfemReal *deta = shape.dN_deta();
        const CfemReal *dzeta = shape.dN_dzeta();

        for (CfemInt i = 0; i < n; ++i)
        {
            const Vector3D &X = physicalNodesCoords[i];
            j00 += X.x * dxi[i];
            j10 += X.y * dxi[i];
            j20 += X.z * dxi[i];
            j01 += X.x * deta[i];
            j11 += X.y * deta[i];
            j21 += X.z * deta[i];
            j02 += X.x * dzeta[i];
            j12 += X.y * dzeta[i];
            j22 += X.z * dzeta[i];
        }

        out.jacobian.data[Tensor3D::idx_xx] = j00;
        out.jacobian.data[Tensor3D::idx_xy] = j01;
        out.jacobian.data[Tensor3D::idx_xz] = j02;
        out.jacobian.data[Tensor3D::idx_yx] = j10;
        out.jacobian.data[Tensor3D::idx_yy] = j11;
        out.jacobian.data[Tensor3D::idx_yz] = j12;
        out.jacobian.data[Tensor3D::idx_zx] = j20;
        out.jacobian.data[Tensor3D::idx_zy] = j21;
        out.jacobian.data[Tensor3D::idx_zz] = j22;
    }

    inline void ReferenceElement::computeJacobianBatch(const ShapeInfoBatch &shape,
                                                       const DynamicVector<Vector3D> &physicalNodesCoords,
                                                       MappingInfoBatch &out) const
    {
        CFEM_ASSERT(shape.numQuadraturePoints == out.numQuadraturePoints && shape.has(ShapeEvaluationFlags::FirstDerivatives), "Shapes derivatives not precomputed?");
        const CfemInt nQp = shape.numQuadraturePoints;
        const CfemInt n = shape.numNodes;

        const CfemInt qEnd = out.isCellAffine ? 1 : nQp;

        for (CfemInt q = 0; q < qEnd; ++q)
        {
            CfemReal j00 = 0.0, j01 = 0.0, j02 = 0.0;
            CfemReal j10 = 0.0, j11 = 0.0, j12 = 0.0;
            CfemReal j20 = 0.0, j21 = 0.0, j22 = 0.0;

            // Pointeurs SoA pour ce point de quadrature précis
            const CfemReal *dxi = shape.dN_dxi(q);
            const CfemReal *deta = shape.dN_deta(q);
            const CfemReal *dzeta = shape.dN_dzeta(q);

            for (CfemInt i = 0; i < n; ++i)
            {
                const Vector3D &X = physicalNodesCoords[i];
                j00 += X.x * dxi[i];
                j10 += X.y * dxi[i];
                j20 += X.z * dxi[i];
                j01 += X.x * deta[i];
                j11 += X.y * deta[i];
                j21 += X.z * deta[i];
                j02 += X.x * dzeta[i];
                j12 += X.y * dzeta[i];
                j22 += X.z * dzeta[i];
            }

            out.jacobians[q].xx() = j00;
            out.jacobians[q].xy() = j01;
            out.jacobians[q].xz() = j02;
            out.jacobians[q].yx() = j10;
            out.jacobians[q].yy() = j11;
            out.jacobians[q].yz() = j12;
            out.jacobians[q].zx() = j20;
            out.jacobians[q].zy() = j21;
            out.jacobians[q].zz() = j22;
        }

        if (out.isCellAffine)
        {
            for (CfemInt q = 1; q < nQp; ++q)
            {
                out.jacobians[q] = out.jacobians[0];
            }
        }
    }

    inline void ReferenceElement::computeDeterminant(MappingInfo &out) const
    {
        if (m_dimension == 3)
        {
            out.detJ = det(out.jacobian);
        }
        else if (m_dimension == 2)
        {
            // Métrique de surface
            const CfemReal j00 = out.jacobian.data[Tensor3D::idx_xx], j01 = out.jacobian.data[Tensor3D::idx_xy];
            const CfemReal j10 = out.jacobian.data[Tensor3D::idx_yx], j11 = out.jacobian.data[Tensor3D::idx_yy];
            const CfemReal j20 = out.jacobian.data[Tensor3D::idx_zx], j21 = out.jacobian.data[Tensor3D::idx_zy];

            CfemReal g11 = j00 * j00 + j10 * j10 + j20 * j20;
            CfemReal g22 = j01 * j01 + j11 * j11 + j21 * j21;
            CfemReal g12 = j00 * j01 + j10 * j11 + j20 * j21;

            out.metric2D = SymTensor2D(g11, g22, g12);
            out.detJ = std::sqrt(std::max<CfemReal>(0.0, g11 * g22 - g12 * g12));
        }
        else if (m_dimension == 1)
        {
            const CfemReal j00 = out.jacobian.data[Tensor3D::idx_xx];
            const CfemReal j10 = out.jacobian.data[Tensor3D::idx_yx];
            const CfemReal j20 = out.jacobian.data[Tensor3D::idx_zx];
            out.metric1D = j00 * j00 + j10 * j10 + j20 * j20;
            out.detJ = std::sqrt(out.metric1D);
        }
    }

    inline void ReferenceElement::computeDeterminantBatch(MappingInfoBatch &out) const
    {
        

        const CfemInt nQp = out.numQuadraturePoints;
        const CfemInt qEnd = out.isCellAffine ? 1 : nQp;

        if (m_dimension == 3)
        {
            for (CfemInt q = 0; q < qEnd; ++q)
            {
                out.detJs[q] = det(out.jacobians[q]);
            }

            if (out.isCellAffine)
            {
                for (CfemInt q = 1; q < nQp; ++q)
                    out.detJs[q] = out.detJs[0];
            }
        }
        else if (m_dimension == 2)
        {
            // Métrique de surface
            for (CfemInt q = 0; q < qEnd; ++q)
            {
                const CfemReal j00 = out.jacobians[q].data[Tensor3D::idx_xx], j01 = out.jacobians[q].data[Tensor3D::idx_xy];
                const CfemReal j10 = out.jacobians[q].data[Tensor3D::idx_yx], j11 = out.jacobians[q].data[Tensor3D::idx_yy];
                const CfemReal j20 = out.jacobians[q].data[Tensor3D::idx_zx], j21 = out.jacobians[q].data[Tensor3D::idx_zy];

                CfemReal g11 = j00 * j00 + j10 * j10 + j20 * j20;
                CfemReal g22 = j01 * j01 + j11 * j11 + j21 * j21;
                CfemReal g12 = j00 * j01 + j10 * j11 + j20 * j21;

                out.metric2Ds[q] = SymTensor2D{g11, g22, g12};
                out.detJs[q] = std::sqrt(std::max<CfemReal>(0.0, g11 * g22 - g12 * g12));
            }

            if (out.isCellAffine)
            {
                for (CfemInt q = 1; q < nQp; ++q)
                {
                    out.metric2Ds[q] = out.metric2Ds[0];
                    out.detJs[q] = out.detJs[0];
                }
            }
        }
        else if (m_dimension == 1)
        {
            for (CfemInt q = 0; q < qEnd; ++q)
            {
                const CfemReal j00 = out.jacobians[q].data[Tensor3D::idx_xx];
                const CfemReal j10 = out.jacobians[q].data[Tensor3D::idx_yx];
                const CfemReal j20 = out.jacobians[q].data[Tensor3D::idx_zx];
                CfemReal g11 = j00 * j00 + j10 * j10 + j20 * j20;

                out.metric1Ds[q] = g11;
                out.detJs[q] = std::sqrt(g11);
            }

            if (out.isCellAffine)
            {
                for (CfemInt q = 1; q < nQp; ++q)
                {
                    out.metric1Ds[q] = out.metric1Ds[0];
                    out.detJs[q] = out.detJs[0];
                }
            }
        }
    }

    /**
     * @brief Computes the contravariant geometric transformation operator G.
     *
     * Maps reference-space derivatives to physical-space derivatives:
     *      \$ \nabla_x(u) = G^T \nabla_{\xi}(u) \$
     *
     * G is defined as the Moore-Penrose pseudo-inverse of the Jacobian J,
     * depending on the manifold dimension:
     *  - 3D (Volume)  : \$ G = J^{-1} \$
     *  - 2D (Surface) : \$ G = (J^T J)^{-1} J^T \$
     *  - 1D (Curve)   : \$ G = J^T / ||J||^2 \$
     *
     * @pre `out.jacobian`, `out.detJ`, and the appropriate metric tensor
     *      (`metric2D` or `metric1D`) must be precomputed.
     * @note Stored in `out.invJacobian` for historical compatibility.
     *
     * @param[in]  shape               Element shape information (unused).
     * @param[in]  physicalNodesCoords Physical node coordinates (unused).
     * @param[out] out                 Mapping structure receiving the operator.
     */
    inline void ReferenceElement::computeInverseJacobian(MappingInfo &out) const
    {
        CFEM_ASSERT(std::abs(out.detJ) > cfem::ZERO_TOLERANCE, "Degenerate Jacobian detected");

        Tensor3D &G = out.invJacobian;

        // ================================================================
        // 3D volumetric element  G = J^{-1}
        // ================================================================
        if (m_dimension == 3)
        {
            G = inverse(out.jacobian, out.detJ);
            return;
        }

        // ================================================================
        // 2D embedded manifold G = (J^T J)^{-1} J^T
        // ================================================================
        if (m_dimension == 2)
        {
            const CfemReal g11 = out.metric2D.xx();
            const CfemReal g22 = out.metric2D.yy();
            const CfemReal g12 = out.metric2D.xy();

            // det(G_metric) = det(J)^2
            const CfemReal invDetMetric = 1.0 / (out.detJ * out.detJ);

            // Inverse metric tensor coefficients
            const CfemReal m11 = g22 * invDetMetric;
            const CfemReal m22 = g11 * invDetMetric;
            const CfemReal m12 = -g12 * invDetMetric;

            // Jacobian columns (surface tangents)
            const CfemReal j00 = out.jacobian.xx();
            const CfemReal j10 = out.jacobian.yx();
            const CfemReal j20 = out.jacobian.zx();

            const CfemReal j01 = out.jacobian.xy();
            const CfemReal j11 = out.jacobian.yy();
            const CfemReal j21 = out.jacobian.zy();

            // First contravariant basis vector
            G.xx() = m11 * j00 + m12 * j01;
            G.xy() = m11 * j10 + m12 * j11;
            G.xz() = m11 * j20 + m12 * j21;

            // Second contravariant basis vector
            G.yx() = m12 * j00 + m22 * j01;
            G.yy() = m12 * j10 + m22 * j11;
            G.yz() = m12 * j20 + m22 * j21;

            // Third row unused for 2D manifolds
            G.zx() = 0.0;
            G.zy() = 0.0;
            G.zz() = 0.0;

            return;
        }

        // ================================================================
        // 1D embedded manifold  G = J^T / ||J||^2
        // ================================================================
        if (m_dimension == 1)
        {
            const CfemReal invMetric = 1.0 / out.metric1D;

            G(0, 0) = out.jacobian(0, 0) * invMetric;
            G(0, 1) = out.jacobian(1, 0) * invMetric;
            G(0, 2) = out.jacobian(2, 0) * invMetric;

            G(1, 0) = 0.0;
            G(1, 1) = 0.0;
            G(1, 2) = 0.0;

            G(2, 0) = 0.0;
            G(2, 1) = 0.0;
            G(2, 2) = 0.0;

            return;
        }

        CFEM_ASSERT(false, "Unsupported manifold dimension");
    }

    inline void ReferenceElement::computeInverseJacobianBatch(MappingInfoBatch &out) const
    {
        const CfemInt nQp = out.numQuadraturePoints;
        const CfemInt qEnd = out.isCellAffine ? 1 : nQp;

        if (m_dimension == 3)
        {
            for (CfemInt q = 0; q < qEnd; ++q)
            {
                CFEM_ASSERT(std::abs(out.detJs[q]) > cfem::ZERO_TOLERANCE, "Degenerate Jacobian detected");
                out.invJacobians[q] = inverse(out.jacobians[q], out.detJs[q]);
            }
        }
        else if (m_dimension == 2)
        {
            for (CfemInt q = 0; q < qEnd; ++q)
            {
                CFEM_ASSERT(std::abs(out.detJs[q]) > cfem::ZERO_TOLERANCE, "Degenerate Jacobian detected");
                const CfemReal g11 = out.metric2Ds[q].xx();
                const CfemReal g22 = out.metric2Ds[q].yy();
                const CfemReal g12 = out.metric2Ds[q].xy();

                // OPTIMISATION : Réutilisation du déterminant (sauve 2 multiplications, 1 soustraction)
                const CfemReal invDetMetric = 1.0 / (out.detJs[q] * out.detJs[q]);

                const CfemReal m11 = g22 * invDetMetric;
                const CfemReal m22 = g11 * invDetMetric;
                const CfemReal m12 = -g12 * invDetMetric;

                const auto &jac_q = out.jacobians[q];
                const CfemReal j00 = jac_q.xx(), j10 = jac_q.yx(), j20 = jac_q.zx();
                const CfemReal j01 = jac_q.xy(), j11 = jac_q.yy(), j21 = jac_q.zy();

                // OPTIMISATION : Déroulage manuel total (zéro boucle)
                auto &G = out.invJacobians[q];
                G.xx() = m11 * j00 + m12 * j01;
                G.xy() = m11 * j10 + m12 * j11;
                G.xz() = m11 * j20 + m12 * j21;
                G.yx() = m12 * j00 + m22 * j01;
                G.yy() = m12 * j10 + m22 * j11;
                G.yz() = m12 * j20 + m22 * j21;
                G.zx() = 0.0;
                G.zy() = 0.0;
                G.zz() = 0.0;
            }
        }
        else if (m_dimension == 1)
        {
            for (CfemInt q = 0; q < qEnd; ++q)
            {
                CFEM_ASSERT(std::abs(out.detJs[q]) > cfem::ZERO_TOLERANCE, "Degenerate Jacobian detected");
                const CfemReal invMetric = 1.0 / out.metric1Ds[q];
                const auto &jac_q = out.jacobians[q];
                auto &G = out.invJacobians[q];

                G.xx() = jac_q.xx() * invMetric;
                G.xy() = jac_q.yx() * invMetric;
                G.xz() = jac_q.zx() * invMetric;
                G.yx() = 0.0;
                G.yy() = 0.0;
                G.yz() = 0.0;
                G.zx() = 0.0;
                G.zy() = 0.0;
                G.zz() = 0.0;
            }
        }

        if (out.isCellAffine)
        {
            for (CfemInt q = 1; q < nQp; ++q)
            {
                out.invJacobians[q] = out.invJacobians[0];
            }
        }
    }

    inline void ReferenceElement::computePhysicalFirstDerivatives(const ShapeInfo &shape,
                                                                  MappingInfo &out) const
    {
        CFEM_ASSERT(shape.has(ShapeEvaluationFlags::FirstDerivatives), "Shapes derivatives not pre-computed");
        CFEM_ASSERT(out.has(MappingEvaluationFlags::InverseJacobian), "Inverse jacobian not pre-computed");
        const CfemInt n = shape.numNodes;

        // J^{-T} appliqué au vecteur local.
        const CfemReal i00 = out.invJacobian.xx(), i01 = out.invJacobian.xy(), i02 = out.invJacobian.xz();
        const CfemReal i10 = out.invJacobian.yx(), i11 = out.invJacobian.yy(), i12 = out.invJacobian.yz();
        const CfemReal i20 = out.invJacobian.zx(), i21 = out.invJacobian.zy(), i22 = out.invJacobian.zz();

        // Pointeurs SoA (lecture)
        const CfemReal *__restrict dxi = shape.dN_dxi();
        const CfemReal *__restrict deta = shape.dN_deta();
        const CfemReal *__restrict dzeta = shape.dN_dzeta();

        // Pointeurs SoA (écriture)
        CfemReal *__restrict dx = out.dN_dx();
        CfemReal *__restrict dy = out.dN_dy();
        CfemReal *__restrict dz = out.dN_dz();

        // Boucle SIMD ultra-rapide
        for (CfemInt a = 0; a < n; ++a)
        {
            const CfemReal gx = dxi[a];
            const CfemReal gy = deta[a];
            const CfemReal gz = dzeta[a];

            dx[a] = i00 * gx + i10 * gy + i20 * gz;
            dy[a] = i01 * gx + i11 * gy + i21 * gz;
            dz[a] = i02 * gx + i12 * gy + i22 * gz;
        }

        out.validFlags |= MappingEvaluationFlags::PhysicalFirstDerivatives;
    }

    /**
     * @brief Computes physical gradients for the entire batch.
     * @note Hyper-optimized utilizing SIMD friendly contiguous memory rows. 
     * Correctly handles embedded manifolds (e.g., 1D/2D elements in 3D ambient space).
     */
    inline void ReferenceElement::computePhysicalFirstDerivativesBatch(const ShapeInfoBatch &shape,
                                                                       MappingInfoBatchView &out) const
    {
        CFEM_ASSERT(shape.has(ShapeEvaluationFlags::FirstDerivatives), "Shapes derivatives not pre-computed");

        const CfemInt nQp = out.numQuadraturePoints;
        const CfemInt n = out.numNodes;
        
        // La dimension de référence (l'élément lui-même : 1 pour Ligne, 2 pour Triangle/Quad, 3 pour Tet/Hex)
        const CfemInt refDim = m_dimension; 
        
        // La dimension de l'espace ambiant déduite intelligemment grâce aux allocations de l'Arena
        const CfemInt spatialDim = (out.dN_dz_data != nullptr) ? 3 : ((out.dN_dy_data != nullptr) ? 2 : 1);

        CfemReal i00 = 0, i01 = 0, i02 = 0, i10 = 0, i11 = 0, i12 = 0, i20 = 0, i21 = 0, i22 = 0;
        
        if (out.isCellAffine && nQp > 0)
        {
            i00 = out.invJacobians[0].xx(); i01 = out.invJacobians[0].xy(); i02 = out.invJacobians[0].xz();
            i10 = out.invJacobians[0].yx(); i11 = out.invJacobians[0].yy(); i12 = out.invJacobians[0].yz();
            i20 = out.invJacobians[0].zx(); i21 = out.invJacobians[0].zy(); i22 = out.invJacobians[0].zz();
        }

        for (CfemInt q = 0; q < nQp; ++q)
        {
            if (!out.isCellAffine)
            {
                i00 = out.invJacobians[q].xx(); i01 = out.invJacobians[q].xy(); i02 = out.invJacobians[q].xz();
                i10 = out.invJacobians[q].yx(); i11 = out.invJacobians[q].yy(); i12 = out.invJacobians[q].yz();
                i20 = out.invJacobians[q].zx(); i21 = out.invJacobians[q].zy(); i22 = out.invJacobians[q].zz();
            }

            // --- LECTURE (Basée sur refDim) ---
            const CfemReal *__restrict dxi = shape.dN_dxi(q);
            const CfemReal *__restrict deta = (refDim >= 2) ? shape.dN_deta(q) : nullptr;
            const CfemReal *__restrict dzeta = (refDim == 3) ? shape.dN_dzeta(q) : nullptr;

            // --- ÉCRITURE (Basée sur spatialDim) ---
            CfemReal *__restrict dx = out.dN_dx_mut(q);
            CfemReal *__restrict dy = (spatialDim >= 2) ? out.dN_dy_mut(q) : nullptr;
            CfemReal *__restrict dz = (spatialDim == 3) ? out.dN_dz_mut(q) : nullptr;

            for (CfemInt a = 0; a < n; ++a)
            {
                const CfemReal gx = dxi[a];
                const CfemReal gy = (refDim >= 2) ? deta[a] : 0.0;
                const CfemReal gz = (refDim == 3) ? dzeta[a] : 0.0;

                // Calcul matriciel optimisé 
                // Même si refDim=2 (gz=0), si spatialDim=3, on écrit bien dans dz !
                dx[a] = i00 * gx + i10 * gy + i20 * gz;
                if (spatialDim >= 2) dy[a] = i01 * gx + i11 * gy + i21 * gz;
                if (spatialDim == 3) dz[a] = i02 * gx + i12 * gy + i22 * gz;
            }
        }

        out.validFlags |= MappingEvaluationFlags::PhysicalFirstDerivatives;
    }

    /**
     * @brief Computes physical gradients for the entire batch.
     * @note Hyper-optimized utilizing SIMD friendly contiguous memory rows.
     */
    inline void ReferenceElement::computePhysicalFirstDerivativesBatch(const ShapeInfoBatch &shape,
                                                                       MappingInfoBatch &out) const
    {
        CFEM_ASSERT(shape.has(ShapeEvaluationFlags::FirstDerivatives), "Shapes derivatives not pre-computed");
        CFEM_ASSERT(out.has(MappingEvaluationFlags::InverseJacobian), "Inverse jacobian not pre-computed");

        const CfemInt nQp = shape.numQuadraturePoints;
        const CfemInt n = shape.numNodes;

        CfemReal i00, i01, i02, i10, i11, i12, i20, i21, i22;
        if (out.isCellAffine && out.numQuadraturePoints > 0)
        {
            i00 = out.invJacobians[0].xx();
            i01 = out.invJacobians[0].xy();
            i02 = out.invJacobians[0].xz();
            i10 = out.invJacobians[0].yx();
            i11 = out.invJacobians[0].yy();
            i12 = out.invJacobians[0].yz();
            i20 = out.invJacobians[0].zx();
            i21 = out.invJacobians[0].zy();
            i22 = out.invJacobians[0].zz();
        }

        for (CfemInt q = 0; q < out.numQuadraturePoints; ++q)
        {
            // Si la cellule n'est PAS affine, on charge le Jacobien spécifique à ce QP
            if (!out.isCellAffine)
            {
                i00 = out.invJacobians[q].xx();
                i01 = out.invJacobians[q].xy();
                i02 = out.invJacobians[q].xz();
                i10 = out.invJacobians[q].yx();
                i11 = out.invJacobians[q].yy();
                i12 = out.invJacobians[q].yz();
                i20 = out.invJacobians[q].zx();
                i21 = out.invJacobians[q].zy();
                i22 = out.invJacobians[q].zz();
            }

            // Pointeurs SoA vers les lignes de la matrice de référence (Lecture)
            const CfemReal *__restrict dxi = shape.dN_dxi(q);
            const CfemReal *__restrict deta = shape.dN_deta(q);
            const CfemReal *__restrict dzeta = shape.dN_dzeta(q);

            // Pointeurs SoA vers les lignes de la matrice physique (Écriture)
            CfemReal *__restrict dx = out.dN_dx(q);
            CfemReal *__restrict dy = out.dN_dy(q);
            CfemReal *__restrict dz = out.dN_dz(q);

            // Boucle sur les nœuds interne : vectorisation SIMD maximale
            for (CfemInt a = 0; a < n; ++a)
            {
                const CfemReal gx = dxi[a];
                const CfemReal gy = deta[a];
                const CfemReal gz = dzeta[a];

                dx[a] = i00 * gx + i10 * gy + i20 * gz;
                dy[a] = i01 * gx + i11 * gy + i21 * gz;
                dz[a] = i02 * gx + i12 * gy + i22 * gz;
            }
        }

        out.validFlags |= MappingEvaluationFlags::PhysicalFirstDerivatives;
    }

    inline void ReferenceElement::computePhysicalSecondDerivatives(const DynamicVector<Vector3D> &physicalNodesCoords,
                                                                   const ShapeInfo &shape,
                                                                   MappingInfo &out) const
    {
        CFEM_ASSERT(shape.has(ShapeEvaluationFlags::FirstDerivatives), "Shapes derivatives not pre-computed");
        CFEM_ASSERT(out.has(MappingEvaluationFlags::InverseJacobian), "Inverse jacobian not pre-computed");

        const CfemInt n = shape.numNodes;

        // Étape 1 : Calcul de la courbure géométrique (H_geom)
        CfemReal hX_xx = 0, hX_yy = 0, hX_zz = 0, hX_xy = 0, hX_xz = 0, hX_yz = 0;
        CfemReal hY_xx = 0, hY_yy = 0, hY_zz = 0, hY_xy = 0, hY_xz = 0, hY_yz = 0;
        CfemReal hZ_xx = 0, hZ_yy = 0, hZ_zz = 0, hZ_xy = 0, hZ_xz = 0, hZ_yz = 0;

        // Pointeurs SoA pour les Hessiens de référence
        const CfemReal *__restrict d2xi2 = shape.d2N_dxi2();
        const CfemReal *__restrict d2eta2 = shape.d2N_deta2();
        const CfemReal *__restrict d2zeta2 = shape.d2N_dzeta2();
        const CfemReal *__restrict d2xideta = shape.d2N_dxideta();
        const CfemReal *__restrict d2xidzeta = shape.d2N_dxidzeta();
        const CfemReal *__restrict d2etadzeta = shape.d2N_detadzeta();

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

            hX_xx += hxx * pos.x;
            hX_yy += hyy * pos.x;
            hX_zz += hzz * pos.x;
            hX_xy += hxy * pos.x;
            hX_xz += hxz * pos.x;
            hX_yz += hyz * pos.x;

            hY_xx += hxx * pos.y;
            hY_yy += hyy * pos.y;
            hY_zz += hzz * pos.y;
            hY_xy += hxy * pos.y;
            hY_xz += hxz * pos.y;
            hY_yz += hyz * pos.y;

            hZ_xx += hxx * pos.z;
            hZ_yy += hyy * pos.z;
            hZ_zz += hzz * pos.z;
            hZ_xy += hxy * pos.z;
            hZ_xz += hxz * pos.z;
            hZ_yz += hyz * pos.z;
        }

        const CfemReal i00 = out.invJacobian.xx(), i01 = out.invJacobian.xy(), i02 = out.invJacobian.xz();
        const CfemReal i10 = out.invJacobian.yx(), i11 = out.invJacobian.yy(), i12 = out.invJacobian.yz();
        const CfemReal i20 = out.invJacobian.zx(), i21 = out.invJacobian.zy(), i22 = out.invJacobian.zz();

        // Pointeurs SoA pour la lecture des gradients physiques
        const CfemReal *__restrict dNdx = out.dN_dx();
        const CfemReal *__restrict dNdy = out.dN_dy();
        const CfemReal *__restrict dNdz = out.dN_dz();

        // Pointeurs SoA pour l'écriture des Hessiens physiques
        CfemReal *__restrict d2x2 = out.d2N_dx2();
        CfemReal *__restrict d2y2 = out.d2N_dy2();
        CfemReal *__restrict d2z2 = out.d2N_dz2();
        CfemReal *__restrict dxdy = out.d2N_dxdy();
        CfemReal *__restrict dxdz = out.d2N_dxdz();
        CfemReal *__restrict dydz = out.d2N_dydz();

        // --- Étape 2 : Matrice R et Transformation vers l'espace physique ---
        for (CfemInt a = 0; a < n; ++a)
        {
            // =========================================================
            // Construction de R
            // =========================================================

            const CfemReal gx = dNdx[a];
            const CfemReal gy = dNdy[a];
            const CfemReal gz = dNdz[a];

            const CfemReal rxx = d2xi2[a] - (gx * hX_xx + gy * hY_xx + gz * hZ_xx);
            const CfemReal ryy = d2eta2[a] - (gx * hX_yy + gy * hY_yy + gz * hZ_yy);
            const CfemReal rzz = d2zeta2[a] - (gx * hX_zz + gy * hY_zz + gz * hZ_zz);
            const CfemReal rxy = d2xideta[a] - (gx * hX_xy + gy * hY_xy + gz * hZ_xy);
            const CfemReal rxz = d2xidzeta[a] - (gx * hX_xz + gy * hY_xz + gz * hZ_xz);
            const CfemReal ryz = d2etadzeta[a] - (gx * hX_yz + gy * hY_yz + gz * hZ_yz);

            // =========================================================
            // T = R * invJ
            // =========================================================

            const CfemReal t00 = rxx * i00 + rxy * i10 + rxz * i20;
            const CfemReal t01 = rxx * i01 + rxy * i11 + rxz * i21;
            const CfemReal t02 = rxx * i02 + rxy * i12 + rxz * i22;

            const CfemReal t10 = rxy * i00 + ryy * i10 + ryz * i20;
            const CfemReal t11 = rxy * i01 + ryy * i11 + ryz * i21;
            const CfemReal t12 = rxy * i02 + ryy * i12 + ryz * i22;

            const CfemReal t20 = rxz * i00 + ryz * i10 + rzz * i20;
            const CfemReal t21 = rxz * i01 + ryz * i11 + rzz * i21;
            const CfemReal t22 = rxz * i02 + ryz * i12 + rzz * i22;

            // =========================================================
            // H = invJ^T * T
            // =========================================================

            d2x2[a] = i00 * t00 + i10 * t10 + i20 * t20;
            d2y2[a] = i01 * t01 + i11 * t11 + i21 * t21;
            d2z2[a] = i02 * t02 + i12 * t12 + i22 * t22;

            dxdy[a] = i00 * t01 + i10 * t11 + i20 * t21;
            dxdz[a] = i00 * t02 + i10 * t12 + i20 * t22;
            dydz[a] = i01 * t02 + i11 * t12 + i21 * t22;
        }
    }

    inline void ReferenceElement::computePhysicalSecondDerivativesBatch(const DynamicVector<Vector3D> &physicalNodesCoords,
                                                                        const ShapeInfoBatch &shape,
                                                                        MappingInfoBatch &out) const
    {
        CFEM_ASSERT(hasFlag(shape.validFlags, ShapeEvaluationFlags::SecondDerivatives), "Second derivative not computed.");
        CFEM_ASSERT(hasFlag(out.validFlags, MappingEvaluationFlags::InverseJacobian), "Inverse Jacobian not computed.");
        CFEM_ASSERT(hasFlag(out.validFlags, MappingEvaluationFlags::PhysicalFirstDerivatives), "Physical First Derivatives not computed.");

        const CfemInt nQp = shape.numQuadraturePoints;
        const CfemInt n = shape.numNodes;

        for (CfemInt q = 0; q < nQp; ++q)
        {
            // --- Étape 1 : Calcul de la courbure géométrique pour le point q ---
            CfemReal hX_xx = 0, hX_yy = 0, hX_zz = 0, hX_xy = 0, hX_xz = 0, hX_yz = 0;
            CfemReal hY_xx = 0, hY_yy = 0, hY_zz = 0, hY_xy = 0, hY_xz = 0, hY_yz = 0;
            CfemReal hZ_xx = 0, hZ_yy = 0, hZ_zz = 0, hZ_xy = 0, hZ_xz = 0, hZ_yz = 0;

            const CfemReal *__restrict d2xi2 = shape.d2N_dxi2(q);
            const CfemReal *__restrict d2eta2 = shape.d2N_deta2(q);
            const CfemReal *__restrict d2zeta2 = shape.d2N_dzeta2(q);
            const CfemReal *__restrict d2xideta = shape.d2N_dxideta(q);
            const CfemReal *__restrict d2xidzeta = shape.d2N_dxidzeta(q);
            const CfemReal *__restrict d2etadzeta = shape.d2N_detadzeta(q);

            for (CfemInt a = 0; a < n; ++a)
            {
                const Vector3D &pos = physicalNodesCoords[a];
                CfemReal hxx = d2xi2[a];
                CfemReal hyy = d2eta2[a];
                CfemReal hzz = d2zeta2[a];
                CfemReal hxy = d2xideta[a];
                CfemReal hxz = d2xidzeta[a];
                CfemReal hyz = d2etadzeta[a];

                hX_xx += hxx * pos.x;
                hX_yy += hyy * pos.x;
                hX_zz += hzz * pos.x;
                hX_xy += hxy * pos.x;
                hX_xz += hxz * pos.x;
                hX_yz += hyz * pos.x;

                hY_xx += hxx * pos.y;
                hY_yy += hyy * pos.y;
                hY_zz += hzz * pos.y;
                hY_xy += hxy * pos.y;
                hY_xz += hxz * pos.y;
                hY_yz += hyz * pos.y;

                hZ_xx += hxx * pos.z;
                hZ_yy += hyy * pos.z;
                hZ_zz += hzz * pos.z;
                hZ_xy += hxy * pos.z;
                hZ_xz += hxz * pos.z;
                hZ_yz += hyz * pos.z;
            }

            // Récupération de l'inverse du Jacobien pour le point q
            const auto &invJ = out.invJacobians[q];

            // =========================================================
            // PRÉ-CHARGEMENT : Extraction dans les registres locaux
            // =========================================================
            const CfemReal i00 = invJ.xx(), i01 = invJ.xy(), i02 = invJ.xz();
            const CfemReal i10 = invJ.yx(), i11 = invJ.yy(), i12 = invJ.yz();
            const CfemReal i20 = invJ.zx(), i21 = invJ.zy(), i22 = invJ.zz();

            // Recouvrement des pointeurs de lignes physiques (Lecture des gradients, Écriture des Hessiens)
            const CfemReal *__restrict dNdx = out.dN_dx(q);
            const CfemReal *__restrict dNdy = out.dN_dy(q);
            const CfemReal *__restrict dNdz = out.dN_dz(q);

            CfemReal *__restrict d2x2 = out.d2N_dx2(q);
            CfemReal *__restrict d2y2 = out.d2N_dy2(q);
            CfemReal *__restrict d2z2 = out.d2N_dz2(q);
            CfemReal *__restrict dxdy = out.d2N_dxdy(q);
            CfemReal *__restrict dxdz = out.d2N_dxdz(q);
            CfemReal *__restrict dydz = out.d2N_dydz(q);

            // Étape 2 : Calcul de R et application de la congruence pour chaque nœud
            for (CfemInt a = 0; a < n; ++a)
            {
                // =========================================================
                // Construction de R
                // =========================================================
                const CfemReal gx = dNdx[a];
                const CfemReal gy = dNdy[a];
                const CfemReal gz = dNdz[a];

                const CfemReal rxx = d2xi2[a] - (gx * hX_xx + gy * hY_xx + gz * hZ_xx);
                const CfemReal ryy = d2eta2[a] - (gx * hX_yy + gy * hY_yy + gz * hZ_yy);
                const CfemReal rzz = d2zeta2[a] - (gx * hX_zz + gy * hY_zz + gz * hZ_zz);
                const CfemReal rxy = d2xideta[a] - (gx * hX_xy + gy * hY_xy + gz * hZ_xy);
                const CfemReal rxz = d2xidzeta[a] - (gx * hX_xz + gy * hY_xz + gz * hZ_xz);
                const CfemReal ryz = d2etadzeta[a] - (gx * hX_yz + gy * hY_yz + gz * hZ_yz);

                // =========================================================
                // T = R * invJ
                // =========================================================
                const CfemReal t00 = rxx * i00 + rxy * i10 + rxz * i20;
                const CfemReal t01 = rxx * i01 + rxy * i11 + rxz * i21;
                const CfemReal t02 = rxx * i02 + rxy * i12 + rxz * i22;

                const CfemReal t10 = rxy * i00 + ryy * i10 + ryz * i20;
                const CfemReal t11 = rxy * i01 + ryy * i11 + ryz * i21;
                const CfemReal t12 = rxy * i02 + ryy * i12 + ryz * i22;

                const CfemReal t20 = rxz * i00 + ryz * i10 + rzz * i20;
                const CfemReal t21 = rxz * i01 + ryz * i11 + rzz * i21;
                const CfemReal t22 = rxz * i02 + ryz * i12 + rzz * i22;

                // =========================================================
                // H = invJ^T * T
                // =========================================================
                d2x2[a] = i00 * t00 + i10 * t10 + i20 * t20;
                d2y2[a] = i01 * t01 + i11 * t11 + i21 * t21;
                d2z2[a] = i02 * t02 + i12 * t12 + i22 * t22;

                dxdy[a] = i00 * t01 + i10 * t11 + i20 * t21;
                dxdz[a] = i00 * t02 + i10 * t12 + i20 * t22;
                dydz[a] = i01 * t02 + i11 * t12 + i21 * t22;
            }
        }
    }

    inline void ReferenceElement::computeGeometricNormal(MappingInfo &out) const
    {
        CFEM_ASSERT(out.has(MappingEvaluationFlags::Jacobian), "Jacobian not pre-computed");

        if (m_dimension == 2)
        {
            // Variété 2D (Surface, Triangle/Quad) dans l'espace 3D
            // Les deux vecteurs tangents sont les deux premières colonnes du Jacobien
            Vector3D t1(out.jacobian.xx(), out.jacobian.yx(), out.jacobian.zx());
            Vector3D t2(out.jacobian.xy(), out.jacobian.yy(), out.jacobian.zy());

            out.normal = (t1.cross(t2)).normalized();
        }
        else if (m_dimension == 1)
        {
            // Variété 1D (Ligne) dans un espace 2D (ou 3D z=0)
            CfemReal dx = out.jacobian.xx();
            CfemReal dy = out.jacobian.yx();
            CfemReal norm = std::sqrt(dx * dx + dy * dy);

            if (norm > cfem::ZERO_TOLERANCE)
            {
                // Règle de la main droite (rotation de -90 degrés du vecteur tangent)
                // L'orientation extérieure dépend de la convention de rotation de ton maillage
                out.normal = Vector3D(dy / norm, -dx / norm, 0.0);
            }
            else
            {
                out.normal = Vector3D(0.0, 0.0, 0.0);
            }
        }
        else
        {
            // Dimension 3 (Volume). Geometric normal undefined inside a 3D cell.
            out.normal = Vector3D(0.0, 0.0, 0.0);
        }
    }

    inline void ReferenceElement::computeGeometricNormalBatch(MappingInfoBatch &outBatch) const
    {
        CFEM_ASSERT(outBatch.has(MappingEvaluationFlags::Jacobian), "Jacobian not pre-computed");

        const CfemInt nQp = outBatch.numQuadraturePoints;
        const CfemInt qEnd = outBatch.isCellAffine ? 1 : nQp;

        if (m_dimension == 2)
        {
            // Variété 2D (Surface, Triangle/Quad) dans l'espace 3D
            // Les deux vecteurs tangents sont les deux premières colonnes du Jacobien
            for (CfemInt q = 0; q < qEnd; ++q)
            {
                Vector3D t1(outBatch.jacobians[q].xx(), outBatch.jacobians[q].yx(), outBatch.jacobians[q].zx());
                Vector3D t2(outBatch.jacobians[q].xy(), outBatch.jacobians[q].yy(), outBatch.jacobians[q].zy());

                outBatch.normals[q] = (t1.cross(t2)).normalized();
            }
        }
        else if (m_dimension == 1)
        {
            for (CfemInt q = 0; q < qEnd; ++q)
            {
                // Variété 1D (Ligne) dans un espace 2D (ou 3D z=0)
                CfemReal dx = outBatch.jacobians[q].xx();
                CfemReal dy = outBatch.jacobians[q].yx();
                CfemReal norm = std::sqrt(dx * dx + dy * dy);

                if (norm > cfem::ZERO_TOLERANCE)
                {
                    // Règle de la main droite (rotation de -90 degrés du vecteur tangent)
                    // L'orientation extérieure dépend de la convention de rotation de ton maillage
                    outBatch.normals[q] = Vector3D(dy / norm, -dx / norm, 0.0);
                }
                else
                {
                    outBatch.normals[q] = Vector3D(0.0, 0.0, 0.0);
                }
            }
        }
        else
        {
            // Dimension 3 (Volume). Geometric normal undefined inside a 3D cell.
            for (CfemInt q = 0; q < qEnd; ++q)
            {
                outBatch.normals[q] = Vector3D(0.0, 0.0, 0.0);
            }
        }

        if (outBatch.isCellAffine)
        {
            for (CfemInt q = 1; q < nQp; ++q)
            {
                outBatch.normals[q] = outBatch.normals[0];
            }
        }
    }
    // ----------------------------------------------------------------------
    // ------------------------------- interpolation ------------------------
    template <typename T>
    inline T ReferenceElement::interpolate(const ShapeInfo &shape,
                                           const DynamicVector<T> &nodalValues) const
    {
        const CfemInt n = shape.numNodes;
        CFEM_ASSERT(shape.has(ShapeEvaluationFlags::Value), "Shape values not computed");
        CFEM_ASSERT(nodalValues.size() == n, "Nodal values size mismatch");

        T result{};
        const CfemReal *N = shape.values();

        for (CfemInt i = 0; i < n; ++i)
        {
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
        const CfemReal *dx = map.dN_dx();
        const CfemReal *dy = map.dN_dy();
        const CfemReal *dz = map.dN_dz();

        for (CfemInt i = 0; i < n; ++i)
        {
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

        const CfemReal *d2x = map.d2N_dx2();
        const CfemReal *d2y = map.d2N_dy2();
        const CfemReal *d2z = map.d2N_dz2();
        const CfemReal *dxdy = map.d2N_dxdy();
        const CfemReal *dxdz = map.d2N_dxdz();
        const CfemReal *dydz = map.d2N_dydz();

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

} // namespace cfem