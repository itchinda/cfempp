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

#ifndef CFEM_REFERENCE_SEGMENT_H
#define CFEM_REFERENCE_SEGMENT_H

#include "ReferenceElement.h"
#include <array>
#include <span>

namespace cfem::reference_elem::utils {

    inline Vector3D projectFacetRefCoordToSegmentRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) {
        return Vector3D((localFacetId == 0) ? -1.0 : 1.0, 0.0, 0.0);
    }

    static const std::array<std::array<CfemInt, 1>, 2> line_face_nodes = {{
        {0}, // Face 0
        {1}  // Face 1
    }};

    static const std::array<std::array<CfemInt, 1>, 2> line1_face_nodes = {{
        {0}, // Face 0 (x = -1) -> Node 0
        {0}  // Face 1 (x =  1) -> Node 1
    }};
}

namespace cfem::reference_elem::kernels {

    // =========================================================================
    // P1 SEGMENT KERNEL (2 nodes, linear)
    // =========================================================================
    struct P1SegmentKernel {
        [[gnu::always_inline]] static inline void values(CfemReal x, CfemReal* __restrict N) noexcept {
            N[0] = 0.5 * (1.0 - x);
            N[1] = 0.5 * (1.0 + x);
        }

        [[gnu::always_inline]] static inline void grads(CfemReal* __restrict dx, 
                                                        CfemReal* __restrict dy, 
                                                        CfemReal* __restrict dz) noexcept {
            dx[0] = -0.5;
            dx[1] =  0.5;
            dy[0] = dy[1] = 0.0;
            dz[0] = dz[1] = 0.0;
        }

        [[gnu::always_inline]] static inline void hess(CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                                       CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept {
            for (CfemInt i = 0; i < 2; ++i) {
                d2x[i] = 0.0; d2y[i] = 0.0; d2z[i] = 0.0;
                dxy[i] = 0.0; dxz[i] = 0.0; dyz[i] = 0.0;
            }
        }
    };

    // =========================================================================
    // P2 SEGMENT KERNEL (3 nodes, quadratic)
    // =========================================================================
    struct P2SegmentKernel {
        [[gnu::always_inline]] static inline void values(CfemReal x, CfemReal* __restrict N) noexcept {
            const CfemReal x2 = x * x;
            N[0] = 0.5 * (x2 - x);
            N[1] = 0.5 * (x2 + x);
            N[2] = 1.0 - x2;
        }

        [[gnu::always_inline]] static inline void grads(CfemReal x, 
                                                        CfemReal* __restrict dx, 
                                                        CfemReal* __restrict dy, 
                                                        CfemReal* __restrict dz) noexcept {
            dx[0] = x - 0.5;
            dx[1] = x + 0.5;
            dx[2] = -2.0 * x;

            for (CfemInt i = 0; i < 3; ++i) {
                dy[i] = 0.0;
                dz[i] = 0.0;
            }
        }

        [[gnu::always_inline]] static inline void hess(CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                                       CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept {
            d2x[0] =  1.0;
            d2x[1] =  1.0;
            d2x[2] = -2.0;

            for (CfemInt i = 0; i < 3; ++i) {
                d2y[i] = 0.0; d2z[i] = 0.0;
                dxy[i] = 0.0; dxz[i] = 0.0; dyz[i] = 0.0;
            }
        }
    };

    // =========================================================================
    // KERNEL SEGMENT P0 (1 nœud, constant)
    // =========================================================================
    struct P0SegmentKernel {
        [[gnu::always_inline]] static inline void values(CfemReal* __restrict N) noexcept {
            N[0] = 1.0;
        }

        [[gnu::always_inline]] static inline void grads(CfemReal* __restrict dx, 
                                                        CfemReal* __restrict dy, 
                                                        CfemReal* __restrict dz) noexcept {
            dx[0] = 0.0; dy[0] = 0.0; dz[0] = 0.0;
        }

        [[gnu::always_inline]] static inline void hess(CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                                       CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept {
            d2x[0] = 0.0; d2y[0] = 0.0; d2z[0] = 0.0;
            dxy[0] = 0.0; dxz[0] = 0.0; dyz[0] = 0.0;
        }
    };
}

namespace cfem {

    /**
     * @brief P1 Reference Segment (2 nodes)
     */
    class P1ReferenceSegment : public ReferenceElement {
    public:
        P1ReferenceSegment() {
            m_numNodes = 2;
            m_dimension = 1;
            m_nodes.resize(m_numNodes);
            m_nodes[0] = Vector3D{-1.0, 0.0, 0.0};
            m_nodes[1] = Vector3D{ 1.0, 0.0, 0.0};
        }

        void evaluateShapesValues(const Vector3D &xi, ShapeInfo& info) const override {
            reference_elem::kernels::P1SegmentKernel::values(xi.x, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D &xi, ShapeInfo& info) const override {
            reference_elem::kernels::P1SegmentKernel::grads(info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xi, ShapeInfo& info) const override {
            reference_elem::kernels::P1SegmentKernel::hess(info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(),
                                                           info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P1SegmentKernel::values(xiPoints[q].x, shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P1SegmentKernel::grads(shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P1SegmentKernel::hess(shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                               shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToSegmentRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::line_face_nodes[iFace];
        }
    };


    /**
     * @brief P2 Reference Segment (3 nodes)
     */
    class P2ReferenceSegment : public ReferenceElement {
    public:
        P2ReferenceSegment() {
            m_numNodes = 3;
            m_dimension = 1;
            m_nodes.resize(m_numNodes);
            // Domaine de référence [-1, 1], ordre VTK : extrémités puis milieu
            m_nodes[0] = Vector3D{-1.0, 0.0, 0.0};
            m_nodes[1] = Vector3D{ 1.0, 0.0, 0.0};
            m_nodes[2] = Vector3D{ 0.0, 0.0, 0.0}; // Nœud central en x=0
        }

        void evaluateShapesValues(const Vector3D &xi, ShapeInfo& info) const override {
            reference_elem::kernels::P2SegmentKernel::values(xi.x, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D &xi, ShapeInfo& info) const override {
            reference_elem::kernels::P2SegmentKernel::grads(xi.x, info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xi, ShapeInfo& info) const override {
            reference_elem::kernels::P2SegmentKernel::hess(info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(),
                                                           info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P2SegmentKernel::values(xiPoints[q].x, shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P2SegmentKernel::grads(xiPoints[q].x, shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P2SegmentKernel::hess(shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                               shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToSegmentRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::line_face_nodes[iFace];
        }
    };


    /**
     * @brief P0 Reference Segment (1 node, constant)
     */
    class P0ReferenceSegment : public ReferenceElement {
    public:
        P0ReferenceSegment() {
            m_numNodes = 1;
            m_dimension = 1;
            m_nodes.resize(m_numNodes);
            m_nodes[0] = Vector3D{0.0, 0.0, 0.0};
        }

        void evaluateShapesValues(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::P0SegmentKernel::values(info.values());
        }

        void evaluateShapesDerivatives(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::P0SegmentKernel::grads(info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::P0SegmentKernel::hess(info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(),
                                                           info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P0SegmentKernel::values(shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P0SegmentKernel::grads(shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P0SegmentKernel::hess(shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                               shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToSegmentRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::line1_face_nodes[iFace];
        }
    };
}

#endif