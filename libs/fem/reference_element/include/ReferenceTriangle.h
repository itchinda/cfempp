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

#ifndef CFEM_REFERENCE_TRIANGLE
#define CFEM_REFERENCE_TRIANGLE

#include "ReferenceElement.h"
#include "Vector3D.h"

#include <array>
#include <span>
#include <cstring>

namespace cfem::reference_elem::utils{

    inline Vector3D projectFacetRefCoordToTriRefCoord(CfemInt localFacetId, 
                                                      const Vector3D& xiFacet) 
    {
        const CfemReal u = xiFacet.x; // [-1, 1] pour une ligne
        CfemReal n1 = 0.5 * (1.0 - u);
        CfemReal n2 = 0.5 * (1.0 + u);

            switch (localFacetId) {
                case 0: return Vector3D(n2, 0.0, 0.0); // Face {0,1}
                case 1: return Vector3D(n1, n2, 0.0); // Face {1,2}
                case 2: return Vector3D(0.0, n1, 0.0); // Face {2,0}
                default: CFEM_ERROR << "Invalid facet ID for Triangle";
            }
        return Vector3D{}; // Make the compiler happy
    }
   
    static constexpr std::array<std::array<CfemInt, 1>, 3> tri1_face_nodes = {{
        {0}, {0}, {0}
    }};

    static constexpr std::array<std::array<CfemInt, 2>, 3> tri3_face_nodes = {{
        {0, 1}, {1, 2}, {2, 0}
    }};
    
    static constexpr std::array<std::array<CfemInt, 3>, 3> tri6_face_nodes = {{
        {0, 1, 3}, {1, 2, 4}, {2, 0, 5}
    }};

    inline void computeTriBarycentric(const Vector3D& xi,
                                      CfemReal& L0,
                                      CfemReal& L1,
                                      CfemReal& L2) noexcept
    {
        L1 = xi.x;
        L2 = xi.y;
        L0 = 1.0 - L1 - L2;
    }

}

namespace cfem::reference_elem::kernels {

    // =========================================================================
    // P0 TRIANGLE KERNEL (1 node, constant)
    // =========================================================================
    struct P0TriKernel {
        [[gnu::always_inline]] static inline void values(CfemReal* __restrict N) noexcept {
            N[0] = 1.0;
        }

        [[gnu::always_inline]] static inline void grads(CfemReal* __restrict dx, 
                                 CfemReal* __restrict dy, 
                                 CfemReal* __restrict dz) noexcept {
            dx[0] = dy[0] = dz[0] = 0.0;
        }

        [[gnu::always_inline]] static inline void hess(CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept {
            d2x[0] = d2y[0] = d2z[0] = dxy[0] = dxz[0] = dyz[0] = 0.0;
        }
    };

    // =========================================================================
    // P1 TRIANGLE KERNEL (3 nodes, linear)
    // =========================================================================
    struct P1TriKernel {
        [[gnu::always_inline]] static inline void values(CfemReal L0, CfemReal L1, CfemReal L2, 
                                  CfemReal* __restrict N) noexcept {
            N[0] = L0; N[1] = L1; N[2] = L2;
        }

        [[gnu::always_inline]] static inline void grads(CfemReal* __restrict dx, 
                                 CfemReal* __restrict dy, 
                                 CfemReal* __restrict dz) noexcept {
            dx[0] = -1.0; dy[0] = -1.0; dz[0] = 0.0;
            dx[1] =  1.0; dy[1] =  0.0; dz[1] = 0.0;
            dx[2] =  0.0; dy[2] =  1.0; dz[2] = 0.0;
        }

        [[gnu::always_inline]] static inline void hess(CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept {
            for (CfemInt i = 0; i < 3; ++i) {
                d2x[i] = d2y[i] = d2z[i] = dxy[i] = dxz[i] = dyz[i] = 0.0;
            }
        }
    };

    // =========================================================================
    // P1-MINI KERNEL TRIANGLE (4 nodes : 3 vertices + 1 bubble)
    // =========================================================================
    struct MiniTriKernel {
        [[gnu::always_inline]] static inline void values(CfemReal L0, CfemReal L1, CfemReal L2, 
                                  CfemReal* __restrict N) noexcept {
            const CfemReal b = 27.0 * L0 * L1 * L2;
            N[3] = b;
            N[0] = L0 - b / 3.0;
            N[1] = L1 - b / 3.0;
            N[2] = L2 - b / 3.0;
        }

        [[gnu::always_inline]] static inline void grads(CfemReal L0, CfemReal L1, CfemReal L2, 
                                 CfemReal* __restrict dx, 
                                 CfemReal* __restrict dy, 
                                 CfemReal* __restrict dz) noexcept {
            const CfemReal db_dxi  = 27.0 * L2 * (L0 - L1);
            const CfemReal db_deta = 27.0 * L1 * (L0 - L2);

            dx[3] = db_dxi;              dy[3] = db_deta;             dz[3] = 0.0;
            dx[0] = -1.0 - db_dxi / 3.0; dy[0] = -1.0 - db_deta / 3.0; dz[0] = 0.0;
            dx[1] =  1.0 - db_dxi / 3.0; dy[1] =  0.0 - db_deta / 3.0; dz[1] = 0.0;
            dx[2] =  0.0 - db_dxi / 3.0; dy[2] =  1.0 - db_deta / 3.0; dz[2] = 0.0;
        }

        [[gnu::always_inline]] static inline void hess(CfemReal L0, CfemReal L1, CfemReal L2, 
                                CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept {
            const CfemReal d2b_dxi2    = -54.0 * L2;
            const CfemReal d2b_deta2   = -54.0 * L1;
            const CfemReal d2b_dxideta =  27.0 * (L0 - L1 - L2);

            d2x[3] = d2b_dxi2;    d2y[3] = d2b_deta2;    dxy[3] = d2b_dxideta;
            d2z[3] = dxz[3] = dyz[3] = 0.0;

            for (CfemInt i = 0; i < 3; ++i) {
                d2x[i] = -d2b_dxi2 / 3.0;
                d2y[i] = -d2b_deta2 / 3.0;
                dxy[i] = -d2b_dxideta / 3.0;
                d2z[i] = dxz[i] = dyz[i] = 0.0;
            }
        }
    };

    // =========================================================================
    // P2 TRIANGLE KERNEL (6 nodes, quadratic)
    // =========================================================================
    struct P2TriKernel {
        [[gnu::always_inline]] static inline void values(CfemReal L0, CfemReal L1, CfemReal L2, 
                                  CfemReal* __restrict N) noexcept {
            N[0] = L0 * (2.0 * L0 - 1.0);
            N[1] = L1 * (2.0 * L1 - 1.0);
            N[2] = L2 * (2.0 * L2 - 1.0);
            N[3] = 4.0 * L0 * L1;
            N[4] = 4.0 * L1 * L2;
            N[5] = 4.0 * L2 * L0;
        }

        [[gnu::always_inline]] static inline void grads(CfemReal L0, CfemReal L1, CfemReal L2, 
                                 CfemReal* __restrict dx, 
                                 CfemReal* __restrict dy, 
                                 CfemReal* __restrict dz) noexcept {
            dx[0] = 1.0 - 4.0 * L0;      dy[0] = 1.0 - 4.0 * L0;
            dx[1] = 4.0 * L1 - 1.0;      dy[1] = 0.0;
            dx[2] = 0.0;                 dy[2] = 4.0 * L2 - 1.0;
            dx[3] = 4.0 * (L0 - L1);     dy[3] = -4.0 * L1;
            dx[4] = 4.0 * L2;            dy[4] = 4.0 * L1;
            dx[5] = -4.0 * L2;           dy[5] = 4.0 * (L0 - L2);
            
            for (CfemInt i=0; i<6; ++i) dz[i] = 0.0;
        }

        [[gnu::always_inline]] static inline void hess(CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept {
            d2x[0] =  4.0; d2y[0] =  4.0; dxy[0] =  4.0;
            d2x[1] =  4.0; d2y[1] =  0.0; dxy[1] =  0.0;
            d2x[2] =  0.0; d2y[2] =  4.0; dxy[2] =  0.0;
            d2x[3] = -8.0; d2y[3] =  0.0; dxy[3] = -4.0;
            d2x[4] =  0.0; d2y[4] =  0.0; dxy[4] =  4.0;
            d2x[5] =  0.0; d2y[5] = -8.0; dxy[5] = -4.0;

            for (CfemInt i=0; i<6; ++i) {
                d2z[i] = dxz[i] = dyz[i] = 0.0;
            }
        }
    };
}

namespace cfem {

    /**
     * @brief Constant Reference Triangle (1 node, P0)
     */
    class P0ReferenceTriangle : public ReferenceElement {
    public:
        P0ReferenceTriangle() {
            m_numNodes = 1;
            m_dimension = 2;
            m_nodes = { Vector3D{1.0/3.0, 1.0/3.0, 0.0} };
        }

        void evaluateShapesValues(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::P0TriKernel::values(info.values());
        }

        void evaluateShapesDerivatives(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::P0TriKernel::grads(info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::P0TriKernel::hess(info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(), 
                                                       info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            if (nQp == 0) return;
            reference_elem::kernels::P0TriKernel::values(shapeBatch.values(0));
            const size_t bytesToCopy = m_numNodes * sizeof(CfemReal);
            for (CfemInt q = 1; q < nQp; ++q) {
                std::memcpy(shapeBatch.values(q), shapeBatch.values(0), bytesToCopy);
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> & xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            if (nQp == 0) return;
            reference_elem::kernels::P0TriKernel::grads(shapeBatch.dN_dxi(0), shapeBatch.dN_deta(0), shapeBatch.dN_dzeta(0));
            const size_t bytesToCopy = m_numNodes * sizeof(CfemReal);
            for (CfemInt q = 1; q < nQp; ++q) {
                std::memcpy(shapeBatch.dN_dxi(q),   shapeBatch.dN_dxi(0),   bytesToCopy);
                std::memcpy(shapeBatch.dN_deta(q),  shapeBatch.dN_deta(0),  bytesToCopy);
                std::memcpy(shapeBatch.dN_dzeta(q), shapeBatch.dN_dzeta(0), bytesToCopy);
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            if (nQp == 0) return;
            reference_elem::kernels::P0TriKernel::hess(shapeBatch.d2N_dxi2(0), shapeBatch.d2N_deta2(0), shapeBatch.d2N_dzeta2(0),
                                                       shapeBatch.d2N_dxideta(0), shapeBatch.d2N_dxidzeta(0), shapeBatch.d2N_detadzeta(0));
            const size_t bytesToCopy = m_numNodes * sizeof(CfemReal);
            for (CfemInt q = 1; q < nQp; ++q) {
                std::memcpy(shapeBatch.d2N_dxi2(q),      shapeBatch.d2N_dxi2(0),      bytesToCopy);
                std::memcpy(shapeBatch.d2N_deta2(q),     shapeBatch.d2N_deta2(0),     bytesToCopy);
                std::memcpy(shapeBatch.d2N_dzeta2(q),    shapeBatch.d2N_dzeta2(0),    bytesToCopy);
                std::memcpy(shapeBatch.d2N_dxideta(q),   shapeBatch.d2N_dxideta(0),   bytesToCopy);
                std::memcpy(shapeBatch.d2N_dxidzeta(q),  shapeBatch.d2N_dxidzeta(0),  bytesToCopy);
                std::memcpy(shapeBatch.d2N_detadzeta(q), shapeBatch.d2N_detadzeta(0), bytesToCopy);
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet);
        }
        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tri1_face_nodes[iFace]; 
        }
    };

    /**
     * @brief Linear Reference Triangle (3 nodes, P1)
     */
    class P1ReferenceTriangle : public ReferenceElement {
    public:
        P1ReferenceTriangle() {
            m_numNodes = 3;
            m_dimension = 2;
            m_nodes = { Vector3D{0.0, 0.0, 0.0}, Vector3D{1.0, 0.0, 0.0}, Vector3D{0.0, 1.0, 0.0} };
        }

        void evaluateShapesValues(const Vector3D& xi, ShapeInfo& info) const override {
            CfemReal L0, L1, L2;
            reference_elem::utils::computeTriBarycentric(xi, L0, L1, L2);
            reference_elem::kernels::P1TriKernel::values(L0, L1, L2, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::P1TriKernel::grads(info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::P1TriKernel::hess(info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(), 
                                                       info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                CfemReal L0, L1, L2;
                reference_elem::utils::computeTriBarycentric(xiPoints[q], L0, L1, L2);
                reference_elem::kernels::P1TriKernel::values(L0, L1, L2, shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P1TriKernel::grads(shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            if (nQp == 0) return;
            reference_elem::kernels::P1TriKernel::hess(shapeBatch.d2N_dxi2(0), shapeBatch.d2N_deta2(0), shapeBatch.d2N_dzeta2(0),
                                                       shapeBatch.d2N_dxideta(0), shapeBatch.d2N_dxidzeta(0), shapeBatch.d2N_detadzeta(0));
            const size_t bytesToCopy = m_numNodes * sizeof(CfemReal);
            for (CfemInt q = 1; q < nQp; ++q) {
                std::memcpy(shapeBatch.d2N_dxi2(q),      shapeBatch.d2N_dxi2(0),      bytesToCopy);
                std::memcpy(shapeBatch.d2N_deta2(q),     shapeBatch.d2N_deta2(0),     bytesToCopy);
                std::memcpy(shapeBatch.d2N_dzeta2(q),    shapeBatch.d2N_dzeta2(0),    bytesToCopy);
                std::memcpy(shapeBatch.d2N_dxideta(q),   shapeBatch.d2N_dxideta(0),   bytesToCopy);
                std::memcpy(shapeBatch.d2N_dxidzeta(q),  shapeBatch.d2N_dxidzeta(0),  bytesToCopy);
                std::memcpy(shapeBatch.d2N_detadzeta(q), shapeBatch.d2N_detadzeta(0), bytesToCopy);
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet); 
        }
        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tri3_face_nodes[iFace];
        }
    };

    /**
     * @brief P1-Mini Reference Triangle (4 nodes: 3 vertices + 1 bubble)
     */
    class MiniReferenceTriangle : public ReferenceElement {
    public:
        MiniReferenceTriangle() {
            m_numNodes = 4;
            m_dimension = 2;
            m_nodes = { 
                Vector3D{0.0, 0.0, 0.0}, Vector3D{1.0, 0.0, 0.0}, Vector3D{0.0, 1.0, 0.0}, 
                Vector3D{1.0/3.0, 1.0/3.0, 0.0} 
            };
        }

        void evaluateShapesValues(const Vector3D& xi, ShapeInfo& info) const override {
            CfemReal L0, L1, L2;
            reference_elem::utils::computeTriBarycentric(xi, L0, L1, L2);
            reference_elem::kernels::MiniTriKernel::values(L0, L1, L2, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D& xi, ShapeInfo& info) const override {
            CfemReal L0, L1, L2;
            reference_elem::utils::computeTriBarycentric(xi, L0, L1, L2);
            reference_elem::kernels::MiniTriKernel::grads(L0, L1, L2, info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D& xi, ShapeInfo& info) const override {
            CfemReal L0, L1, L2;
            reference_elem::utils::computeTriBarycentric(xi, L0, L1, L2);
            reference_elem::kernels::MiniTriKernel::hess(L0, L1, L2, info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(), 
                                                         info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                CfemReal L0, L1, L2;
                reference_elem::utils::computeTriBarycentric(xiPoints[q], L0, L1, L2);
                reference_elem::kernels::MiniTriKernel::values(L0, L1, L2, shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                CfemReal L0, L1, L2;
                reference_elem::utils::computeTriBarycentric(xiPoints[q], L0, L1, L2);
                reference_elem::kernels::MiniTriKernel::grads(L0, L1, L2, shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            // Dépendance de \xi pour le hessien de la bulle (pas de memcpy)
            for (CfemInt q = 0; q < nQp; ++q) {
                CfemReal L0, L1, L2;
                reference_elem::utils::computeTriBarycentric(xiPoints[q], L0, L1, L2);
                reference_elem::kernels::MiniTriKernel::hess(L0, L1, L2, shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                             shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet); 
        }
        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tri3_face_nodes[iFace]; 
        }
    };

    /**
     * @brief Quadratic Reference Triangle (6 nodes, P2)
     */
    class P2ReferenceTriangle : public ReferenceElement {
    public:
        P2ReferenceTriangle() {
            m_numNodes = 6;
            m_dimension = 2;
            m_nodes = { 
                Vector3D{0.0, 0.0, 0.0}, Vector3D{1.0, 0.0, 0.0}, Vector3D{0.0, 1.0, 0.0},
                Vector3D{0.5, 0.0, 0.0}, Vector3D{0.5, 0.5, 0.0}, Vector3D{0.0, 0.5, 0.0}
            };
        }

        void evaluateShapesValues(const Vector3D& xi, ShapeInfo& info) const override {
            CfemReal L0, L1, L2;
            reference_elem::utils::computeTriBarycentric(xi, L0, L1, L2);
            reference_elem::kernels::P2TriKernel::values(L0, L1, L2, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D& xi, ShapeInfo& info) const override {
            CfemReal L0, L1, L2;
            reference_elem::utils::computeTriBarycentric(xi, L0, L1, L2);
            reference_elem::kernels::P2TriKernel::grads(L0, L1, L2, info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::P2TriKernel::hess(info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(), 
                                                       info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                CfemReal L0, L1, L2;
                reference_elem::utils::computeTriBarycentric(xiPoints[q], L0, L1, L2);
                reference_elem::kernels::P2TriKernel::values(L0, L1, L2, shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                CfemReal L0, L1, L2;
                reference_elem::utils::computeTriBarycentric(xiPoints[q], L0, L1, L2);
                reference_elem::kernels::P2TriKernel::grads(L0, L1, L2, shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            if (nQp == 0) return;

            // Hessien constant (indépendant de xi), on utilise le broadcast Memcpy
            reference_elem::kernels::P2TriKernel::hess(shapeBatch.d2N_dxi2(0), shapeBatch.d2N_deta2(0), shapeBatch.d2N_dzeta2(0),
                                                       shapeBatch.d2N_dxideta(0), shapeBatch.d2N_dxidzeta(0), shapeBatch.d2N_detadzeta(0));

            const size_t bytesToCopy = m_numNodes * sizeof(CfemReal);
            for (CfemInt q = 1; q < nQp; ++q) {
                std::memcpy(shapeBatch.d2N_dxi2(q),      shapeBatch.d2N_dxi2(0),      bytesToCopy);
                std::memcpy(shapeBatch.d2N_deta2(q),     shapeBatch.d2N_deta2(0),     bytesToCopy);
                std::memcpy(shapeBatch.d2N_dzeta2(q),    shapeBatch.d2N_dzeta2(0),    bytesToCopy);
                std::memcpy(shapeBatch.d2N_dxideta(q),   shapeBatch.d2N_dxideta(0),   bytesToCopy);
                std::memcpy(shapeBatch.d2N_dxidzeta(q),  shapeBatch.d2N_dxidzeta(0),  bytesToCopy);
                std::memcpy(shapeBatch.d2N_detadzeta(q), shapeBatch.d2N_detadzeta(0), bytesToCopy);
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet); 
        }
        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tri6_face_nodes[iFace];
        }
    };
}


#endif