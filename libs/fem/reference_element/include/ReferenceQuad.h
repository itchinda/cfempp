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

#ifndef CFEM_REFERENCE_QUAD
#define CFEM_REFERENCE_QUAD

#include "ReferenceElement.h"
#include <array>
#include <span>

namespace cfem::reference_elem::utils {

    inline Vector3D projectFacetRefCoordToQuadRefCoord(CfemInt localFacetId, const Vector3D &xiFacet)
    {
        const CfemReal u = xiFacet.x;

        switch (localFacetId)
        {
            case 0: return Vector3D(u, -1.0, 0.0);  // Face {0,1} : y=-1
            case 1: return Vector3D(1.0, u, 0.0);   // Face {1,2} : x=1
            case 2: return Vector3D(-u, 1.0, 0.0);  // Face {2,3} : y=1
            case 3: return Vector3D(-1.0, -u, 0.0); // Face {3,0} : x=-1
            default: CFEM_ERROR << "Invalid facet ID for Quad cell";
        }
        return Vector3D{}; // Make the compiler happy
    }

    static const std::array<std::array<CfemInt, 2>, 4> quad4_face_nodes = {{
        {0, 1}, {1, 2}, {2, 3}, {3, 0}
    }};

    static const std::array<std::array<CfemInt, 3>, 4> quad9_face_nodes = {{
        {0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7}
    }};

    static const std::array<std::array<CfemInt, 1>, 4> quad1_face_nodes = {{
        {0}, {0}, {0}, {0}
    }};
}

namespace cfem::reference_elem::kernels {

    // =========================================================================
    // KERNEL QUAD Q1 (4 nœuds, bilinéaire)
    // =========================================================================
    struct Q1QuadKernel {
        static constexpr CfemReal o4 = 0.25;

        [[gnu::always_inline]] static inline void values(CfemReal x, CfemReal y, CfemReal z, 
                                                         CfemReal* __restrict N) noexcept 
        {
            const CfemReal xm = 1.0 - x;
            const CfemReal xp = 1.0 + x;
            const CfemReal ym = 1.0 - y;
            const CfemReal yp = 1.0 + y;

            N[0] = o4 * xm * ym;
            N[1] = o4 * xp * ym;
            N[2] = o4 * xp * yp;
            N[3] = o4 * xm * yp;
        }

        [[gnu::always_inline]] static inline void grads(CfemReal x, CfemReal y, CfemReal z, 
                                                        CfemReal* __restrict dx, 
                                                        CfemReal* __restrict dy, 
                                                        CfemReal* __restrict dz) noexcept 
        {
            const CfemReal xm = 1.0 - x;
            const CfemReal xp = 1.0 + x;
            const CfemReal ym = 1.0 - y;
            const CfemReal yp = 1.0 + y;

            dx[0] = -o4 * ym;
            dx[1] =  o4 * ym;
            dx[2] =  o4 * yp;
            dx[3] = -o4 * yp;

            dy[0] = -o4 * xm;
            dy[1] = -o4 * xp;
            dy[2] =  o4 * xp;
            dy[3] =  o4 * xm;

            dz[0] = dz[1] = dz[2] = dz[3] = 0.0;
        }

        [[gnu::always_inline]] static inline void hess(CfemReal x, CfemReal y, CfemReal z, 
                                                       CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                                       CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept 
        {
            for (CfemInt i = 0; i < 4; ++i) {
                d2x[i] = 0.0;
                d2y[i] = 0.0;
                d2z[i] = 0.0;
                dxz[i] = 0.0;
                dyz[i] = 0.0;
            }

            dxy[0] =  o4;
            dxy[1] = -o4;
            dxy[2] =  o4;
            dxy[3] = -o4;
        }
    };

    // =========================================================================
    // KERNEL QUAD Q2 (9 nœuds, biquadratique)
    // =========================================================================
    struct Q2QuadKernel {
        static constexpr CfemInt node_map[3][3] = {
            {0, 7, 3}, // X = -1 | Y = -1, 0, 1
            {4, 8, 6}, // X =  0 | Y = -1, 0, 1
            {1, 5, 2}  // X = +1 | Y = -1, 0, 1
        };

        [[gnu::always_inline]] static inline void eval1D_hessian(CfemReal u, CfemReal* N, CfemReal* dN, CfemReal* d2N) noexcept
        {
            const CfemReal u2 = u * u;
            N[0] = 0.5 * (u2 - u); dN[0] = u - 0.5;   d2N[0] = 1.0;
            N[1] = 1.0 - u2;       dN[1] = -2.0 * u;  d2N[1] = -2.0;
            N[2] = 0.5 * (u2 + u); dN[2] = u + 0.5;   d2N[2] = 1.0;
        }

        [[gnu::always_inline]] static inline void values(CfemReal x, CfemReal y, CfemReal z, 
                                                         CfemReal* __restrict N) noexcept 
        {
            CfemReal Nx[3], Ny[3], dNx[3], dNy[3], d2Nx[3], d2Ny[3];
            eval1D_hessian(x, Nx, dNx, d2Nx);
            eval1D_hessian(y, Ny, dNy, d2Ny);

            for (CfemInt i = 0; i < 3; ++i) {
                for (CfemInt j = 0; j < 3; ++j) {
                    N[node_map[i][j]] = Nx[i] * Ny[j];
                }
            }
        }

        [[gnu::always_inline]] static inline void grads(CfemReal x, CfemReal y, CfemReal z, 
                                                        CfemReal* __restrict dx, 
                                                        CfemReal* __restrict dy, 
                                                        CfemReal* __restrict dz) noexcept 
        {
            CfemReal Nx[3], Ny[3], dNx[3], dNy[3], d2Nx[3], d2Ny[3];
            eval1D_hessian(x, Nx, dNx, d2Nx);
            eval1D_hessian(y, Ny, dNy, d2Ny);

            for (CfemInt i = 0; i < 3; ++i) {
                for (CfemInt j = 0; j < 3; ++j) {
                    const int node = node_map[i][j];
                    dx[node] = dNx[i] * Ny[j];
                    dy[node] = Nx[i] * dNy[j];
                    dz[node] = 0.0;
                }
            }
        }

        [[gnu::always_inline]] static inline void hess(CfemReal x, CfemReal y, CfemReal z, 
                                                       CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                                       CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept 
        {
            CfemReal Nx[3], Ny[3], dNx[3], dNy[3], d2Nx[3], d2Ny[3];
            eval1D_hessian(x, Nx, dNx, d2Nx);
            eval1D_hessian(y, Ny, dNy, d2Ny);

            for (CfemInt i = 0; i < 3; ++i) {
                for (CfemInt j = 0; j < 3; ++j) {
                    const int node = node_map[i][j];
                    d2x[node] = d2Nx[i] * Ny[j];
                    d2y[node] = Nx[i] * d2Ny[j];
                    d2z[node] = 0.0;

                    dxy[node] = dNx[i] * dNy[j];
                    dxz[node] = 0.0;
                    dyz[node] = 0.0;
                }
            }
        }
    };

    // =========================================================================
    // KERNEL QUAD MINI (5 nœuds : 4 sommets + 1 bulle)
    // =========================================================================
    struct MiniQuadKernel {
        static constexpr CfemReal o4 = 0.25;

        [[gnu::always_inline]] static inline void values(CfemReal x, CfemReal y, CfemReal z, 
                                                         CfemReal* __restrict N) noexcept 
        {
            // Bulle Centrale
            const CfemReal b = (1.0 - x*x) * (1.0 - y*y);
            N[4] = b;

            // Sommets
            const CfemReal correction = b * o4;
            N[0] = (o4 * (1.0 - x) * (1.0 - y)) - correction;
            N[1] = (o4 * (1.0 + x) * (1.0 - y)) - correction;
            N[2] = (o4 * (1.0 + x) * (1.0 + y)) - correction;
            N[3] = (o4 * (1.0 - x) * (1.0 + y)) - correction;
        }

        [[gnu::always_inline]] static inline void grads(CfemReal x, CfemReal y, CfemReal z, 
                                                        CfemReal* __restrict dx, 
                                                        CfemReal* __restrict dy, 
                                                        CfemReal* __restrict dz) noexcept 
        {
            const CfemReal db_dx = -2.0 * x * (1.0 - y*y);
            const CfemReal db_dy = -2.0 * y * (1.0 - x*x);

            dx[4] = db_dx;
            dy[4] = db_dy;
            dz[4] = 0.0;

            const CfemReal c_dx = db_dx * o4;
            const CfemReal c_dy = db_dy * o4;

            dx[0] = (-o4 * (1.0 - y)) - c_dx; dy[0] = (-o4 * (1.0 - x)) - c_dy; dz[0] = 0.0;
            dx[1] = ( o4 * (1.0 - y)) - c_dx; dy[1] = (-o4 * (1.0 + x)) - c_dy; dz[1] = 0.0;
            dx[2] = ( o4 * (1.0 + y)) - c_dx; dy[2] = ( o4 * (1.0 + x)) - c_dy; dz[2] = 0.0;
            dx[3] = (-o4 * (1.0 + y)) - c_dx; dy[3] = ( o4 * (1.0 - x)) - c_dy; dz[3] = 0.0;
        }

        [[gnu::always_inline]] static inline void hess(CfemReal x, CfemReal y, CfemReal z, 
                                                       CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                                       CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept 
        {
            // Hessiens de la bulle
            const CfemReal d2b_dx2 = -2.0 * (1.0 - y*y);
            const CfemReal d2b_dy2 = -2.0 * (1.0 - x*x);
            const CfemReal d2b_dxy =  4.0 * x * y;

            d2x[4] = d2b_dx2;
            d2y[4] = d2b_dy2;
            dxy[4] = d2b_dxy;
            d2z[4] = 0.0; dxz[4] = 0.0; dyz[4] = 0.0;

            // Sommets (Correction de la composante purement nulle pour xx et yy du Q1)
            const CfemReal c_dx2 = d2b_dx2 * o4;
            const CfemReal c_dy2 = d2b_dy2 * o4;
            const CfemReal c_dxy = d2b_dxy * o4;

            for (CfemInt i = 0; i < 4; ++i) {
                d2x[i] = -c_dx2;
                d2y[i] = -c_dy2;
                d2z[i] = 0.0;
                dxz[i] = 0.0;
                dyz[i] = 0.0;
            }

            dxy[0] =  o4 - c_dxy;
            dxy[1] = -o4 - c_dxy;
            dxy[2] =  o4 - c_dxy;
            dxy[3] = -o4 - c_dxy;
        }
    };

    // =========================================================================
    // KERNEL QUAD P0 (1 nœud, constant par morceaux)
    // =========================================================================
    struct Q0QuadKernel {
        [[gnu::always_inline]] static inline void values(CfemReal* __restrict N) noexcept {
            N[0] = 1.0;
        }

        [[gnu::always_inline]] static inline void grads(CfemReal* __restrict dx, 
                                                        CfemReal* __restrict dy, 
                                                        CfemReal* __restrict dz) noexcept {
            dx[0] = 0.0;
            dy[0] = 0.0;
            dz[0] = 0.0;
        }

        [[gnu::always_inline]] static inline void hess(CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                                       CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept {
            d2x[0] = 0.0;
            d2y[0] = 0.0;
            d2z[0] = 0.0;
            dxy[0] = 0.0;
            dxz[0] = 0.0;
            dyz[0] = 0.0;
        }
    };
}

namespace cfem
{

    /**
     * @brief Q1 Reference Quadrilateral (4 nodes, bilinear)
     */
    class Q1ReferenceQuad : public ReferenceElement {
    public:
        Q1ReferenceQuad() {
            m_numNodes = 4;
            m_dimension = 2;
            m_nodes.resize(m_numNodes);

            m_nodes[0] = Vector3D{-1.0, -1.0, 0.0};
            m_nodes[1] = Vector3D{ 1.0, -1.0, 0.0};
            m_nodes[2] = Vector3D{ 1.0,  1.0, 0.0};
            m_nodes[3] = Vector3D{-1.0,  1.0, 0.0};
        }

        // =========================================================================
        // SCALAR INTERFACES
        // =========================================================================

        void evaluateShapesValues(const Vector3D &xi, ShapeInfo& info) const override {
            reference_elem::kernels::Q1QuadKernel::values(xi.x, xi.y, xi.z, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D &xi, ShapeInfo& info) const override {
            reference_elem::kernels::Q1QuadKernel::grads(xi.x, xi.y, xi.z, 
                                                         info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xi, ShapeInfo& info) const override {
            reference_elem::kernels::Q1QuadKernel::hess(xi.x, xi.y, xi.z, 
                                                        info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(),
                                                        info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        // =====================================================================
        // HPC BATCH INTERFACES
        // =====================================================================

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::Q1QuadKernel::values(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                              shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::Q1QuadKernel::grads(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                             shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::Q1QuadKernel::hess(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                            shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                            shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToQuadRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::quad4_face_nodes[iFace];
        }
    };


    /**
     * @brief Q2 Reference Quadrilateral (9 nodes, biquadratic)
     */
    class Q2ReferenceQuad : public ReferenceElement {
    public:
        Q2ReferenceQuad() {
            m_numNodes = 9;
            m_dimension = 2;
            m_nodes.resize(m_numNodes);
            
            // Sommets
            m_nodes[0] = Vector3D{-1.0, -1.0, 0.0};
            m_nodes[1] = Vector3D{ 1.0, -1.0, 0.0};
            m_nodes[2] = Vector3D{ 1.0,  1.0, 0.0};
            m_nodes[3] = Vector3D{-1.0,  1.0, 0.0};

            // Milieux d'arêtes
            m_nodes[4] = Vector3D{ 0.0, -1.0, 0.0};
            m_nodes[5] = Vector3D{ 1.0,  0.0, 0.0};
            m_nodes[6] = Vector3D{ 0.0,  1.0, 0.0};
            m_nodes[7] = Vector3D{-1.0,  0.0, 0.0};

            // Centre
            m_nodes[8] = Vector3D{ 0.0,  0.0, 0.0};
        }

        // =========================================================================
        // SCALAR INTERFACES
        // =========================================================================

        void evaluateShapesValues(const Vector3D &xi, ShapeInfo& info) const override {
            reference_elem::kernels::Q2QuadKernel::values(xi.x, xi.y, xi.z, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D &xi, ShapeInfo& info) const override {
            reference_elem::kernels::Q2QuadKernel::grads(xi.x, xi.y, xi.z, 
                                                         info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xi, ShapeInfo& info) const override {
            reference_elem::kernels::Q2QuadKernel::hess(xi.x, xi.y, xi.z, 
                                                        info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(),
                                                        info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        // =====================================================================
        // HPC BATCH INTERFACES
        // =====================================================================

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::Q2QuadKernel::values(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                              shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::Q2QuadKernel::grads(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                             shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::Q2QuadKernel::hess(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                            shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                            shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToQuadRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::quad9_face_nodes[iFace];
        }
    };

    /**
     * @brief Mini Reference Quadrilateral (5 nodes: 4 vertices + 1 bubble)
     */
    class MiniReferenceQuad : public ReferenceElement {
    public:
        MiniReferenceQuad() {
            m_numNodes = 5;
            m_dimension = 2;
            m_nodes.resize(m_numNodes);
            m_nodes[0] = {-1.0, -1.0, 0.0};
            m_nodes[1] = { 1.0, -1.0, 0.0};
            m_nodes[2] = { 1.0,  1.0, 0.0};
            m_nodes[3] = {-1.0,  1.0, 0.0};
            m_nodes[4] = { 0.0,  0.0, 0.0}; // Bulle
        }

        // =========================================================================
        // SCALAR INTERFACES
        // =========================================================================

        void evaluateShapesValues(const Vector3D& xi, ShapeInfo& info) const override {
            reference_elem::kernels::MiniQuadKernel::values(xi.x, xi.y, xi.z, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D& xi, ShapeInfo& info) const override {
            reference_elem::kernels::MiniQuadKernel::grads(xi.x, xi.y, xi.z, 
                                                           info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D& xi, ShapeInfo& info) const override {
            reference_elem::kernels::MiniQuadKernel::hess(xi.x, xi.y, xi.z, 
                                                          info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(),
                                                          info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        // =====================================================================
        // HPC BATCH INTERFACES
        // =====================================================================

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::MiniQuadKernel::values(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                                shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::MiniQuadKernel::grads(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                               shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::MiniQuadKernel::hess(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                              shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                              shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override {
            return reference_elem::utils::projectFacetRefCoordToQuadRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return reference_elem::utils::quad4_face_nodes[iFace];
        }
    };


    /**
     * @brief Q0 Reference Quadrilateral (1 node, piecewise constant)
     */
    class Q0ReferenceQuad : public ReferenceElement {
    public:
        Q0ReferenceQuad() {
            m_numNodes = 1;
            m_dimension = 2;

            // Un seul nœud au centre du carré de référence [-1, 1]x[-1, 1]
            m_nodes = {Vector3D{0.0, 0.0, 0.0}};
        }

        // =========================================================================
        // SCALAR INTERFACES
        // =========================================================================

        void evaluateShapesValues(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::Q0QuadKernel::values(info.values());
        }

        void evaluateShapesDerivatives(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::Q0QuadKernel::grads(info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::Q0QuadKernel::hess(info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(),
                                                        info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        // =====================================================================
        // HPC BATCH INTERFACES
        // =====================================================================

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::Q0QuadKernel::values(shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::Q0QuadKernel::grads(shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::Q0QuadKernel::hess(shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                            shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToQuadRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::quad1_face_nodes[iFace];
        }
    };
}
#endif