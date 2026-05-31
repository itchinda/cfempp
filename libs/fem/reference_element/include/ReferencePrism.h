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

#ifndef CFEM_REFERENCE_PRISM
#define CFEM_REFERENCE_PRISM

#include "ReferenceElement.h"
#include <span>
#include <vector>

namespace cfem::reference_elem::utils{

    inline Vector3D projectFacetRefCoordToPrismRefCoord(CfemInt localFacetId, 
                                                      const Vector3D& xiFacet) 
    {
        const CfemReal u = xiFacet.x;
        const CfemReal v = xiFacet.y;

        switch (localFacetId) {
            // Bases Triangles (u,v) dans [0,1]^2
            case 0: return Vector3D(v, u, -1.0);               // Face {0,2,1} : z=-1
            case 1: return Vector3D(u, v, 1.0);                // Face {3,4,5} : z=1
            
            // Faces Quads verticales (u,v) dans [-1,1]^2. Coordonnée z du prisme = v du quad.
            case 2: return Vector3D(0.5*(1.0+u), 0.0, v);      // Face {0,1,4,3} : y=0, x in [0,1]
            case 3: return Vector3D(0.5*(1.0-u), 0.5*(1.0+u), v); // Face {1,2,5,4} : oblique x+y=1
            case 4: return Vector3D(0.0, 0.5*(1.0-u), v);      // Face {2,0,3,5} : x=0, y in [0,1]
            default: CFEM_ERROR << "Invalid facet ID for Edge (Prism)";
        }
        return Vector3D{}; // Make the compiler happy
    }

    static const std::vector<std::vector<CfemInt>> prism6_face_nodes = {
        {0, 2, 1},       // Face 0 (Tri)
        {3, 4, 5},       // Face 1 (Tri)
        {0, 1, 4, 3},    // Face 2 (Quad)
        {1, 2, 5, 4},    // Face 3 (Quad)
        {2, 0, 3, 5}     // Face 4 (Quad)
    };

    static const std::vector<std::vector<CfemInt>> prism18_face_nodes = {
        {0, 2, 1, 8, 7, 6},                         // Face 0 (Tri6)
        {3, 4, 5, 9, 10, 11},                       // Face 1 (Tri6)
        {0, 1, 4, 3, 6, 13, 9, 12, 15},             // Face 2 (Quad9)
        {1, 2, 5, 4, 7, 14, 10, 13, 16},            // Face 3 (Quad9)
        {2, 0, 3, 5, 8, 12, 11, 14, 17}             // Face 4 (Quad9)
    };

    static const std::vector<std::vector<CfemInt>> prism1_face_nodes = {
        {0}, {0}, {0}, {0}, {0}
    };
}

namespace cfem::reference_elem::kernels {

    struct P1PrismKernel {
        static constexpr CfemReal o2 = 0.5;

        [[gnu::always_inline]] static inline void values(CfemReal x, CfemReal y, CfemReal z, 
                                                         CfemReal* __restrict N) noexcept 
        {
            const CfemReal Nt[3] = {1.0 - x - y, x, y};
            const CfemReal Nl[2] = {o2 * (1.0 - z), o2 * (1.0 + z)};

            for (CfemInt t = 0; t < 3; ++t) {
                for (CfemInt l = 0; l < 2; ++l) {
                    N[t + l * 3] = Nt[t] * Nl[l];
                }
            }
        }

        [[gnu::always_inline]] static inline void grads(CfemReal x, CfemReal y, CfemReal z, 
                                                        CfemReal* __restrict dx, 
                                                        CfemReal* __restrict dy, 
                                                        CfemReal* __restrict dz) noexcept 
        {
            const CfemReal Nt[3]      = {1.0 - x - y, x, y};
            const CfemReal dNt_dx[3]  = {-1.0, 1.0, 0.0};
            const CfemReal dNt_dy[3]  = {-1.0, 0.0, 1.0};
            
            const CfemReal Nl[2]      = {o2 * (1.0 - z), o2 * (1.0 + z)};
            const CfemReal dNl_dz[2]  = {-o2, o2};

            for (CfemInt t = 0; t < 3; ++t) {
                for (CfemInt l = 0; l < 2; ++l) {
                    const int node = t * 2 + l;
                    dx[node] = dNt_dx[t] * Nl[l];
                    dy[node] = dNt_dy[t] * Nl[l];
                    dz[node] = Nt[t] * dNl_dz[l];
                }
            }
        }

        [[gnu::always_inline]] static inline void hess(CfemReal x, CfemReal y, CfemReal z, 
                                                       CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                                       CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept 
        {
            const CfemReal dNt_dx[3]  = {-1.0, 1.0, 0.0};
            const CfemReal dNt_dy[3]  = {-1.0, 0.0, 1.0};
            const CfemReal dNl_dz[2]  = {-o2, o2};

            for (CfemInt t = 0; t < 3; ++t) {
                for (CfemInt l = 0; l < 2; ++l) {
                    const int node = t * 2 + l;
                    // Pure second derivatives are zero
                    d2x[node] = d2y[node] = d2z[node] = dxy[node] = 0.0;
                    // Mixed derivatives with z
                    dxz[node] = dNt_dx[t] * dNl_dz[l];
                    dyz[node] = dNt_dy[t] * dNl_dz[l];
                }
            }
        }
    };

    struct P2PrismKernel {
        // Mapping strict selon la connectivité de P2ReferencePrism
        // t (Triangle P2) : 0..5
        // l (Ligne P2)    : 0 (-1), 1 (0), 2 (+1)
        static constexpr CfemInt node_map[6][3] = {
            {0, 12, 3}, {1, 13, 4}, {2, 14, 5},
            {6, 15, 9}, {7, 16, 10}, {8, 17, 11}
        };

        [[gnu::always_inline]] static inline void eval2D_tri_all(CfemReal x, CfemReal y, 
                                          CfemReal* N, CfemReal* dNx, CfemReal* dNy, 
                                          CfemReal* d2Nx, CfemReal* d2Ny, CfemReal* d2Nxy) noexcept
        {
            const CfemReal L1 = x, L2 = y, L0 = 1.0 - x - y;

            N[0] = L0 * (2.0*L0 - 1.0);  N[1] = L1 * (2.0*L1 - 1.0);  N[2] = L2 * (2.0*L2 - 1.0);
            N[3] = 4.0 * L0 * L1;        N[4] = 4.0 * L1 * L2;        N[5] = 4.0 * L2 * L0;

            dNx[0] = -3.0 + 4.0*L1 + 4.0*L2; dNy[0] = -3.0 + 4.0*L1 + 4.0*L2;
            dNx[1] =  4.0*L1 - 1.0;          dNy[1] =  0.0;
            dNx[2] =  0.0;                   dNy[2] =  4.0*L2 - 1.0;
            dNx[3] =  4.0 - 8.0*L1 - 4.0*L2; dNy[3] = -4.0*L1;
            dNx[4] =  4.0*L2;                dNy[4] =  4.0*L1;
            dNx[5] = -4.0*L2;                dNy[5] =  4.0 - 4.0*L1 - 8.0*L2;

            d2Nx[0] =  4.0; d2Ny[0] =  4.0; d2Nxy[0] =  4.0;
            d2Nx[1] =  4.0; d2Ny[1] =  0.0; d2Nxy[1] =  0.0;
            d2Nx[2] =  0.0; d2Ny[2] =  4.0; d2Nxy[2] =  0.0;
            d2Nx[3] = -8.0; d2Ny[3] =  0.0; d2Nxy[3] = -4.0;
            d2Nx[4] =  0.0; d2Ny[4] =  0.0; d2Nxy[4] =  4.0;
            d2Nx[5] =  0.0; d2Ny[5] = -8.0; d2Nxy[5] = -4.0;
        }

        [[gnu::always_inline]] static inline void eval1D_line_all(CfemReal z, CfemReal* N, CfemReal* dNz, CfemReal* d2Nz) noexcept
        {
            const CfemReal z2 = z * z;
            N[0] = 0.5 * (z2 - z); dNz[0] = z - 0.5;  d2Nz[0] = 1.0;
            N[1] = 1.0 - z2;       dNz[1] = -2.0 * z; d2Nz[1] = -2.0;
            N[2] = 0.5 * (z2 + z); dNz[2] = z + 0.5;  d2Nz[2] = 1.0;
        }

        [[gnu::always_inline]] static inline void values(CfemReal x, CfemReal y, CfemReal z, 
                                                         CfemReal* __restrict N) noexcept 
        {
            CfemReal Nt[6], dNtx[6], dNty[6], d2Ntx[6], d2Nty[6], d2Ntxy[6];
            CfemReal Nl[3], dNlz[3], d2Nlz[3];

            eval2D_tri_all(x, y, Nt, dNtx, dNty, d2Ntx, d2Nty, d2Ntxy);
            eval1D_line_all(z, Nl, dNlz, d2Nlz);

            for (CfemInt t = 0; t < 6; ++t) {
                for (CfemInt l = 0; l < 3; ++l) {
                    N[node_map[t][l]] = Nt[t] * Nl[l];
                }
            }
        }

        [[gnu::always_inline]] static inline void grads(CfemReal x, CfemReal y, CfemReal z, 
                                                        CfemReal* __restrict dx, 
                                                        CfemReal* __restrict dy, 
                                                        CfemReal* __restrict dz) noexcept 
        {
            CfemReal Nt[6], dNtx[6], dNty[6], d2Ntx[6], d2Nty[6], d2Ntxy[6];
            CfemReal Nl[3], dNlz[3], d2Nlz[3];

            eval2D_tri_all(x, y, Nt, dNtx, dNty, d2Ntx, d2Nty, d2Ntxy);
            eval1D_line_all(z, Nl, dNlz, d2Nlz);

            for (CfemInt t = 0; t < 6; ++t) {
                for (CfemInt l = 0; l < 3; ++l) {
                    const int node = node_map[t][l];
                    dx[node] = dNtx[t] * Nl[l];
                    dy[node] = dNty[t] * Nl[l];
                    dz[node] = Nt[t]   * dNlz[l];
                }
            }
        }

        [[gnu::always_inline]] static inline void hess(CfemReal x, CfemReal y, CfemReal z, 
                                                       CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                                       CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept 
        {
            CfemReal Nt[6], dNtx[6], dNty[6], d2Ntx[6], d2Nty[6], d2Ntxy[6];
            CfemReal Nl[3], dNlz[3], d2Nlz[3];

            eval2D_tri_all(x, y, Nt, dNtx, dNty, d2Ntx, d2Nty, d2Ntxy);
            eval1D_line_all(z, Nl, dNlz, d2Nlz);

            for (CfemInt t = 0; t < 6; ++t) {
                for (CfemInt l = 0; l < 3; ++l) {
                    const int node = node_map[t][l];
                    
                    d2x[node] = d2Ntx[t]  * Nl[l];
                    d2y[node] = d2Nty[t]  * Nl[l];
                    d2z[node] = Nt[t]     * d2Nlz[l];
                    
                    dxy[node] = d2Ntxy[t] * Nl[l];
                    dxz[node] = dNtx[t]   * dNlz[l];
                    dyz[node] = dNty[t]   * dNlz[l];
                }
            }
        }
    };

    struct MiniPrismKernel {
        [[gnu::always_inline]] static inline void values(CfemReal x, CfemReal y, CfemReal z, 
                                                         CfemReal* __restrict N) noexcept 
        {
            constexpr CfemReal o6 = 1.0 / 6.0;

            const CfemReal L[3] = {1.0 - x - y, x, y};
            const CfemReal M[2] = {0.5 * (1.0 - z), 0.5 * (1.0 + z)};

            // Calcul de la Bulle (Nœud 6)
            const CfemReal triB = 27.0 * L[0] * L[1] * L[2];
            const CfemReal zB = 1.0 - z * z;
            const CfemReal B = triB * zB;

            N[6] = B;

            // Fonctions des Sommets (0-5) avec correction
            const CfemReal correction = B * o6;
            for (CfemInt i = 0; i < 6; ++i) {
                const int iT = i % 3;
                const int iZ = i / 3;
                N[i] = L[iT] * M[iZ] - correction;
            }
        }

        [[gnu::always_inline]] static inline void grads(CfemReal x, CfemReal y, CfemReal z, 
                                                        CfemReal* __restrict dx, 
                                                        CfemReal* __restrict dy, 
                                                        CfemReal* __restrict dz) noexcept 
        {
            constexpr CfemReal o6 = 1.0 / 6.0;

            const CfemReal L[3] = {1.0 - x - y, x, y};
            const CfemReal M[2] = {0.5 * (1.0 - z), 0.5 * (1.0 + z)};
            const CfemReal dM[2] = {-0.5, 0.5};

            // Gradients de la bulle
            const CfemReal triB = 27.0 * L[0] * L[1] * L[2];
            const CfemReal zB = 1.0 - z * z;

            const CfemReal gBx = 27.0 * y * (1.0 - 2.0*x - y) * zB;
            const CfemReal gBy = 27.0 * x * (1.0 - x - 2.0*y) * zB;
            const CfemReal gBz = triB * (-2.0 * z);

            dx[6] = gBx;
            dy[6] = gBy;
            dz[6] = gBz;

            // Gradients des sommets avec correction
            const CfemReal c_dx = gBx * o6;
            const CfemReal c_dy = gBy * o6;
            const CfemReal c_dz = gBz * o6;

            for (CfemInt i = 0; i < 6; ++i) {
                const int iT = i % 3;
                const int iZ = i / 3;
                const CfemReal drL = (iT == 0 ? -1.0 : (iT == 1 ? 1.0 : 0.0));
                const CfemReal dsL = (iT == 0 ? -1.0 : (iT == 2 ? 1.0 : 0.0));

                dx[i] = drL * M[iZ] - c_dx;
                dy[i] = dsL * M[iZ] - c_dy;
                dz[i] = L[iT] * dM[iZ] - c_dz;
            }
        }

        [[gnu::always_inline]] static inline void hess(CfemReal x, CfemReal y, CfemReal z, 
                                                       CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                                       CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept 
        {
            constexpr CfemReal o6 = 1.0 / 6.0;

            const CfemReal L[3] = {1.0 - x - y, x, y};
            const CfemReal dM[2] = {-0.5, 0.5};

            const CfemReal triB = 27.0 * L[0] * L[1] * L[2];
            const CfemReal zB = 1.0 - z * z;

            // Hessiens de la bulle
            const CfemReal hBxx = -54.0 * y * zB;
            const CfemReal hByy = -54.0 * x * zB;
            const CfemReal hBzz = -2.0 * triB;
            const CfemReal hBxy =  27.0 * (1.0 - 2.0*x - 2.0*y) * zB;
            const CfemReal hBxz = -54.0 * y * z * (1.0 - 2.0*x - y); 
            const CfemReal hByz = -54.0 * x * z * (1.0 - x - 2.0*y); 

            d2x[6] = hBxx;
            d2y[6] = hByy;
            d2z[6] = hBzz;
            dxy[6] = hBxy;
            dxz[6] = hBxz;
            dyz[6] = hByz;

            // Hessiens des sommets avec correction
            const CfemReal c_xx = hBxx * o6;
            const CfemReal c_yy = hByy * o6;
            const CfemReal c_zz = hBzz * o6;
            const CfemReal c_xy = hBxy * o6;
            const CfemReal c_xz = hBxz * o6;
            const CfemReal c_yz = hByz * o6;

            for (CfemInt i = 0; i < 6; ++i) {
                const int iT = i % 3;
                const int iZ = i / 3;
                const CfemReal drL = (iT == 0 ? -1.0 : (iT == 1 ? 1.0 : 0.0));
                const CfemReal dsL = (iT == 0 ? -1.0 : (iT == 2 ? 1.0 : 0.0));

                d2x[i] = -c_xx;
                d2y[i] = -c_yy;
                d2z[i] = -c_zz;
                dxy[i] = -c_xy;
                dxz[i] = drL * dM[iZ] - c_xz;
                dyz[i] = dsL * dM[iZ] - c_yz;
            }
        }
    };

    // =========================================================================
    // KERNEL PRISME P0 (1 nœud, constant)
    // =========================================================================
    struct P0PrismKernel {
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

}

namespace cfem{

    /**
     * @brief P1 Reference Prism (Wedge) (6 nodes, P1)
     */
    class P1ReferencePrism : public ReferenceElement {
    public:
        P1ReferencePrism() {
            m_numNodes = 6;
            m_dimension = 3;
            m_nodes.resize(m_numNodes);
            // VTK order : Tri Z=-1 (0,1,2), puis Z=1 (3,4,5)
            m_nodes[0] = {0,0,-1}; m_nodes[1] = {1,0,-1}; m_nodes[2] = {0,1,-1};
            m_nodes[3] = {0,0, 1}; m_nodes[4] = {1,0, 1}; m_nodes[5] = {0,1, 1};
        }

        // =========================================================================
        // SCALAR INTERFACES
        // =========================================================================

        void evaluateShapesValues(const Vector3D& xi, ShapeInfo& info) const override {
            reference_elem::kernels::P1PrismKernel::values(xi.x, xi.y, xi.z, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D& xi, ShapeInfo& info) const override {
            reference_elem::kernels::P1PrismKernel::grads(xi.x, xi.y, xi.z, 
                                                          info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D& xi, ShapeInfo& info) const override {
            reference_elem::kernels::P1PrismKernel::hess(xi.x, xi.y, xi.z, 
                                                         info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(), 
                                                         info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        // =====================================================================
        // HPC BATCH INTERFACES
        // =====================================================================

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P1PrismKernel::values(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                                shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P1PrismKernel::grads(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                              shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P1PrismKernel::hess(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                             shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                             shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToPrismRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::prism6_face_nodes[iFace];
        }
    };


    /**
     * @brief P2 Reference Prism (Wedge) (18 nodes, P2)
     */
    class P2ReferencePrism : public ReferenceElement {
    public:
        P2ReferencePrism() {
            m_numNodes = 18;
            m_dimension = 3;
            m_nodes.resize(m_numNodes);
            
            // 1. Sommets (0-5)
            m_nodes[0] = {0.0, 0.0, -1.0}; m_nodes[1] = {1.0, 0.0, -1.0}; m_nodes[2] = {0.0, 1.0, -1.0};
            m_nodes[3] = {0.0, 0.0,  1.0}; m_nodes[4] = {1.0, 0.0,  1.0}; m_nodes[5] = {0.0, 1.0,  1.0};

            // 2. Milieux d'arêtes horizontales (6-11)
            m_nodes[6] = {0.5, 0.0, -1.0}; m_nodes[7] = {0.5, 0.5, -1.0}; m_nodes[8] = {0.0, 0.5, -1.0};
            m_nodes[9] = {0.5, 0.0,  1.0}; m_nodes[10]= {0.5, 0.5,  1.0}; m_nodes[11]= {0.0, 0.5,  1.0};

            // 3. Milieux d'arêtes verticales (12-14)
            m_nodes[12]= {0.0, 0.0,  0.0}; m_nodes[13]= {1.0, 0.0,  0.0}; m_nodes[14]= {0.0, 1.0,  0.0};

            // 4. Centres des faces quadrilatérales (15-17)
            m_nodes[15]= {0.5, 0.0,  0.0}; m_nodes[16]= {0.5, 0.5,  0.0}; m_nodes[17]= {0.0, 0.5,  0.0};
        }

        // =========================================================================
        // SCALAR INTERFACES
        // =========================================================================

        void evaluateShapesValues(const Vector3D& xi, ShapeInfo& info) const override {
            reference_elem::kernels::P2PrismKernel::values(xi.x, xi.y, xi.z, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D& xi, ShapeInfo& info) const override {
            reference_elem::kernels::P2PrismKernel::grads(xi.x, xi.y, xi.z, 
                                                          info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D& xi, ShapeInfo& info) const override {
            reference_elem::kernels::P2PrismKernel::hess(xi.x, xi.y, xi.z, 
                                                         info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(), 
                                                         info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        // =====================================================================
        // HPC BATCH INTERFACES
        // =====================================================================

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P2PrismKernel::values(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                               shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P2PrismKernel::grads(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                              shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P2PrismKernel::hess(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                             shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                             shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToPrismRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::prism18_face_nodes[iFace];
        }
    };

    /**
     * @brief Mini Reference Prism (7 nodes: 6 vertices + 1 bubble)
     */
    class MiniReferencePrism : public ReferenceElement {
    public:
        MiniReferencePrism() {
            m_numNodes = 7;
            m_dimension = 3;
            m_nodes.resize(m_numNodes);
            m_nodes[0] = {0.0, 0.0, -1.0}; m_nodes[1] = {1.0, 0.0, -1.0}; m_nodes[2] = {0.0, 1.0, -1.0};
            m_nodes[3] = {0.0, 0.0,  1.0}; m_nodes[4] = {1.0, 0.0,  1.0}; m_nodes[5] = {0.0, 1.0,  1.0};
            m_nodes[6] = {1.0/3.0, 1.0/3.0, 0.0}; // Bulle au centre de gravité
        }

        // =========================================================================
        // SCALAR INTERFACES
        // =========================================================================

        void evaluateShapesValues(const Vector3D& xi, ShapeInfo& info) const override {
            reference_elem::kernels::MiniPrismKernel::values(xi.x, xi.y, xi.z, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D& xi, ShapeInfo& info) const override {
            reference_elem::kernels::MiniPrismKernel::grads(xi.x, xi.y, xi.z, 
                                                            info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D& xi, ShapeInfo& info) const override {
            reference_elem::kernels::MiniPrismKernel::hess(xi.x, xi.y, xi.z, 
                                                           info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(), 
                                                           info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        // =====================================================================
        // HPC BATCH INTERFACES
        // =====================================================================

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::MiniPrismKernel::values(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                                 shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::MiniPrismKernel::grads(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                                shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::MiniPrismKernel::hess(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                               shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                               shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToPrismRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::prism6_face_nodes[iFace];
        }
    };

    /**
     * @brief P0 Reference Prism (1 node, piecewise constant)
     */
    class P0ReferencePrism : public ReferenceElement {
    public:
        P0ReferencePrism() {
            m_numNodes = 1;
            m_dimension = 3;
            // Un seul nœud au centre de gravité
            m_nodes = { Vector3D{1.0/3.0, 1.0/3.0, 0.0} };
        }

        // =========================================================================
        // SCALAR INTERFACES
        // =========================================================================

        void evaluateShapesValues(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::P0PrismKernel::values(info.values());
        }

        void evaluateShapesDerivatives(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::P0PrismKernel::grads(info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D&, ShapeInfo& info) const override {
            reference_elem::kernels::P0PrismKernel::hess(info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(), 
                                                         info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        // =====================================================================
        // HPC BATCH INTERFACES
        // =====================================================================

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P0PrismKernel::values(shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P0PrismKernel::grads(shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, ShapeInfoBatch& shapeBatch) const override {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P0PrismKernel::hess(shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                             shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToPrismRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::prism1_face_nodes[iFace];
        }
    };
}
#endif

