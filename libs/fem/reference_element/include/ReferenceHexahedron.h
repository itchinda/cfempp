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

#ifndef CFEM_REFERENCE_HEXAHEDRON
#define CFEM_REFERENCE_HEXAHEDRON

#include <array>
#include <cstring>
#include <span>

#include "ReferenceElement.h"

namespace cfem::reference_elem::utils
{

    inline Vector3D projectFacetRefCoordToHexaRefCoord(CfemInt localFacetId,
                                                       const Vector3D &xiFacet)
    {
        const CfemReal u = xiFacet.x;
        const CfemReal v = xiFacet.y;

        switch (localFacetId)
        {
        case 0:
            return Vector3D(v, u, -1.0); // Face {0,3,2,1} : z=-1
        case 1:
            return Vector3D(u, v, 1.0); // Face {4,5,6,7} : z=1
        case 2:
            return Vector3D(u, -1.0, v); // Face {0,1,5,4} : y=-1
        case 3:
            return Vector3D(1.0, u, v); // Face {1,2,6,5} : x=1
        case 4:
            return Vector3D(-u, 1.0, v); // Face {2,3,7,6} : y=1
        case 5:
            return Vector3D(-1.0, -u, v); // Face {3,0,4,7} : x=-1
        default:
            CFEM_ERROR << "Invalid facet ID for Hexa";
        }
        return Vector3D{}; // Make the compiler happy
    }

    static constexpr std::array<std::array<int, 4>, 6> hex8_face_nodes = {{
        {0, 3, 2, 1}, // Face 0 (Bottom)
        {4, 5, 6, 7}, // Face 1 (Top)
        {0, 1, 5, 4}, // Face 2 (Front)
        {1, 2, 6, 5}, // Face 3 (Right)
        {2, 3, 7, 6}, // Face 4 (Back)
        {3, 0, 4, 7}  // Face 5 (Left)
    }};

    static constexpr std::array<std::array<CfemInt, 9>, 6> hex27_face_nodes = {{
        {0, 3, 2, 1, 11, 10, 9, 8, 24},   // Face 0
        {4, 5, 6, 7, 12, 13, 14, 15, 25}, // Face 1
        {0, 1, 5, 4, 8, 17, 12, 16, 22},  // Face 2
        {1, 2, 6, 5, 9, 18, 13, 17, 21},  // Face 3
        {2, 3, 7, 6, 10, 19, 14, 18, 23}, // Face 4
        {3, 0, 4, 7, 11, 16, 15, 19, 20}  // Face 5
    }};

    static constexpr std::array<std::array<CfemInt, 1>, 6> hex1_face_nodes = {{
        {0}, {0}, {0}, {0}, {0}, {0}
    }};
}

namespace cfem::reference_elem::kernels {
    // =========================================================================
    // Q0 HEX KERNEL (1 node, constant)
    // =========================================================================
    struct Q0HexaKernel {
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
    // MINI HEX KERNEL  (9 nodes : 8 vertices Q1 + 1 bubble)
    // =========================================================================
    struct MiniHexaKernel {
    private:
        static constexpr CfemReal xr[8] = {-1.0, 1.0, 1.0, -1.0, -1.0,  1.0, 1.0, -1.0};
        static constexpr CfemReal ys[8] = {-1.0, -1.0, 1.0,  1.0, -1.0, -1.0, 1.0,  1.0};
        static constexpr CfemReal zt[8] = {-1.0, -1.0,-1.0, -1.0,  1.0,  1.0, 1.0,  1.0};
        static constexpr CfemReal o8 = 0.125;

    public:
        [[gnu::always_inline]] static inline void values(CfemReal x, CfemReal y, CfemReal z, 
                                                         CfemReal* __restrict N) noexcept 
        {
            const CfemReal b = (1.0 - x * x) * (1.0 - y * y) * (1.0 - z * z);
            N[8] = b;

            const CfemReal correction = b * o8;
            for (CfemInt i = 0; i < 8; ++i) {
                const CfemReal px = 1.0 + xr[i] * x;
                const CfemReal py = 1.0 + ys[i] * y;
                const CfemReal pz = 1.0 + zt[i] * z;
                N[i] = (o8 * px * py * pz) - correction;
            }
        }

        [[gnu::always_inline]] static inline void grads(CfemReal x, CfemReal y, CfemReal z, 
                                                        CfemReal* __restrict dx, 
                                                        CfemReal* __restrict dy, 
                                                        CfemReal* __restrict dz) noexcept 
        {
            const CfemReal mx2 = 1.0 - x * x;
            const CfemReal my2 = 1.0 - y * y;
            const CfemReal mz2 = 1.0 - z * z;

            const CfemReal db_dx = -2.0 * x * my2 * mz2;
            const CfemReal db_dy = -2.0 * y * mx2 * mz2;
            const CfemReal db_dz = -2.0 * z * mx2 * my2;

            dx[8] = db_dx; dy[8] = db_dy; dz[8] = db_dz;

            const CfemReal corr_dx = db_dx * o8;
            const CfemReal corr_dy = db_dy * o8;
            const CfemReal corr_dz = db_dz * o8;

            for (CfemInt i = 0; i < 8; ++i) {
                const CfemReal px = 1.0 + xr[i] * x;
                const CfemReal py = 1.0 + ys[i] * y;
                const CfemReal pz = 1.0 + zt[i] * z;

                dx[i] = (o8 * xr[i] * py * pz) - corr_dx;
                dy[i] = (o8 * ys[i] * px * pz) - corr_dy;
                dz[i] = (o8 * zt[i] * px * py) - corr_dz;
            }
        }

        [[gnu::always_inline]] static inline void hess(CfemReal x, CfemReal y, CfemReal z, 
                                                       CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                                       CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept 
        {
            const CfemReal mx2 = 1.0 - x * x;
            const CfemReal my2 = 1.0 - y * y;
            const CfemReal mz2 = 1.0 - z * z;

            const CfemReal d2b_dx2 = -2.0 * my2 * mz2;
            const CfemReal d2b_dy2 = -2.0 * mx2 * mz2;
            const CfemReal d2b_dz2 = -2.0 * mx2 * my2;
            const CfemReal d2b_dxy =  4.0 * x * y * mz2;
            const CfemReal d2b_dxz =  4.0 * x * z * my2;
            const CfemReal d2b_dyz =  4.0 * y * z * mx2;

            d2x[8] = d2b_dx2; d2y[8] = d2b_dy2; d2z[8] = d2b_dz2;
            dxy[8] = d2b_dxy; dxz[8] = d2b_dxz; dyz[8] = d2b_dyz;

            const CfemReal corr_dx2 = d2b_dx2 * o8;
            const CfemReal corr_dy2 = d2b_dy2 * o8;
            const CfemReal corr_dz2 = d2b_dz2 * o8;
            const CfemReal corr_dxy = d2b_dxy * o8;
            const CfemReal corr_dxz = d2b_dxz * o8;
            const CfemReal corr_dyz = d2b_dyz * o8;

            for (CfemInt i = 0; i < 8; ++i) {
                const CfemReal px = 1.0 + xr[i] * x;
                const CfemReal py = 1.0 + ys[i] * y;
                const CfemReal pz = 1.0 + zt[i] * z;

                d2x[i] = -corr_dx2; // Pure Q1 second derivatives are 0
                d2y[i] = -corr_dy2;
                d2z[i] = -corr_dz2;
                dxy[i] = (o8 * xr[i] * ys[i] * pz) - corr_dxy;
                dxz[i] = (o8 * xr[i] * py * zt[i]) - corr_dxz;
                dyz[i] = (o8 * px * ys[i] * zt[i]) - corr_dyz;
            }
        }
    };

    // =========================================================================
    // Q1 HEX KERNEL (8 nodes, trilinear)
    // =========================================================================
    struct Q1HexaKernel {
        [[gnu::always_inline]] static inline void values(CfemReal x, CfemReal y, CfemReal z, 
                                                         CfemReal* __restrict N) noexcept 
        {
            const CfemReal xm = 1.0 - x; const CfemReal xp = 1.0 + x;
            const CfemReal ym = 1.0 - y; const CfemReal yp = 1.0 + y;
            const CfemReal zm = 1.0 - z; const CfemReal zp = 1.0 + z;

            constexpr CfemReal o8 = 0.125; // 1.0 / 8.0

            N[0] = o8 * xm * ym * zm;
            N[1] = o8 * xp * ym * zm;
            N[2] = o8 * xp * yp * zm;
            N[3] = o8 * xm * yp * zm;
            N[4] = o8 * xm * ym * zp;
            N[5] = o8 * xp * ym * zp;
            N[6] = o8 * xp * yp * zp;
            N[7] = o8 * xm * yp * zp;
        }

        [[gnu::always_inline]] static inline void grads(CfemReal x, CfemReal y, CfemReal z, 
                                                        CfemReal* __restrict dx, 
                                                        CfemReal* __restrict dy, 
                                                        CfemReal* __restrict dz) noexcept 
        {
            const CfemReal xm = 1.0 - x; const CfemReal xp = 1.0 + x;
            const CfemReal ym = 1.0 - y; const CfemReal yp = 1.0 + y;
            const CfemReal zm = 1.0 - z; const CfemReal zp = 1.0 + z;

            constexpr CfemReal o8 = 0.125;

            // d / dx
            dx[0] = -o8 * ym * zm; dx[1] =  o8 * ym * zm;
            dx[2] =  o8 * yp * zm; dx[3] = -o8 * yp * zm;
            dx[4] = -o8 * ym * zp; dx[5] =  o8 * ym * zp;
            dx[6] =  o8 * yp * zp; dx[7] = -o8 * yp * zp;

            // d / dy
            dy[0] = -o8 * xm * zm; dy[1] = -o8 * xp * zm;
            dy[2] =  o8 * xp * zm; dy[3] =  o8 * xm * zm;
            dy[4] = -o8 * xm * zp; dy[5] = -o8 * xp * zp;
            dy[6] =  o8 * xp * zp; dy[7] =  o8 * xm * zp;

            // d / dz
            dz[0] = -o8 * xm * ym; dz[1] = -o8 * xp * ym;
            dz[2] = -o8 * xp * yp; dz[3] = -o8 * xm * yp;
            dz[4] =  o8 * xm * ym; dz[5] =  o8 * xp * ym;
            dz[6] =  o8 * xp * yp; dz[7] =  o8 * xm * yp;
        }

        [[gnu::always_inline]] static inline void hess(CfemReal x, CfemReal y, CfemReal z, 
                                                       CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                                       CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept 
        {
            const CfemReal xm = 1.0 - x; const CfemReal xp = 1.0 + x;
            const CfemReal ym = 1.0 - y; const CfemReal yp = 1.0 + y;
            const CfemReal zm = 1.0 - z; const CfemReal zp = 1.0 + z;

            constexpr CfemReal o8 = 0.125;

            // Pure second derivatives are strictly zero for trilinear fields
            for (CfemInt i = 0; i < 8; ++i) {
                d2x[i] = 0.0;
                d2y[i] = 0.0;
                d2z[i] = 0.0;
            }

            // Cross derivatives: d2N / (dx * dy)
            dxy[0] =  o8 * zm; dxy[1] = -o8 * zm; dxy[2] =  o8 * zm; dxy[3] = -o8 * zm;
            dxy[4] =  o8 * zp; dxy[5] = -o8 * zp; dxy[6] =  o8 * zp; dxy[7] = -o8 * zp;

            // Cross derivatives: d2N / (dx * dz)
            dxz[0] =  o8 * ym; dxz[1] = -o8 * ym; dxz[2] = -o8 * yp; dxz[3] =  o8 * yp;
            dxz[4] = -o8 * ym; dxz[5] =  o8 * ym; dxz[6] =  o8 * yp; dxz[7] = -o8 * yp;

            // Cross derivatives: d2N / (dy * dz)
            dyz[0] =  o8 * xm; dyz[1] =  o8 * xp; dyz[2] = -o8 * xp; dyz[3] = -o8 * xm;
            dyz[4] = -o8 * xm; dyz[5] = -o8 * xp; dyz[6] =  o8 * xp; dyz[7] =  o8 * xm;
        }
    };

    // =========================================================================
    // Q2 HEX KERNEL (27 nodes, tri-quadratic)
    // =========================================================================
    struct Q2HexaKernel {
    public:
        // i (X) parcourt : 0 (-1), 1 (0), 2 (+1)
        // j (Y) parcourt : 0 (-1), 1 (0), 2 (+1)
        // k (Z) parcourt : 0 (-1), 1 (0), 2 (+1)
        static constexpr CfemInt node_map[3][3][3] = {
            { {0, 16, 4},   {11, 20, 15}, {3, 19, 7} },
            { {8, 22, 12},  {24, 26, 25}, {10, 23, 14} },
            { {1, 17, 5},   {9, 21, 13},  {2, 18, 6} }
        };

    private:
        static inline void eval1D_shape(CfemReal u, CfemReal *N) noexcept {
            const CfemReal u2 = u * u;
            N[0] = 0.5 * (u2 - u);
            N[1] = 1.0 - u2;
            N[2] = 0.5 * (u2 + u);
        }

        static inline void eval1D_deriv(CfemReal u, CfemReal *N, CfemReal *dN) noexcept {
            const CfemReal u2 = u * u;
            N[0] = 0.5 * (u2 - u); dN[0] = u - 0.5;
            N[1] = 1.0 - u2;       dN[1] = -2.0 * u;
            N[2] = 0.5 * (u2 + u); dN[2] = u + 0.5;
        }

        static inline void eval1D_hessian(CfemReal u, CfemReal *N, CfemReal *dN, CfemReal *d2N) noexcept {
            const CfemReal u2 = u * u;
            N[0] = 0.5 * (u2 - u); dN[0] = u - 0.5;  d2N[0] = 1.0;
            N[1] = 1.0 - u2;       dN[1] = -2.0 * u; d2N[1] = -2.0;
            N[2] = 0.5 * (u2 + u); dN[2] = u + 0.5;  d2N[2] = 1.0;
        }

    public:
        [[gnu::always_inline]] static inline void values(CfemReal x, CfemReal y, CfemReal z, 
                                                         CfemReal* __restrict N) noexcept 
        {
            CfemReal Nx[3], Ny[3], Nz[3];
            eval1D_shape(x, Nx);
            eval1D_shape(y, Ny);
            eval1D_shape(z, Nz);

            for (CfemInt i = 0; i < 3; ++i) {
                for (CfemInt j = 0; j < 3; ++j) {
                    for (CfemInt k = 0; k < 3; ++k) {
                        N[node_map[i][j][k]] = Nx[i] * Ny[j] * Nz[k];
                    }
                }
            }
        }

        [[gnu::always_inline]] static inline void grads(CfemReal x, CfemReal y, CfemReal z, 
                                                        CfemReal* __restrict dx, 
                                                        CfemReal* __restrict dy, 
                                                        CfemReal* __restrict dz) noexcept 
        {
            CfemReal Nx[3], Ny[3], Nz[3], dNx[3], dNy[3], dNz[3];
            eval1D_deriv(x, Nx, dNx);
            eval1D_deriv(y, Ny, dNy);
            eval1D_deriv(z, Nz, dNz);

            for (CfemInt i = 0; i < 3; ++i) {
                for (CfemInt j = 0; j < 3; ++j) {
                    for (CfemInt k = 0; k < 3; ++k) {
                        const CfemInt node = node_map[i][j][k];
                        dx[node] = dNx[i] * Ny[j] * Nz[k];
                        dy[node] = Nx[i] * dNy[j] * Nz[k];
                        dz[node] = Nx[i] * Ny[j] * dNz[k];
                    }
                }
            }
        }

        [[gnu::always_inline]] static inline void hess(CfemReal x, CfemReal y, CfemReal z, 
                                                       CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                                       CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept 
        {
            CfemReal Nx[3], Ny[3], Nz[3], dNx[3], dNy[3], dNz[3], d2Nx[3], d2Ny[3], d2Nz[3];
            eval1D_hessian(x, Nx, dNx, d2Nx);
            eval1D_hessian(y, Ny, dNy, d2Ny);
            eval1D_hessian(z, Nz, dNz, d2Nz);

            for (CfemInt i = 0; i < 3; ++i) {
                for (CfemInt j = 0; j < 3; ++j) {
                    for (CfemInt k = 0; k < 3; ++k) {
                        const CfemInt node = node_map[i][j][k];
                        
                        d2x[node] = d2Nx[i] * Ny[j] * Nz[k];
                        d2y[node] = Nx[i] * d2Ny[j] * Nz[k];
                        d2z[node] = Nx[i] * Ny[j] * d2Nz[k];

                        dxy[node] = dNx[i] * dNy[j] * Nz[k];
                        dxz[node] = dNx[i] * Ny[j] * dNz[k];
                        dyz[node] = Nx[i] * dNy[j] * dNz[k];
                    }
                }
            }
        }
    };
}

namespace cfem
{
    /**
     * @brief Constant Reference Hexahedron (1 node, Q0)
     */
    class Q0ReferenceHexahedron : public ReferenceElement
    {
    public:
        Q0ReferenceHexahedron()
        {
            m_numNodes = 1;
            m_dimension = 3;
            // Single node at the reference center [-1, 1]^3
            m_nodes = {Vector3D{0.0, 0.0, 0.0}};
        }

        // =========================================================================
        // SCALAR INTERFACES
        // =========================================================================

        void evaluateShapesValues(const Vector3D &, ShapeInfo& info) const override
        {
            reference_elem::kernels::Q0HexaKernel::values(info.values());
        }

        void evaluateShapesDerivatives(const Vector3D &, ShapeInfo& info) const override
        {
            reference_elem::kernels::Q0HexaKernel::grads(info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D &, ShapeInfo& info) const override
        {
            reference_elem::kernels::Q0HexaKernel::hess(info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(),
                                                        info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        // =====================================================================
        // HPC BATCH INTERFACES (Full Memcpy Broadcasting)
        // =====================================================================

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override
        {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            if (nQp == 0) return;

            reference_elem::kernels::Q0HexaKernel::values(shapeBatch.values(0));

            const size_t bytesToCopy = m_numNodes * sizeof(CfemReal);
            for (CfemInt q = 1; q < nQp; ++q) {
                std::memcpy(shapeBatch.values(q), shapeBatch.values(0), bytesToCopy);
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override
        {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            if (nQp == 0) return;

            reference_elem::kernels::Q0HexaKernel::grads(shapeBatch.dN_dxi(0), shapeBatch.dN_deta(0), shapeBatch.dN_dzeta(0));

            const size_t bytesToCopy = m_numNodes * sizeof(CfemReal);
            for (CfemInt q = 1; q < nQp; ++q) {
                std::memcpy(shapeBatch.dN_dxi(q),   shapeBatch.dN_dxi(0),   bytesToCopy);
                std::memcpy(shapeBatch.dN_deta(q),  shapeBatch.dN_deta(0),  bytesToCopy);
                std::memcpy(shapeBatch.dN_dzeta(q), shapeBatch.dN_dzeta(0), bytesToCopy);
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override
        {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            if (nQp == 0) return;

            reference_elem::kernels::Q0HexaKernel::hess(shapeBatch.d2N_dxi2(0), shapeBatch.d2N_deta2(0), shapeBatch.d2N_dzeta2(0),
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

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override
        {
            return cfem::reference_elem::utils::projectFacetRefCoordToHexaRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override
        {
            return cfem::reference_elem::utils::hex1_face_nodes[iFace];
        }
    };

    /**
     * @brief Trilinear Reference Hexahedron (8 nodes, Q1)
     */
    class Q1ReferenceHexahedron : public ReferenceElement
    {
    public:
        Q1ReferenceHexahedron()
        {
            m_numNodes = 8;
            m_dimension = 3;
            m_nodes.resize(m_numNodes);
            // VTK order : Face Z=-1 (0,1,2,3), puis Face Z=1 (4,5,6,7)
            m_nodes[0] = {-1.0, -1.0, -1.0};
            m_nodes[1] = { 1.0, -1.0, -1.0};
            m_nodes[2] = { 1.0,  1.0, -1.0};
            m_nodes[3] = {-1.0,  1.0, -1.0};
            m_nodes[4] = {-1.0, -1.0,  1.0};
            m_nodes[5] = { 1.0, -1.0,  1.0};
            m_nodes[6] = { 1.0,  1.0,  1.0};
            m_nodes[7] = {-1.0,  1.0,  1.0};
        }

        // =========================================================================
        // SCALAR INTERFACES
        // =========================================================================
        
        void evaluateShapesValues(const Vector3D &xi, ShapeInfo& info) const override
        {
            reference_elem::kernels::Q1HexaKernel::values(xi.x, xi.y, xi.z, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D &xi, ShapeInfo& info) const override
        {
            reference_elem::kernels::Q1HexaKernel::grads(xi.x, xi.y, xi.z, 
                                                         info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xi, ShapeInfo& info) const override
        {
            reference_elem::kernels::Q1HexaKernel::hess(xi.x, xi.y, xi.z, 
                                                        info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(),
                                                        info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        // =====================================================================
        // HPC BATCH INTERFACES
        // =====================================================================
        
        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override
        {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q)
            {
                reference_elem::kernels::Q1HexaKernel::values(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                              shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override
        {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q)
            {
                reference_elem::kernels::Q1HexaKernel::grads(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                             shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override
        {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q)
            {
                reference_elem::kernels::Q1HexaKernel::hess(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                            shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                            shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override
        {
            return cfem::reference_elem::utils::projectFacetRefCoordToHexaRefCoord(localFacetId, xiFacet); 
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override
        {
            return cfem::reference_elem::utils::hex8_face_nodes[iFace];
        }
    };

    /**
     * @brief Q1-Mini Reference Hexahedron (9 nodes: 8 vertices + 1 bubble)
     */
    class MiniReferenceHexahedron : public ReferenceElement
    {
    public:
        MiniReferenceHexahedron()
        {
            m_numNodes = 9; // 8 sommets + 1 bulle
            m_dimension = 3;
            m_nodes.resize(m_numNodes);

            // VTK order Q1
            m_nodes[0] = {-1.0, -1.0, -1.0};
            m_nodes[1] = { 1.0, -1.0, -1.0};
            m_nodes[2] = { 1.0,  1.0, -1.0};
            m_nodes[3] = {-1.0,  1.0, -1.0};
            m_nodes[4] = {-1.0, -1.0,  1.0};
            m_nodes[5] = { 1.0, -1.0,  1.0};
            m_nodes[6] = { 1.0,  1.0,  1.0};
            m_nodes[7] = {-1.0,  1.0,  1.0};
            // Bubble node (center)
            m_nodes[8] = { 0.0,  0.0,  0.0};
        }

        // =========================================================================
        // SCALAR INTERFACES
        // =========================================================================

        void evaluateShapesValues(const Vector3D &xi, ShapeInfo& info) const override
        {
            reference_elem::kernels::MiniHexaKernel::values(xi.x, xi.y, xi.z, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D &xi, ShapeInfo& info) const override
        {
            reference_elem::kernels::MiniHexaKernel::grads(xi.x, xi.y, xi.z, 
                                                           info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xi, ShapeInfo& info) const override
        {
            reference_elem::kernels::MiniHexaKernel::hess(xi.x, xi.y, xi.z, 
                                                          info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(),
                                                          info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        // =====================================================================
        // HPC BATCH INTERFACES
        // =====================================================================

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override
        {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::MiniHexaKernel::values(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                                shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override
        {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::MiniHexaKernel::grads(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                               shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override
        {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            // Les dérivées dépendent de la position, pas de memcpy
            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::MiniHexaKernel::hess(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                              shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                              shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override
        {
            return cfem::reference_elem::utils::projectFacetRefCoordToHexaRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override
        {
            return cfem::reference_elem::utils::hex8_face_nodes[iFace]; // Les faces n'incluent pas la bulle (noeud 8)
        }
    };

    /**
     * @brief Tri-Quadratic Reference Hexahedron (27 nodes, Q2)
     */
    class Q2ReferenceHexahedron : public ReferenceElement
    {
    public:
        Q2ReferenceHexahedron()
        {
            m_numNodes = 27;
            m_dimension = 3;
            m_nodes.resize(m_numNodes);
            
            // Define 1D positions r, s, t : -1.0, 0.0, 1.0
            CfemReal pos[3] = {-1.0, 0.0, 1.0};

            for (CfemInt i = 0; i < 3; ++i) {
                for (CfemInt j = 0; j < 3; ++j) {
                    for (CfemInt k = 0; k < 3; ++k) {
                        m_nodes[reference_elem::kernels::Q2HexaKernel::node_map[i][j][k]] = {pos[i], pos[j], pos[k]};
                    }
                }
            }
        }

        // =========================================================================
        // SCALAR INTERFACES
        // =========================================================================

        void evaluateShapesValues(const Vector3D &xi, ShapeInfo& info) const override
        {
            reference_elem::kernels::Q2HexaKernel::values(xi.x, xi.y, xi.z, info.values());
        }

        void evaluateShapesDerivatives(const Vector3D &xi, ShapeInfo& info) const override
        {
            reference_elem::kernels::Q2HexaKernel::grads(xi.x, xi.y, xi.z, 
                                                         info.dN_dxi(), info.dN_deta(), info.dN_dzeta());
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xi, ShapeInfo& info) const override
        {
            reference_elem::kernels::Q2HexaKernel::hess(xi.x, xi.y, xi.z, 
                                                        info.d2N_dxi2(), info.d2N_deta2(), info.d2N_dzeta2(),
                                                        info.d2N_dxideta(), info.d2N_dxidzeta(), info.d2N_detadzeta());
        }

        // =====================================================================
        // HPC BATCH INTERFACES
        // =====================================================================

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override
        {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q)
            {
                reference_elem::kernels::Q2HexaKernel::values(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                              shapeBatch.values(q));
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override
        {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q)
            {
                reference_elem::kernels::Q2HexaKernel::grads(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                             shapeBatch.dN_dxi(q), shapeBatch.dN_deta(q), shapeBatch.dN_dzeta(q));
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints, ShapeInfoBatch& shapeBatch) const override
        {
            const CfemInt nQp = shapeBatch.numQuadraturePoints;
            for (CfemInt q = 0; q < nQp; ++q)
            {
                reference_elem::kernels::Q2HexaKernel::hess(xiPoints[q].x, xiPoints[q].y, xiPoints[q].z, 
                                                            shapeBatch.d2N_dxi2(q), shapeBatch.d2N_deta2(q), shapeBatch.d2N_dzeta2(q),
                                                            shapeBatch.d2N_dxideta(q), shapeBatch.d2N_dxidzeta(q), shapeBatch.d2N_detadzeta(q));
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override
        {
            return cfem::reference_elem::utils::projectFacetRefCoordToHexaRefCoord(localFacetId, xiFacet); 
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override
        {
            return cfem::reference_elem::utils::hex27_face_nodes[iFace];
        }
    };

}

#endif