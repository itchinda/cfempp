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

    static const std::vector<std::vector<int>> hex8_face_nodes = {
        {0, 3, 2, 1}, // Face 0 (Bottom)
        {4, 5, 6, 7}, // Face 1 (Top)
        {0, 1, 5, 4}, // Face 2 (Front)
        {1, 2, 6, 5}, // Face 3 (Right)
        {2, 3, 7, 6}, // Face 4 (Back)
        {3, 0, 4, 7}  // Face 5 (Left)
    };

    static const std::vector<std::vector<int>> hex27_face_nodes = {
        {0, 3, 2, 1, 11, 10, 9, 8, 24},   // Face 0
        {4, 5, 6, 7, 12, 13, 14, 15, 25}, // Face 1
        {0, 1, 5, 4, 8, 17, 12, 16, 22},  // Face 2
        {1, 2, 6, 5, 9, 18, 13, 17, 21},  // Face 3
        {2, 3, 7, 6, 10, 19, 14, 18, 23}, // Face 4
        {3, 0, 4, 7, 11, 16, 15, 19, 20}  // Face 5
    };

    static const std::vector<std::vector<int>> hex1_face_nodes = {
        {0}, {0}, {0}, {0}, {0}, {0}};
}

namespace cfem
{

    class Q1ReferenceHexahedron : public ReferenceElement
    {
    public:
        Q1ReferenceHexahedron()
        {
            m_numNodes = 8;
            m_dimension = 3;
            m_nodes.resize(m_numNodes);
            // Ordre VTK : Face Z=-1 (0,1,2,3), puis Face Z=1 (4,5,6,7)
            m_nodes[0] = {-1, -1, -1};
            m_nodes[1] = {1, -1, -1};
            m_nodes[2] = {1, 1, -1};
            m_nodes[3] = {-1, 1, -1};
            m_nodes[4] = {-1, -1, 1};
            m_nodes[5] = {1, -1, 1};
            m_nodes[6] = {1, 1, 1};
            m_nodes[7] = {-1, 1, 1};
        }

        void evaluateShapesValues(const Vector3D &xi, CfemReal *values) const override
        {
            const CfemReal xm = 1.0 - xi.x;
            const CfemReal xp = 1.0 + xi.x;
            const CfemReal ym = 1.0 - xi.y;
            const CfemReal yp = 1.0 + xi.y;
            const CfemReal zm = 1.0 - xi.z;
            const CfemReal zp = 1.0 + xi.z;

            constexpr CfemReal o8 = 0.125; // 1.0 / 8.0

            values[0] = o8 * xm * ym * zm;
            values[1] = o8 * xp * ym * zm;
            values[2] = o8 * xp * yp * zm;
            values[3] = o8 * xm * yp * zm;
            values[4] = o8 * xm * ym * zp;
            values[5] = o8 * xp * ym * zp;
            values[6] = o8 * xp * yp * zp;
            values[7] = o8 * xm * yp * zp;
        }

        void evaluateShapesDerivatives(const Vector3D &xi,
                                 CfemReal *dN_dxi,
                                 CfemReal *dN_deta,
                                 CfemReal *dN_dzeta) const override
        {
            const CfemReal xm = 1.0 - xi.x;
            const CfemReal xp = 1.0 + xi.x;
            const CfemReal ym = 1.0 - xi.y;
            const CfemReal yp = 1.0 + xi.y;
            const CfemReal zm = 1.0 - xi.z;
            const CfemReal zp = 1.0 + xi.z;

            constexpr CfemReal o8 = 0.125;

            // First derivatives with respect to xi (x)
            dN_dxi[0] = -o8 * ym * zm;
            dN_dxi[1] = o8 * ym * zm;
            dN_dxi[2] = o8 * yp * zm;
            dN_dxi[3] = -o8 * yp * zm;
            dN_dxi[4] = -o8 * ym * zp;
            dN_dxi[5] = o8 * ym * zp;
            dN_dxi[6] = o8 * yp * zp;
            dN_dxi[7] = -o8 * yp * zp;

            // First derivatives with respect to eta (y)
            dN_deta[0] = -o8 * xm * zm;
            dN_deta[1] = -o8 * xp * zm;
            dN_deta[2] = o8 * xp * zm;
            dN_deta[3] = o8 * xm * zm;
            dN_deta[4] = -o8 * xm * zp;
            dN_deta[5] = -o8 * xp * zp;
            dN_deta[6] = o8 * xp * zp;
            dN_deta[7] = o8 * xm * zp;

            // First derivatives with respect to zeta (z)
            dN_dzeta[0] = -o8 * xm * ym;
            dN_dzeta[1] = -o8 * xp * ym;
            dN_dzeta[2] = -o8 * xp * yp;
            dN_dzeta[3] = -o8 * xm * yp;
            dN_dzeta[4] = o8 * xm * ym;
            dN_dzeta[5] = o8 * xp * ym;
            dN_dzeta[6] = o8 * xp * yp;
            dN_dzeta[7] = o8 * xm * yp;
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xi,
                                       CfemReal *d2N_dxi2,
                                       CfemReal *d2N_deta2,
                                       CfemReal *d2N_dzeta2,
                                       CfemReal *d2N_dxideta,
                                       CfemReal *d2N_dxidzeta,
                                       CfemReal *d2N_detadzeta) const override
        {
            const CfemReal xm = 1.0 - xi.x;
            const CfemReal xp = 1.0 + xi.x;
            const CfemReal ym = 1.0 - xi.y;
            const CfemReal yp = 1.0 + xi.y;
            const CfemReal zm = 1.0 - xi.z;
            const CfemReal zp = 1.0 + xi.z;

            constexpr CfemReal o8 = 0.125;

            // Pure second derivatives are strictly zero for trilinear fields
            for (CfemInt i = 0; i < 8; ++i)
            {
                d2N_dxi2[i] = 0.0;
                d2N_deta2[i] = 0.0;
                d2N_dzeta2[i] = 0.0;
            }

            // Cross derivatives: d2N / (dxi * deta)
            d2N_dxideta[0] = o8 * zm;
            d2N_dxideta[1] = -o8 * zm;
            d2N_dxideta[2] = o8 * zm;
            d2N_dxideta[3] = -o8 * zm;
            d2N_dxideta[4] = o8 * zp;
            d2N_dxideta[5] = -o8 * zp;
            d2N_dxideta[6] = o8 * zp;
            d2N_dxideta[7] = -o8 * zp;

            // Cross derivatives: d2N / (dxi * dzeta)
            d2N_dxidzeta[0] = o8 * ym;
            d2N_dxidzeta[1] = -o8 * ym;
            d2N_dxidzeta[2] = -o8 * yp;
            d2N_dxidzeta[3] = o8 * yp;
            d2N_dxidzeta[4] = -o8 * ym;
            d2N_dxidzeta[5] = o8 * ym;
            d2N_dxidzeta[6] = o8 * yp;
            d2N_dxidzeta[7] = -o8 * yp;

            // Cross derivatives: d2N / (deta * dzeta)
            d2N_detadzeta[0] = o8 * xm;
            d2N_detadzeta[1] = o8 * xp;
            d2N_detadzeta[2] = -o8 * xp;
            d2N_detadzeta[3] = -o8 * xm;
            d2N_detadzeta[4] = -o8 * xm;
            d2N_detadzeta[5] = -o8 * xp;
            d2N_detadzeta[6] = o8 * xp;
            d2N_detadzeta[7] = o8 * xm;
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints,
                                 DynamicMatrix<CfemReal> &valuesBatch) const override
        {
            constexpr CfemReal o8 = 0.125;

            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                const CfemReal xm = 1.0 - xiPoints[q].x;
                const CfemReal xp = 1.0 + xiPoints[q].x;
                const CfemReal ym = 1.0 - xiPoints[q].y;
                const CfemReal yp = 1.0 + xiPoints[q].y;
                const CfemReal zm = 1.0 - xiPoints[q].z;
                const CfemReal zp = 1.0 + xiPoints[q].z;

                valuesBatch[q][0] = o8 * xm * ym * zm;
                valuesBatch[q][1] = o8 * xp * ym * zm;
                valuesBatch[q][2] = o8 * xp * yp * zm;
                valuesBatch[q][3] = o8 * xm * yp * zm;
                valuesBatch[q][4] = o8 * xm * ym * zp;
                valuesBatch[q][5] = o8 * xp * ym * zp;
                valuesBatch[q][6] = o8 * xp * yp * zp;
                valuesBatch[q][7] = o8 * xm * yp * zp;
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints,
                                      DynamicMatrix<CfemReal> &dN_dxi,
                                      DynamicMatrix<CfemReal> &dN_deta,
                                      DynamicMatrix<CfemReal> &dN_dzeta) const override
        {
            constexpr CfemReal o8 = 0.125;

            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                const Vector3D &xi = xiPoints[q];
                const CfemReal xm = 1.0 - xi.x;
                const CfemReal xp = 1.0 + xi.x;
                const CfemReal ym = 1.0 - xi.y;
                const CfemReal yp = 1.0 + xi.y;
                const CfemReal zm = 1.0 - xi.z;
                const CfemReal zp = 1.0 + xi.z;

                // d/dxi
                dN_dxi[q][0] = -o8 * ym * zm;
                dN_dxi[q][1] = o8 * ym * zm;
                dN_dxi[q][2] = o8 * yp * zm;
                dN_dxi[q][3] = -o8 * yp * zm;
                dN_dxi[q][4] = -o8 * ym * zp;
                dN_dxi[q][5] = o8 * ym * zp;
                dN_dxi[q][6] = o8 * yp * zp;
                dN_dxi[q][7] = -o8 * yp * zp;

                // d/deta
                dN_deta[q][0] = -o8 * xm * zm;
                dN_deta[q][1] = -o8 * xp * zm;
                dN_deta[q][2] = o8 * xp * zm;
                dN_deta[q][3] = o8 * xm * zm;
                dN_deta[q][4] = -o8 * xm * zp;
                dN_deta[q][5] = -o8 * xp * zp;
                dN_deta[q][6] = o8 * xp * zp;
                dN_deta[q][7] = o8 * xm * zp;

                // d/dzeta
                dN_dzeta[q][0] = -o8 * xm * ym;
                dN_dzeta[q][1] = -o8 * xp * ym;
                dN_dzeta[q][2] = -o8 * xp * yp;
                dN_dzeta[q][3] = -o8 * xm * yp;
                dN_dzeta[q][4] = o8 * xm * ym;
                dN_dzeta[q][5] = o8 * xp * ym;
                dN_dzeta[q][6] = o8 * xp * yp;
                dN_dzeta[q][7] = o8 * xm * yp;
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints,
                                            DynamicMatrix<CfemReal> &d2N_dxi2,
                                            DynamicMatrix<CfemReal> &d2N_deta2,
                                            DynamicMatrix<CfemReal> &d2N_dzeta2,
                                            DynamicMatrix<CfemReal> &d2N_dxideta,
                                            DynamicMatrix<CfemReal> &d2N_dxidzeta,
                                            DynamicMatrix<CfemReal> &d2N_detadzeta) const override
        {
            constexpr CfemReal o8 = 0.125;

            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                const Vector3D &xi = xiPoints[q];
                const CfemReal xm = 1.0 - xi.x;
                const CfemReal xp = 1.0 + xi.x;
                const CfemReal ym = 1.0 - xi.y;
                const CfemReal yp = 1.0 + xi.y;
                const CfemReal zm = 1.0 - xi.z;
                const CfemReal zp = 1.0 + xi.z;

                for (CfemInt i = 0; i < 8; ++i)
                {
                    d2N_dxi2[q][i] = 0.0;
                    d2N_deta2[q][i] = 0.0;
                    d2N_dzeta2[q][i] = 0.0;
                }

                // d2N / (dxi * deta)
                d2N_dxideta[q][0] = o8 * zm;
                d2N_dxideta[q][1] = -o8 * zm;
                d2N_dxideta[q][2] = o8 * zm;
                d2N_dxideta[q][3] = -o8 * zm;
                d2N_dxideta[q][4] = o8 * zp;
                d2N_dxideta[q][5] = -o8 * zp;
                d2N_dxideta[q][6] = o8 * zp;
                d2N_dxideta[q][7] = -o8 * zp;

                // d2N / (dxi * dzeta)
                d2N_dxidzeta[q][0] = o8 * ym;
                d2N_dxidzeta[q][1] = -o8 * ym;
                d2N_dxidzeta[q][2] = -o8 * yp;
                d2N_dxidzeta[q][3] = o8 * yp;
                d2N_dxidzeta[q][4] = -o8 * ym;
                d2N_dxidzeta[q][5] = o8 * ym;
                d2N_dxidzeta[q][6] = o8 * yp;
                d2N_dxidzeta[q][7] = -o8 * yp;

                // d2N / (deta * dzeta)
                d2N_detadzeta[q][0] = o8 * xm;
                d2N_detadzeta[q][1] = o8 * xp;
                d2N_detadzeta[q][2] = -o8 * xp;
                d2N_detadzeta[q][3] = -o8 * xm;
                d2N_detadzeta[q][4] = -o8 * xm;
                d2N_detadzeta[q][5] = -o8 * xp;
                d2N_detadzeta[q][6] = o8 * xp;
                d2N_detadzeta[q][7] = o8 * xm;
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override
        {
            return cfem::reference_elem::utils::projectFacetRefCoordToHexaRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override
        {
            return cfem::reference_elem::utils::hex8_face_nodes[iFace];
        }
    };

    class Q2ReferenceHexahedron : public ReferenceElement
    {
    private:
        // N[0] correspond à -1.0 | N[1] correspond à 0.0 | N[2] correspond à +1.0
        static inline void eval1D_shape(CfemReal u, CfemReal *N)
        {
            const CfemReal u2 = u * u;
            N[0] = 0.5 * (u2 - u);
            N[1] = 1.0 - u2;
            N[2] = 0.5 * (u2 + u);
        }

        static inline void eval1D_deriv(CfemReal u, CfemReal *N, CfemReal *dN)
        {
            const CfemReal u2 = u * u;
            N[0] = 0.5 * (u2 - u);
            dN[0] = u - 0.5;
            N[1] = 1.0 - u2;
            dN[1] = -2.0 * u;
            N[2] = 0.5 * (u2 + u);
            dN[2] = u + 0.5;
        }

        static inline void eval1D_hessian(CfemReal u, CfemReal *N, CfemReal *dN, CfemReal *d2N)
        {
            const CfemReal u2 = u * u;
            N[0] = 0.5 * (u2 - u);
            dN[0] = u - 0.5;
            d2N[0] = 1.0;
            N[1] = 1.0 - u2;
            dN[1] = -2.0 * u;
            d2N[1] = -2.0;
            N[2] = 0.5 * (u2 + u);
            dN[2] = u + 0.5;
            d2N[2] = 1.0;
        }

        // i (X) parcourt : 0 (-1), 1 (0), 2 (+1)
        // j (Y) parcourt : 0 (-1), 1 (0), 2 (+1)
        // k (Z) parcourt : 0 (-1), 1 (0), 2 (+1)
        static constexpr CfemInt node_map[3][3][3] = {
            // i = 0 (X = -1.0)
            {
                {0, 16, 4},   // j=0 (Y=-1) | k=0,1,2 (Z=-1,0,1)
                {11, 20, 15}, // j=1 (Y= 0) | k=0,1,2 (Z=-1,0,1)
                {3, 19, 7}    // j=2 (Y=+1) | k=0,1,2 (Z=-1,0,1)
            },
            // i = 1 (X = 0.0)
            {
                {8, 22, 12},  // j=0 (Y=-1) | k=0,1,2 (Z=-1,0,1)
                {24, 26, 25}, // j=1 (Y= 0) | k=0,1,2 (Z=-1,0,1)
                {10, 23, 14}  // j=2 (Y=+1) | k=0,1,2 (Z=-1,0,1)
            },
            // i = 2 (X = +1.0)
            {
                {1, 17, 5},  // j=0 (Y=-1) | k=0,1,2 (Z=-1,0,1)
                {9, 21, 13}, // j=1 (Y= 0) | k=0,1,2 (Z=-1,0,1)
                {2, 18, 6}   // j=2 (Y=+1) | k=0,1,2 (Z=-1,0,1)
            }};

    public:
        Q2ReferenceHexahedron()
        {
            m_numNodes = 27;
            m_dimension = 3;
            m_nodes.resize(m_numNodes);
            // On définit les positions 1D pour r, s, t : -1.0, 0.0, 1.0
            CfemReal pos[3] = {-1.0, 0.0, 1.0};

            for (CfemInt i = 0; i < 3; ++i)
                for (CfemInt j = 0; j < 3; ++j)
                    for (CfemInt k = 0; k < 3; ++k)
                        m_nodes[node_map[i][j][k]] = {pos[i], pos[j], pos[k]};
        }

        void evaluateShapesValues(const Vector3D &xi,
                            CfemReal *values) const override
        {
            CfemReal Nx[3], Ny[3], Nz[3];
            eval1D_shape(xi.x, Nx);
            eval1D_shape(xi.y, Ny);
            eval1D_shape(xi.z, Nz);

            for (CfemInt i = 0; i < 3; ++i)
            {
                for (CfemInt j = 0; j < 3; ++j)
                {
                    for (CfemInt k = 0; k < 3; ++k)
                    {
                        values[node_map[i][j][k]] = Nx[i] * Ny[j] * Nz[k];
                    }
                }
            }
        }

        void evaluateShapesDerivatives(const Vector3D &xi,
                                 CfemReal *dN_dxi,
                                 CfemReal *dN_deta,
                                 CfemReal *dN_dzeta) const override
        {

            CfemReal Nx[3], Ny[3], Nz[3], dNx[3], dNy[3], dNz[3];
            eval1D_deriv(xi.x, Nx, dNx);
            eval1D_deriv(xi.y, Ny, dNy);
            eval1D_deriv(xi.z, Nz, dNz);

            for (CfemInt i = 0; i < 3; ++i)
            {
                for (CfemInt j = 0; j < 3; ++j)
                {
                    for (CfemInt k = 0; k < 3; ++k)
                    {
                        const CfemInt node = node_map[i][j][k];
                        dN_dxi[node] = dNx[i] * Ny[j] * Nz[k];
                        dN_deta[node] = Nx[i] * dNy[j] * Nz[k];
                        dN_dzeta[node] = Nx[i] * Ny[j] * dNz[k];
                    }
                }
            }
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xi,
                                       CfemReal *d2N_dxi2,
                                       CfemReal *d2N_deta2,
                                       CfemReal *d2N_dzeta2,
                                       CfemReal *d2N_dxideta,
                                       CfemReal *d2N_dxidzeta,
                                       CfemReal *d2N_detadzeta) const override
        {
            CfemReal Nx[3], Ny[3], Nz[3], dNx[3], dNy[3], dNz[3], d2Nx[3], d2Ny[3], d2Nz[3];
            eval1D_hessian(xi.x, Nx, dNx, d2Nx);
            eval1D_hessian(xi.y, Ny, dNy, d2Ny);
            eval1D_hessian(xi.z, Nz, dNz, d2Nz);

            for (CfemInt i = 0; i < 3; ++i)
            {
                for (CfemInt j = 0; j < 3; ++j)
                {
                    for (CfemInt k = 0; k < 3; ++k)
                    {
                        const CfemInt node = node_map[i][j][k];
                        d2N_dxi2[node] = d2Nx[i] * Ny[j] * Nz[k];
                        d2N_deta2[node] = Nx[i] * d2Ny[j] * Nz[k];
                        d2N_dzeta2[node] = Nx[i] * Ny[j] * d2Nz[k];

                        d2N_dxideta[node] = dNx[i] * dNy[j] * Nz[k];
                        d2N_dxidzeta[node] = dNx[i] * Ny[j] * dNz[k];
                        d2N_detadzeta[node] = Nx[i] * dNy[j] * dNz[k];
                    }
                }
            }
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints,
                                 DynamicMatrix<CfemReal> &valuesBatch) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                CfemReal Nx[3], Ny[3], Nz[3];
                eval1D_shape(xiPoints[q].x, Nx);
                eval1D_shape(xiPoints[q].y, Ny);
                eval1D_shape(xiPoints[q].z, Nz);

                for (CfemInt i = 0; i < 3; ++i)
                {
                    for (CfemInt j = 0; j < 3; ++j)
                    {
                        for (CfemInt k = 0; k < 3; ++k)
                        {
                            valuesBatch[q][node_map[i][j][k]] = Nx[i] * Ny[j] * Nz[k];
                        }
                    }
                }
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints,
                                      DynamicMatrix<CfemReal> &dN_dxi,
                                      DynamicMatrix<CfemReal> &dN_deta,
                                      DynamicMatrix<CfemReal> &dN_dzeta) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                CfemReal Nx[3], Ny[3], Nz[3], dNx[3], dNy[3], dNz[3];
                eval1D_deriv(xiPoints[q].x, Nx, dNx);
                eval1D_deriv(xiPoints[q].y, Ny, dNy);
                eval1D_deriv(xiPoints[q].z, Nz, dNz);

                for (CfemInt i = 0; i < 3; ++i)
                {
                    for (CfemInt j = 0; j < 3; ++j)
                    {
                        for (CfemInt k = 0; k < 3; ++k)
                        {
                            const CfemInt node = node_map[i][j][k];
                            dN_dxi[q][node] = dNx[i] * Ny[j] * Nz[k];
                            dN_deta[q][node] = Nx[i] * dNy[j] * Nz[k];
                            dN_dzeta[q][node] = Nx[i] * Ny[j] * dNz[k];
                        }
                    }
                }
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints,
                                            DynamicMatrix<CfemReal> &d2N_dxi2,
                                            DynamicMatrix<CfemReal> &d2N_deta2,
                                            DynamicMatrix<CfemReal> &d2N_dzeta2,
                                            DynamicMatrix<CfemReal> &d2N_dxideta,
                                            DynamicMatrix<CfemReal> &d2N_dxidzeta,
                                            DynamicMatrix<CfemReal> &d2N_detadzeta) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                CfemReal Nx[3], Ny[3], Nz[3], dNx[3], dNy[3], dNz[3], d2Nx[3], d2Ny[3], d2Nz[3];
                eval1D_hessian(xiPoints[q].x, Nx, dNx, d2Nx);
                eval1D_hessian(xiPoints[q].y, Ny, dNy, d2Ny);
                eval1D_hessian(xiPoints[q].z, Nz, dNz, d2Nz);

                for (CfemInt i = 0; i < 3; ++i)
                {
                    for (CfemInt j = 0; j < 3; ++j)
                    {
                        for (CfemInt k = 0; k < 3; ++k)
                        {
                            const CfemInt node = node_map[i][j][k];
                            d2N_dxi2[q][node] = d2Nx[i] * Ny[j] * Nz[k];
                            d2N_deta2[q][node] = Nx[i] * d2Ny[j] * Nz[k];
                            d2N_dzeta2[q][node] = Nx[i] * Ny[j] * d2Nz[k];

                            d2N_dxideta[q][node] = dNx[i] * dNy[j] * Nz[k];
                            d2N_dxidzeta[q][node] = dNx[i] * Ny[j] * dNz[k];
                            d2N_detadzeta[q][node] = Nx[i] * dNy[j] * dNz[k];
                        }
                    }
                }
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override
        {
            return cfem::reference_elem::utils::projectFacetRefCoordToHexaRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override
        {
            return cfem::reference_elem::utils::hex27_face_nodes[iFace];
        }
    };

    class MiniReferenceHexahedron : public ReferenceElement
    {
    protected:
        static constexpr CfemReal xr[8] = {-1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0, -1.0};
        static constexpr CfemReal ys[8] = {-1.0, -1.0, 1.0, 1.0, -1.0, -1.0, 1.0, 1.0};
        static constexpr CfemReal zt[8] = {-1.0, -1.0, -1.0, -1.0, 1.0, 1.0, 1.0, 1.0};
        static constexpr CfemReal o8 = 0.125;

    public:
        MiniReferenceHexahedron()
        {
            m_numNodes = 9; // 8 sommets + 1 bulle
            m_dimension = 3;
            m_nodes.resize(m_numNodes);

            // Ordre VTK Q1
            m_nodes[0] = {-1, -1, -1};
            m_nodes[1] = {1, -1, -1};
            m_nodes[2] = {1, 1, -1};
            m_nodes[3] = {-1, 1, -1};
            m_nodes[4] = {-1, -1, 1};
            m_nodes[5] = {1, -1, 1};
            m_nodes[6] = {1, 1, 1};
            m_nodes[7] = {-1, 1, 1};
            // Nœud Bulle au centre
            m_nodes[8] = {0, 0, 0};
        }

        void evaluateShapesValues(const Vector3D &xi, CfemReal *values) const override
        {
            const CfemReal x = xi.x, y = xi.y, z = xi.z;

            // Calcul de la bulle
            const CfemReal b = (1.0 - x * x) * (1.0 - y * y) * (1.0 - z * z);
            values[8] = b;

            // Calcul des sommets Q1 avec correction
            const CfemReal correction = b * o8;
            for (CfemInt i = 0; i < 8; ++i)
            {
                const CfemReal px = 1.0 + xr[i] * x;
                const CfemReal py = 1.0 + ys[i] * y;
                const CfemReal pz = 1.0 + zt[i] * z;
                values[i] = (o8 * px * py * pz) - correction;
            }
        }

        void evaluateShapesDerivatives(const Vector3D &xi,
                                 CfemReal *dN_dxi,
                                 CfemReal *dN_deta,
                                 CfemReal *dN_dzeta) const override
        {
            const CfemReal x = xi.x, y = xi.y, z = xi.z;
            const CfemReal mx2 = 1.0 - x * x;
            const CfemReal my2 = 1.0 - y * y;
            const CfemReal mz2 = 1.0 - z * z;

            // Gradients de la bulle
            const CfemReal db_dx = -2.0 * x * my2 * mz2;
            const CfemReal db_dy = -2.0 * y * mx2 * mz2;
            const CfemReal db_dz = -2.0 * z * mx2 * my2;

            dN_dxi[8] = db_dx;
            dN_deta[8] = db_dy;
            dN_dzeta[8] = db_dz;

            // Gradients Q1 corrigés
            const CfemReal corr_dx = db_dx * o8;
            const CfemReal corr_dy = db_dy * o8;
            const CfemReal corr_dz = db_dz * o8;

            for (CfemInt i = 0; i < 8; ++i)
            {
                const CfemReal px = 1.0 + xr[i] * x;
                const CfemReal py = 1.0 + ys[i] * y;
                const CfemReal pz = 1.0 + zt[i] * z;

                dN_dxi[i] = (o8 * xr[i] * py * pz) - corr_dx;
                dN_deta[i] = (o8 * ys[i] * px * pz) - corr_dy;
                dN_dzeta[i] = (o8 * zt[i] * px * py) - corr_dz;
            }
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xi,
                                       CfemReal *d2N_dxi2, CfemReal *d2N_deta2, CfemReal *d2N_dzeta2,
                                       CfemReal *d2N_dxideta, CfemReal *d2N_dxidzeta, CfemReal *d2N_detadzeta) const override
        {
            const CfemReal x = xi.x, y = xi.y, z = xi.z;
            const CfemReal mx2 = 1.0 - x * x;
            const CfemReal my2 = 1.0 - y * y;
            const CfemReal mz2 = 1.0 - z * z;

            // Hessiens de la bulle
            const CfemReal d2b_dx2 = -2.0 * my2 * mz2;
            const CfemReal d2b_dy2 = -2.0 * mx2 * mz2;
            const CfemReal d2b_dz2 = -2.0 * mx2 * my2;
            const CfemReal d2b_dxy = 4.0 * x * y * mz2;
            const CfemReal d2b_dxz = 4.0 * x * z * my2;
            const CfemReal d2b_dyz = 4.0 * y * z * mx2;

            d2N_dxi2[8] = d2b_dx2;
            d2N_deta2[8] = d2b_dy2;
            d2N_dzeta2[8] = d2b_dz2;
            d2N_dxideta[8] = d2b_dxy;
            d2N_dxidzeta[8] = d2b_dxz;
            d2N_detadzeta[8] = d2b_dyz;

            // Hessiens Q1 corrigés
            const CfemReal corr_dx2 = d2b_dx2 * o8;
            const CfemReal corr_dy2 = d2b_dy2 * o8;
            const CfemReal corr_dz2 = d2b_dz2 * o8;
            const CfemReal corr_dxy = d2b_dxy * o8;
            const CfemReal corr_dxz = d2b_dxz * o8;
            const CfemReal corr_dyz = d2b_dyz * o8;

            for (CfemInt i = 0; i < 8; ++i)
            {
                const CfemReal px = 1.0 + xr[i] * x;
                const CfemReal py = 1.0 + ys[i] * y;
                const CfemReal pz = 1.0 + zt[i] * z;

                d2N_dxi2[i] = -corr_dx2; // La composante pure Q1 est nulle
                d2N_deta2[i] = -corr_dy2;
                d2N_dzeta2[i] = -corr_dz2;
                d2N_dxideta[i] = (o8 * xr[i] * ys[i] * pz) - corr_dxy;
                d2N_dxidzeta[i] = (o8 * xr[i] * py * zt[i]) - corr_dxz;
                d2N_detadzeta[i] = (o8 * px * ys[i] * zt[i]) - corr_dyz;
            }
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints,
                                 DynamicMatrix<CfemReal> &valuesBatch) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                const CfemReal x = xiPoints[q].x, y = xiPoints[q].y, z = xiPoints[q].z;
                const CfemReal b = (1.0 - x * x) * (1.0 - y * y) * (1.0 - z * z);

                valuesBatch[q][8] = b;
                const CfemReal correction = b * o8;

                for (CfemInt i = 0; i < 8; ++i)
                {
                    const CfemReal px = 1.0 + xr[i] * x;
                    const CfemReal py = 1.0 + ys[i] * y;
                    const CfemReal pz = 1.0 + zt[i] * z;
                    valuesBatch[q][i] = (o8 * px * py * pz) - correction;
                }
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints,
                                      DynamicMatrix<CfemReal> &dN_dxi,
                                      DynamicMatrix<CfemReal> &dN_deta,
                                      DynamicMatrix<CfemReal> &dN_dzeta) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                const CfemReal x = xiPoints[q].x, y = xiPoints[q].y, z = xiPoints[q].z;
                const CfemReal mx2 = 1.0 - x * x;
                const CfemReal my2 = 1.0 - y * y;
                const CfemReal mz2 = 1.0 - z * z;

                const CfemReal db_dx = -2.0 * x * my2 * mz2;
                const CfemReal db_dy = -2.0 * y * mx2 * mz2;
                const CfemReal db_dz = -2.0 * z * mx2 * my2;

                dN_dxi[q][8] = db_dx;
                dN_deta[q][8] = db_dy;
                dN_dzeta[q][8] = db_dz;

                const CfemReal corr_dx = db_dx * o8;
                const CfemReal corr_dy = db_dy * o8;
                const CfemReal corr_dz = db_dz * o8;

                for (CfemInt i = 0; i < 8; ++i)
                {
                    const CfemReal px = 1.0 + xr[i] * x;
                    const CfemReal py = 1.0 + ys[i] * y;
                    const CfemReal pz = 1.0 + zt[i] * z;

                    dN_dxi[q][i] = (o8 * xr[i] * py * pz) - corr_dx;
                    dN_deta[q][i] = (o8 * ys[i] * px * pz) - corr_dy;
                    dN_dzeta[q][i] = (o8 * zt[i] * px * py) - corr_dz;
                }
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints,
                                            DynamicMatrix<CfemReal> &d2N_dxi2,
                                            DynamicMatrix<CfemReal> &d2N_deta2,
                                            DynamicMatrix<CfemReal> &d2N_dzeta2,
                                            DynamicMatrix<CfemReal> &d2N_dxideta,
                                            DynamicMatrix<CfemReal> &d2N_dxidzeta,
                                            DynamicMatrix<CfemReal> &d2N_detadzeta) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                const CfemReal x = xiPoints[q].x, y = xiPoints[q].y, z = xiPoints[q].z;
                const CfemReal mx2 = 1.0 - x * x;
                const CfemReal my2 = 1.0 - y * y;
                const CfemReal mz2 = 1.0 - z * z;

                const CfemReal d2b_dx2 = -2.0 * my2 * mz2;
                const CfemReal d2b_dy2 = -2.0 * mx2 * mz2;
                const CfemReal d2b_dz2 = -2.0 * mx2 * my2;
                const CfemReal d2b_dxy = 4.0 * x * y * mz2;
                const CfemReal d2b_dxz = 4.0 * x * z * my2;
                const CfemReal d2b_dyz = 4.0 * y * z * mx2;

                d2N_dxi2[q][8] = d2b_dx2;
                d2N_deta2[q][8] = d2b_dy2;
                d2N_dzeta2[q][8] = d2b_dz2;
                d2N_dxideta[q][8] = d2b_dxy;
                d2N_dxidzeta[q][8] = d2b_dxz;
                d2N_detadzeta[q][8] = d2b_dyz;

                const CfemReal corr_dx2 = d2b_dx2 * o8;
                const CfemReal corr_dy2 = d2b_dy2 * o8;
                const CfemReal corr_dz2 = d2b_dz2 * o8;
                const CfemReal corr_dxy = d2b_dxy * o8;
                const CfemReal corr_dxz = d2b_dxz * o8;
                const CfemReal corr_dyz = d2b_dyz * o8;

                for (CfemInt i = 0; i < 8; ++i)
                {
                    const CfemReal px = 1.0 + xr[i] * x;
                    const CfemReal py = 1.0 + ys[i] * y;
                    const CfemReal pz = 1.0 + zt[i] * z;

                    d2N_dxi2[q][i] = -corr_dx2;
                    d2N_deta2[q][i] = -corr_dy2;
                    d2N_dzeta2[q][i] = -corr_dz2;
                    d2N_dxideta[q][i] = (o8 * xr[i] * ys[i] * pz) - corr_dxy;
                    d2N_dxidzeta[q][i] = (o8 * xr[i] * py * zt[i]) - corr_dxz;
                    d2N_detadzeta[q][i] = (o8 * px * ys[i] * zt[i]) - corr_dyz;
                }
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

    class Q0ReferenceHexahedron : public ReferenceElement
    {
    public:
        Q0ReferenceHexahedron()
        {
            m_numNodes = 1;
            m_dimension = 3;
            // Un seul nœud au centre du cube de référence [-1, 1]^3
            m_nodes = {Vector3D{0.0, 0.0, 0.0}};
        }

        void evaluateShapesValues(const Vector3D &, CfemReal *values) const override
        {
            values[0] = 1.0;
        }

        void evaluateShapesDerivatives(const Vector3D &,
                                 CfemReal *dN_dxi,
                                 CfemReal *dN_deta,
                                 CfemReal *dN_dzeta) const override
        {
            dN_dxi[0] = 0.0;
            dN_deta[0] = 0.0;
            dN_dzeta[0] = 0.0;
        }

        void evaluateShapesSecondDerivatives(const Vector3D &,
                                       CfemReal *d2N_dxi2,
                                       CfemReal *d2N_deta2,
                                       CfemReal *d2N_dzeta2,
                                       CfemReal *d2N_dxideta,
                                       CfemReal *d2N_dxidzeta,
                                       CfemReal *d2N_detadzeta) const override
        {
            d2N_dxi2[0] = 0.0;
            d2N_deta2[0] = 0.0;
            d2N_dzeta2[0] = 0.0;
            d2N_dxideta[0] = 0.0;
            d2N_dxidzeta[0] = 0.0;
            d2N_detadzeta[0] = 0.0;
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, DynamicMatrix<CfemReal> &values) const override
        {
            // Constant piecewise field value
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                values[q][0] = 1.0;
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints,
                                      DynamicMatrix<CfemReal> &dN_dxi,
                                      DynamicMatrix<CfemReal> &dN_deta,
                                      DynamicMatrix<CfemReal> &dN_dzeta) const override
        {
            // Constant function results in flat zero gradients everywhere
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                dN_dxi[q][0] = 0.0;
                dN_deta[q][0] = 0.0;
                dN_dzeta[q][0] = 0.0;
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints,
                                            DynamicMatrix<CfemReal> &d2N_dxi2,
                                            DynamicMatrix<CfemReal> &d2N_deta2,
                                            DynamicMatrix<CfemReal> &d2N_dzeta2,
                                            DynamicMatrix<CfemReal> &d2N_dxideta,
                                            DynamicMatrix<CfemReal> &d2N_dxidzeta,
                                            DynamicMatrix<CfemReal> &d2N_detadzeta) const override
        {
            // Flat null second derivatives
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                d2N_dxi2[q][0] = 0.0;
                d2N_deta2[q][0] = 0.0;
                d2N_dzeta2[q][0] = 0.0;
                d2N_dxideta[q][0] = 0.0;
                d2N_dxidzeta[q][0] = 0.0;
                d2N_detadzeta[q][0] = 0.0;
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

}

#endif