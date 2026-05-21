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

namespace cfem::reference_elem::utils
{

    inline Vector3D projectFacetRefCoordToQuadRefCoord(CfemInt localFacetId,
                                                       const Vector3D &xiFacet)
    {
        const CfemReal u = xiFacet.x;

        switch (localFacetId)
        {
        case 0:
            return Vector3D(u, -1.0, 0.0); // Face {0,1} : y=-1
        case 1:
            return Vector3D(1.0, u, 0.0); // Face {1,2} : x=1
        case 2:
            return Vector3D(-u, 1.0, 0.0); // Face {2,3} : y=1
        case 3:
            return Vector3D(-1.0, -u, 0.0); // Face {3,0} : x=-1
        default:
            CFEM_ERROR << "Invalid facet ID for Quad cell";
        }
        return Vector3D{}; // Make the compiler happy
    }

    static const std::vector<std::vector<int>> quad4_face_nodes = {
        {0, 1}, {1, 2}, {2, 3}, {3, 0}};
    static const std::vector<std::vector<int>> quad9_face_nodes = {
        {0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7}};
    static const std::vector<std::vector<int>> quad1_face_nodes = {
        {0}, {0}, {0}, {0}};
}

namespace cfem
{

    class Q1ReferenceQuad : public ReferenceElement
    {
    public:
        Q1ReferenceQuad()
        {
            m_numNodes = 4;
            m_dimension = 2;
            m_nodes.resize(m_numNodes);

            m_nodes[0] = Vector3D{-1.0, -1.0, 0.0};
            m_nodes[1] = Vector3D{1.0, -1.0, 0.0};
            m_nodes[2] = Vector3D{1.0, 1.0, 0.0};
            m_nodes[3] = Vector3D{-1.0, 1.0, 0.0};
        }

        void evaluateShapesValues(const Vector3D &xi,
                            CfemReal *values) const override
        {
            constexpr CfemReal o4 = 0.25;

            const CfemReal xm = 1.0 - xi.x;
            const CfemReal xp = 1.0 + xi.x;
            const CfemReal ym = 1.0 - xi.y;
            const CfemReal yp = 1.0 + xi.y;

            values[0] = o4 * xm * ym;
            values[1] = o4 * xp * ym;
            values[2] = o4 * xp * yp;
            values[3] = o4 * xm * yp;
        }

        void evaluateShapesDerivatives(const Vector3D &xi,
                                 CfemReal *dN_dxi,
                                 CfemReal *dN_deta,
                                 CfemReal *dN_dzeta) const override
        {
            constexpr CfemReal o4 = 0.25;

            const CfemReal xm = 1.0 - xi.x;
            const CfemReal xp = 1.0 + xi.x;
            const CfemReal ym = 1.0 - xi.y;
            const CfemReal yp = 1.0 + xi.y;

            // Dérivées par rapport à x (xi)
            dN_dxi[0] = -o4 * ym;
            dN_dxi[1] = o4 * ym;
            dN_dxi[2] = o4 * yp;
            dN_dxi[3] = -o4 * yp;

            // Dérivées par rapport à y (eta)
            dN_deta[0] = -o4 * xm;
            dN_deta[1] = -o4 * xp;
            dN_deta[2] = o4 * xp;
            dN_deta[3] = o4 * xm;

            // 2D : Z est nul
            for (CfemInt i = 0; i < 4; ++i)
                dN_dzeta[i] = 0.0;
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xi,
                                       CfemReal *d2N_dxi2,
                                       CfemReal *d2N_deta2,
                                       CfemReal *d2N_dzeta2,
                                       CfemReal *d2N_dxideta,
                                       CfemReal *d2N_dxidzeta,
                                       CfemReal *d2N_detadzeta) const override
        {
            constexpr CfemReal o4 = 0.25;

            for (CfemInt i = 0; i < 4; ++i)
            {
                d2N_dxi2[i] = 0.0;
                d2N_deta2[i] = 0.0;
                d2N_dzeta2[i] = 0.0;
                d2N_dxidzeta[i] = 0.0;
                d2N_detadzeta[i] = 0.0;
            }

            // Seule la dérivée croisée xy survit dans un bilinéaire
            d2N_dxideta[0] = o4;
            d2N_dxideta[1] = -o4;
            d2N_dxideta[2] = o4;
            d2N_dxideta[3] = -o4;
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints,
                                 DynamicMatrix<CfemReal> &valuesBatch) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            constexpr CfemReal o4 = 0.25;

            for (CfemInt q = 0; q < nQp; ++q)
            {
                const CfemReal xm = 1.0 - xiPoints[q].x;
                const CfemReal xp = 1.0 + xiPoints[q].x;
                const CfemReal ym = 1.0 - xiPoints[q].y;
                const CfemReal yp = 1.0 + xiPoints[q].y;

                valuesBatch[q][0] = o4 * xm * ym;
                valuesBatch[q][1] = o4 * xp * ym;
                valuesBatch[q][2] = o4 * xp * yp;
                valuesBatch[q][3] = o4 * xm * yp;
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> &xiPoints,
                                      DynamicMatrix<CfemReal> &dN_dxi,
                                      DynamicMatrix<CfemReal> &dN_deta,
                                      DynamicMatrix<CfemReal> &dN_dzeta) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            constexpr CfemReal o4 = 0.25;

            for (CfemInt q = 0; q < nQp; ++q)
            {
                const CfemReal xm = 1.0 - xiPoints[q].x;
                const CfemReal xp = 1.0 + xiPoints[q].x;
                const CfemReal ym = 1.0 - xiPoints[q].y;
                const CfemReal yp = 1.0 + xiPoints[q].y;

                // Dérivées par rapport à x (xi)
                dN_dxi[q][0] = -o4 * ym;
                dN_dxi[q][1] = o4 * ym;
                dN_dxi[q][2] = o4 * yp;
                dN_dxi[q][3] = -o4 * yp;

                // Dérivées par rapport à y (eta)
                dN_deta[q][0] = -o4 * xm;
                dN_deta[q][1] = -o4 * xp;
                dN_deta[q][2] = o4 * xp;
                dN_deta[q][3] = o4 * xm;

                // 2D : Z est nul
                for (CfemInt i = 0; i < 4; ++i)
                    dN_dzeta[q][i] = 0.0;
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
            constexpr CfemReal o4 = 0.25;

            for (CfemInt q = 0; q < nQp; ++q)
            {
                for (CfemInt i = 0; i < 4; ++i)
                {
                    d2N_dxi2[q][i] = 0.0;
                    d2N_deta2[q][i] = 0.0;
                    d2N_dzeta2[q][i] = 0.0;
                    d2N_dxidzeta[q][i] = 0.0;
                    d2N_detadzeta[q][i] = 0.0;
                }

                // Seule la dérivée croisée xy survit dans un bilinéaire
                d2N_dxideta[q][0] = o4;
                d2N_dxideta[q][1] = -o4;
                d2N_dxideta[q][2] = o4;
                d2N_dxideta[q][3] = -o4;
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override
        {
            return cfem::reference_elem::utils::projectFacetRefCoordToQuadRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override
        {
            return cfem::reference_elem::utils::quad4_face_nodes[iFace];
        }
    };

    class Q2ReferenceQuad : public ReferenceElement
    {
    protected:
        // Helpers 1D P2 (-1.0, 0.0, 1.0)
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

        // Connectivité Q2 stricte (X, Y)
        // i (X) parcourt: 0 (-1), 1 (0), 2 (1)
        // j (Y) parcourt: 0 (-1), 1 (0), 2 (1)
        static constexpr CfemInt node_map[3][3] = {
            {0, 7, 3}, // X = -1 | Y = -1, 0, 1
            {4, 8, 6}, // X =  0 | Y = -1, 0, 1
            {1, 5, 2}  // X = +1 | Y = -1, 0, 1
        };

    public:
        Q2ReferenceQuad()
        {
            m_numNodes = 9;
            m_dimension = 2;
            m_nodes.resize(m_numNodes);
            // Coordonnées de référence sur [-1, 1]
            // Sommets
            m_nodes[0] = Vector3D{-1.0, -1.0, 0.0}; // i=0, j=0
            m_nodes[1] = Vector3D{1.0, -1.0, 0.0};  // i=2, j=0
            m_nodes[2] = Vector3D{1.0, 1.0, 0.0};   // i=2, j=2
            m_nodes[3] = Vector3D{-1.0, 1.0, 0.0};  // i=0, j=2

            // Milieux d'arêtes
            m_nodes[4] = Vector3D{0.0, -1.0, 0.0}; // i=1, j=0
            m_nodes[5] = Vector3D{1.0, 0.0, 0.0};  // i=2, j=1
            m_nodes[6] = Vector3D{0.0, 1.0, 0.0};  // i=1, j=2
            m_nodes[7] = Vector3D{-1.0, 0.0, 0.0}; // i=0, j=1

            // Centre
            m_nodes[8] = Vector3D{0.0, 0.0, 0.0}; // i=1, j=1
        }

        void evaluateShapesValues(const Vector3D &xi,
                            CfemReal *values) const override
        {
            CfemReal Nx[3], Ny[3], dNx[3], dNy[3], d2Nx[3], d2Ny[3];
            eval1D_hessian(xi.x, Nx, dNx, d2Nx);
            eval1D_hessian(xi.y, Ny, dNy, d2Ny);

            for (CfemInt i = 0; i < 3; ++i)
            {
                for (CfemInt j = 0; j < 3; ++j)
                {
                    values[node_map[i][j]] = Nx[i] * Ny[j];
                }
            }
        }

        void evaluateShapesDerivatives(const Vector3D &xi,
                                 CfemReal *dN_dxi,
                                 CfemReal *dN_deta,
                                 CfemReal *dN_dzeta) const override
        {
            CfemReal Nx[3], Ny[3], dNx[3], dNy[3], d2Nx[3], d2Ny[3];
            eval1D_hessian(xi.x, Nx, dNx, d2Nx);
            eval1D_hessian(xi.y, Ny, dNy, d2Ny);

            for (CfemInt i = 0; i < 3; ++i)
            {
                for (CfemInt j = 0; j < 3; ++j)
                {
                    const CfemInt node = node_map[i][j];
                    dN_dxi[node] = dNx[i] * Ny[j];
                    dN_deta[node] = Nx[i] * dNy[j];
                    dN_dzeta[node] = 0.0;
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

            CfemReal Nx[3], Ny[3], dNx[3], dNy[3], d2Nx[3], d2Ny[3];
            eval1D_hessian(xi.x, Nx, dNx, d2Nx);
            eval1D_hessian(xi.y, Ny, dNy, d2Ny);

            for (CfemInt i = 0; i < 3; ++i)
            {
                for (CfemInt j = 0; j < 3; ++j)
                {
                    const CfemInt node = node_map[i][j];
                    d2N_dxi2[node] = d2Nx[i] * Ny[j];
                    d2N_deta2[node] = Nx[i] * d2Ny[j];
                    d2N_dzeta2[node] = 0.0;

                    d2N_dxideta[node] = dNx[i] * dNy[j];
                    d2N_dxidzeta[node] = 0.0;
                    d2N_detadzeta[node] = 0.0;
                }
            }
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints,
                                 DynamicMatrix<CfemReal> &valuesBatch) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                CfemReal Nx[3], Ny[3], dNx[3], dNy[3], d2Nx[3], d2Ny[3];
                eval1D_hessian(xiPoints[q].x, Nx, dNx, d2Nx);
                eval1D_hessian(xiPoints[q].y, Ny, dNy, d2Ny);

                for (CfemInt i = 0; i < 3; ++i)
                {
                    for (CfemInt j = 0; j < 3; ++j)
                    {
                        valuesBatch[q][node_map[i][j]] = Nx[i] * Ny[j];
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
                CfemReal Nx[3], Ny[3], dNx[3], dNy[3], d2Nx[3], d2Ny[3];
                eval1D_hessian(xiPoints[q].x, Nx, dNx, d2Nx);
                eval1D_hessian(xiPoints[q].y, Ny, dNy, d2Ny);

                for (CfemInt i = 0; i < 3; ++i)
                {
                    for (CfemInt j = 0; j < 3; ++j)
                    {
                        const CfemInt node = node_map[i][j];
                        dN_dxi[q][node] = dNx[i] * Ny[j];
                        dN_deta[q][node] = Nx[i] * dNy[j];
                        dN_dzeta[q][node] = 0.0;
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
                CfemReal Nx[3], Ny[3], dNx[3], dNy[3], d2Nx[3], d2Ny[3];
                eval1D_hessian(xiPoints[q].x, Nx, dNx, d2Nx);
                eval1D_hessian(xiPoints[q].y, Ny, dNy, d2Ny);

                for (CfemInt i = 0; i < 3; ++i)
                {
                    for (CfemInt j = 0; j < 3; ++j)
                    {
                        const CfemInt node = node_map[i][j];
                        d2N_dxi2[q][node] = d2Nx[i] * Ny[j];
                        d2N_deta2[q][node] = Nx[i] * d2Ny[j];
                        d2N_dzeta2[q][node] = 0.0;

                        d2N_dxideta[q][node] = dNx[i] * dNy[j];
                        d2N_dxidzeta[q][node] = 0.0;
                        d2N_detadzeta[q][node] = 0.0;
                    }
                }
            }
        }
        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override
        {
            return cfem::reference_elem::utils::projectFacetRefCoordToQuadRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override
        {
            return cfem::reference_elem::utils::quad9_face_nodes[iFace];
        }
    };

    class MiniReferenceQuad : public ReferenceElement
    {
    public:
        MiniReferenceQuad()
        {
            m_numNodes = 5;
            m_dimension = 2;
            m_nodes.resize(5);
            m_nodes[0] = {-1, -1, 0};
            m_nodes[1] = {1, -1, 0};
            m_nodes[2] = {1, 1, 0};
            m_nodes[3] = {-1, 1, 0};
            m_nodes[4] = {0, 0, 0}; // Bulle
        }

        void evaluateShapesValues(const Vector3D& xi, 
                            CfemReal* values) const override 
        {
            constexpr CfemReal o4 = 0.25;

                const CfemReal x = xi.x;
                const CfemReal y = xi.y;

                // Bulle Centrale
                const CfemReal b = (1.0 - x*x) * (1.0 - y*y);
                values[4] = b;

                // Sommets
                const CfemReal correction = b * o4;
                values[0] = (o4 * (1.0 - x) * (1.0 - y)) - correction;
                values[1] = (o4 * (1.0 + x) * (1.0 - y)) - correction;
                values[2] = (o4 * (1.0 + x) * (1.0 + y)) - correction;
                values[3] = (o4 * (1.0 - x) * (1.0 + y)) - correction;
        }

        void evaluateShapesDerivatives(const Vector3D& xi, 
                                      CfemReal* dN_dxi, 
                                      CfemReal* dN_deta, 
                                      CfemReal* dN_dzeta) const override 
        {
            constexpr CfemReal o4 = 0.25;

                const CfemReal x = xi.x;
                const CfemReal y = xi.y;

                const CfemReal db_dx = -2.0 * x * (1.0 - y*y);
                const CfemReal db_dy = -2.0 * y * (1.0 - x*x);

                dN_dxi[4]   = db_dx;
                dN_deta[4]  = db_dy;
                dN_dzeta[4] = 0.0;

                const CfemReal c_dx = db_dx * o4;
                const CfemReal c_dy = db_dy * o4;

                dN_dxi[0] = (-o4 * (1.0 - y)) - c_dx; dN_deta[0] = (-o4 * (1.0 - x)) - c_dy; dN_dzeta[0] = 0.0;
                dN_dxi[1] = ( o4 * (1.0 - y)) - c_dx; dN_deta[1] = (-o4 * (1.0 + x)) - c_dy; dN_dzeta[1] = 0.0;
                dN_dxi[2] = ( o4 * (1.0 + y)) - c_dx; dN_deta[2] = ( o4 * (1.0 + x)) - c_dy; dN_dzeta[2] = 0.0;
                dN_dxi[3] = (-o4 * (1.0 + y)) - c_dx; dN_deta[3] = ( o4 * (1.0 - x)) - c_dy; dN_dzeta[3] = 0.0;
            
        }

        void evaluateShapesSecondDerivatives(const Vector3D& xi,
                                       CfemReal* d2N_dxi2,   
                                       CfemReal* d2N_deta2,  
                                       CfemReal* d2N_dzeta2,
                                       CfemReal* d2N_dxideta,
                                       CfemReal* d2N_dxidzeta,
                                       CfemReal* d2N_detadzeta) const override 
        {

            constexpr CfemReal o4 = 0.25;

                const CfemReal x = xi.x;
                const CfemReal y = xi.y;

                // Hessiens de la bulle
                const CfemReal d2b_dx2 = -2.0 * (1.0 - y*y);
                const CfemReal d2b_dy2 = -2.0 * (1.0 - x*x);
                const CfemReal d2b_dxy =  4.0 * x * y;

                d2N_dxi2[4]      = d2b_dx2;
                d2N_deta2[4]     = d2b_dy2;
                d2N_dxideta[4]   = d2b_dxy;
                d2N_dzeta2[4]    = 0.0; d2N_dxidzeta[4] = 0.0; d2N_detadzeta[4] = 0.0;

                // Sommets (Correction de la composante purement nulle pour xx et yy du Q1)
                const CfemReal c_dx2 = d2b_dx2 * o4;
                const CfemReal c_dy2 = d2b_dy2 * o4;
                const CfemReal c_dxy = d2b_dxy * o4;

                for(CfemInt i=0; i<4; ++i) {
                    d2N_dxi2[i] = -c_dx2;
                    d2N_deta2[i] = -c_dy2;
                    d2N_dzeta2[i] = 0.0;
                    d2N_dxidzeta[i] = 0.0;
                    d2N_detadzeta[i] = 0.0;
                }

                d2N_dxideta[0] =  o4 - c_dxy;
                d2N_dxideta[1] = -o4 - c_dxy;
                d2N_dxideta[2] =  o4 - c_dxy;
                d2N_dxideta[3] = -o4 - c_dxy;
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, 
                                DynamicMatrix<CfemReal>& valuesBatch) const override 
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            constexpr CfemReal o4 = 0.25;

            for (CfemInt q = 0; q < nQp; ++q) {
                const CfemReal x = xiPoints[q].x;
                const CfemReal y = xiPoints[q].y;

                // Bulle Centrale
                const CfemReal b = (1.0 - x*x) * (1.0 - y*y);
                valuesBatch[q][4] = b;

                // Sommets
                const CfemReal correction = b * o4;
                valuesBatch[q][0] = (o4 * (1.0 - x) * (1.0 - y)) - correction;
                valuesBatch[q][1] = (o4 * (1.0 + x) * (1.0 - y)) - correction;
                valuesBatch[q][2] = (o4 * (1.0 + x) * (1.0 + y)) - correction;
                valuesBatch[q][3] = (o4 * (1.0 - x) * (1.0 + y)) - correction;
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, 
                                    DynamicMatrix<CfemReal>& dN_dxi, 
                                    DynamicMatrix<CfemReal>& dN_deta, 
                                    DynamicMatrix<CfemReal>& dN_dzeta) const override 
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            constexpr CfemReal o4 = 0.25;

            for (CfemInt q = 0; q < nQp; ++q) {
                const CfemReal x = xiPoints[q].x;
                const CfemReal y = xiPoints[q].y;

                const CfemReal db_dx = -2.0 * x * (1.0 - y*y);
                const CfemReal db_dy = -2.0 * y * (1.0 - x*x);

                dN_dxi[q][4]   = db_dx;
                dN_deta[q][4]  = db_dy;
                dN_dzeta[q][4] = 0.0;

                const CfemReal c_dx = db_dx * o4;
                const CfemReal c_dy = db_dy * o4;

                dN_dxi[q][0] = (-o4 * (1.0 - y)) - c_dx; dN_deta[q][0] = (-o4 * (1.0 - x)) - c_dy; dN_dzeta[q][0] = 0.0;
                dN_dxi[q][1] = ( o4 * (1.0 - y)) - c_dx; dN_deta[q][1] = (-o4 * (1.0 + x)) - c_dy; dN_dzeta[q][1] = 0.0;
                dN_dxi[q][2] = ( o4 * (1.0 + y)) - c_dx; dN_deta[q][2] = ( o4 * (1.0 + x)) - c_dy; dN_dzeta[q][2] = 0.0;
                dN_dxi[q][3] = (-o4 * (1.0 + y)) - c_dx; dN_deta[q][3] = ( o4 * (1.0 - x)) - c_dy; dN_dzeta[q][3] = 0.0;
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D>& xiPoints,
                                            DynamicMatrix<CfemReal>& d2N_dxi2,   
                                            DynamicMatrix<CfemReal>& d2N_deta2,
                                            DynamicMatrix<CfemReal>& d2N_dzeta2,
                                            DynamicMatrix<CfemReal>& d2N_dxideta,
                                            DynamicMatrix<CfemReal>& d2N_dxidzeta,
                                            DynamicMatrix<CfemReal>& d2N_detadzeta) const override 
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            constexpr CfemReal o4 = 0.25;

            for (CfemInt q = 0; q < nQp; ++q) {
                const CfemReal x = xiPoints[q].x;
                const CfemReal y = xiPoints[q].y;

                // Hessiens de la bulle
                const CfemReal d2b_dx2 = -2.0 * (1.0 - y*y);
                const CfemReal d2b_dy2 = -2.0 * (1.0 - x*x);
                const CfemReal d2b_dxy =  4.0 * x * y;

                d2N_dxi2[q][4]      = d2b_dx2;
                d2N_deta2[q][4]     = d2b_dy2;
                d2N_dxideta[q][4]   = d2b_dxy;
                d2N_dzeta2[q][4]    = 0.0; d2N_dxidzeta[q][4] = 0.0; d2N_detadzeta[q][4] = 0.0;

                // Sommets (Correction de la composante purement nulle pour xx et yy du Q1)
                const CfemReal c_dx2 = d2b_dx2 * o4;
                const CfemReal c_dy2 = d2b_dy2 * o4;
                const CfemReal c_dxy = d2b_dxy * o4;

                for(CfemInt i=0; i<4; ++i) {
                    d2N_dxi2[q][i] = -c_dx2;
                    d2N_deta2[q][i] = -c_dy2;
                    d2N_dzeta2[q][i] = 0.0;
                    d2N_dxidzeta[q][i] = 0.0;
                    d2N_detadzeta[q][i] = 0.0;
                }

                d2N_dxideta[q][0] =  o4 - c_dxy;
                d2N_dxideta[q][1] = -o4 - c_dxy;
                d2N_dxideta[q][2] =  o4 - c_dxy;
                d2N_dxideta[q][3] = -o4 - c_dxy;
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override
        {
            return reference_elem::utils::projectFacetRefCoordToQuadRefCoord(localFacetId, xiFacet);
        }
        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override
        {
            return reference_elem::utils::quad4_face_nodes[iFace];
        }
    };

    class Q0ReferenceQuad : public ReferenceElement
    {
    public:
        Q0ReferenceQuad()
        {
            m_numNodes = 1;
            m_dimension = 2;

            // Un seul nœud au centre du carré de référence [-1, 1]x[-1, 1]
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
            return cfem::reference_elem::utils::projectFacetRefCoordToQuadRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override
        {
            return cfem::reference_elem::utils::quad1_face_nodes[iFace];
        }
    };
}
#endif