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

#ifndef CFEM_REFERENCE_TETRAHEDRON
#define CFEM_REFERENCE_TETRAHEDRON

#include "ReferenceElement.h"

namespace cfem::reference_elem::utils{

    inline Vector3D projectFacetRefCoordToTetraRefCoord(CfemInt localFacetId, 
                                                      const Vector3D& xiFacet) 
    {
        const CfemReal u = xiFacet.x;
        const CfemReal v = xiFacet.y;

        CfemReal n0 = 1.0 - u - v;
        switch (localFacetId) {
            case 0: return Vector3D(u, v, 0.0);                // Face {0,1,2} : z=0
            case 1: return Vector3D(u, 0.0, v);                // Face {0,1,3} : y=0
            case 2: return Vector3D(n0, u, v);                 // Face {1,2,3} : x+y+z=1
            case 3: return Vector3D(0.0, n0, v);               // Face {2,0,3} : x=0
            default: CFEM_ERROR << "Invalid facet ID for Tetra";
        }
        return Vector3D{}; // Make the compiler happy
    }

    static const std::vector<std::vector<int>> tet4_face_nodes = {
        {0, 2, 1}, 
        {0, 1, 3}, 
        {1, 2, 3}, 
        {2, 0, 3}
    };

    static const std::vector<std::vector<int>> tet10_face_nodes = {
        {0, 2, 1, 6, 5, 4}, 
        {0, 1, 3, 4, 8, 7}, 
        {1, 2, 3, 5, 9, 8}, 
        {2, 0, 3, 6, 7, 9}
    };

    static const std::vector<std::vector<int>> tet1_face_nodes = {
        {0}, {0}, {0}, {0}
    };
}

namespace cfem {

    /**
     * @brief Tétraèdre Linéaire (4 nœuds)
     */
    class P1ReferenceTetrahedron : public ReferenceElement {
    public:
        P1ReferenceTetrahedron() {
            m_numNodes = 4;
            m_dimension = 3;
            
            m_nodes = {
                Vector3D{0.0, 0.0, 0.0}, 
                Vector3D{1.0, 0.0, 0.0}, 
                Vector3D{0.0, 1.0, 0.0}, 
                Vector3D{0.0, 0.0, 1.0} 
            };
        }

        void evaluateShapesValues(const Vector3D& xi, CfemReal* values) const override 
        {
            values[0] = 1.0 - xi.x - xi.y - xi.z;
            values[1] = xi.x;
            values[2] = xi.y;
            values[3] = xi.z;
        }

        void evaluateShapesDerivatives(const Vector3D&, 
                                 CfemReal* dN_dxi, 
                                 CfemReal* dN_deta, 
                                 CfemReal* dN_dzeta) const override 
        {
            // d/dxi
            dN_dxi[0] = -1.0; dN_dxi[1] =  1.0; dN_dxi[2] =  0.0; dN_dxi[3] =  0.0;
            // d/deta
            dN_deta[0] = -1.0; dN_deta[1] =  0.0; dN_deta[2] =  1.0; dN_deta[3] =  0.0;
            // d/dzeta
            dN_dzeta[0] = -1.0; dN_dzeta[1] =  0.0; dN_dzeta[2] =  0.0; dN_dzeta[3] =  1.0;
        }

        void evaluateShapesSecondDerivatives(const Vector3D&,
                                       CfemReal* d2N_dxi2, 
                                       CfemReal* d2N_deta2, 
                                       CfemReal* d2N_dzeta2,
                                       CfemReal* d2N_dxideta, 
                                       CfemReal* d2N_dxidzeta, 
                                       CfemReal* d2N_detadzeta) const override 
        {
            for (CfemInt i = 0; i < 4; ++i) {
                d2N_dxi2[i]      = 0.0;
                d2N_deta2[i]     = 0.0;
                d2N_dzeta2[i]    = 0.0;
                d2N_dxideta[i]   = 0.0;
                d2N_dxidzeta[i]  = 0.0;
                d2N_detadzeta[i] = 0.0;
            }
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, DynamicMatrix<CfemReal>& values) const override
        {
            // Constant piecewise field value
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                const auto& xi =  xiPoints[q];
                values[q][0] = 1.0 - xi.x - xi.y - xi.z;
                values[q][1] = xi.x;
                values[q][2] = xi.y;
                values[q][3] = xi.z;
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> & xiPoints, 
                                      DynamicMatrix<CfemReal> &dN_dxi, 
                                      DynamicMatrix<CfemReal> &dN_deta, 
                                      DynamicMatrix<CfemReal> &dN_dzeta) const override
        {
            // Constant function results in flat zero gradients everywhere
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                // d/dxi
                dN_dxi[q][0] = -1.0; dN_dxi[q][1] =  1.0; dN_dxi[q][2] =  0.0; dN_dxi[q][3] =  0.0;
                // d/deta
                dN_deta[q][0] = -1.0; dN_deta[q][1] =  0.0; dN_deta[q][2] =  1.0; dN_deta[q][3] =  0.0;
                // d/dzeta
                dN_dzeta[q][0] = -1.0; dN_dzeta[q][1] =  0.0; dN_dzeta[q][2] =  0.0; dN_dzeta[q][3] =  1.0;
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
            for (CfemInt q = 0; q < nQp; ++q) {
                d2N_dxi2[q][0]      = 0.0;
                d2N_deta2[q][0]     = 0.0;
                d2N_dzeta2[q][0]    = 0.0;
                d2N_dxideta[q][0]   = 0.0;
                d2N_dxidzeta[q][0]  = 0.0;
                d2N_detadzeta[q][0] = 0.0;
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToTetraRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tet4_face_nodes[iFace];
        }
    };

    /**
     * @brief Tétraèdre Quadratique (10 nœuds)
     */
    class P2ReferenceTetrahedron : public ReferenceElement {
    public:
        P2ReferenceTetrahedron() {
            m_numNodes = 10;
            m_dimension = 3;
            
            m_nodes = {
                Vector3D{0.0, 0.0 ,0.0},
                Vector3D{1.0, 0.0 ,0.0}, 
                Vector3D{0.0, 1.0 ,0.0}, 
                Vector3D{0.0, 0.0 ,1.0},
                // Milieux d'arêtes (Standard convention: 01, 12, 20, 03, 13, 23)
                Vector3D{0.5, 0.0, 0.0}, 
                Vector3D{0.5, 0.5, 0.0}, 
                Vector3D{0.0, 0.5, 0.0},
                Vector3D{0.0, 0.0, 0.5}, 
                Vector3D{0.5, 0.0, 0.5}, 
                Vector3D{0.0, 0.5, 0.5}
            };
        }

        void evaluateShapesValues(const Vector3D& xi, CfemReal* values) const override 
        {
            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;
            const CfemReal L3 = xi.z;
            const CfemReal L0 = 1.0 - L1 - L2 - L3;

            // Sommets
            values[0] = L0 * (2.0 * L0 - 1.0);
            values[1] = L1 * (2.0 * L1 - 1.0);
            values[2] = L2 * (2.0 * L2 - 1.0);
            values[3] = L3 * (2.0 * L3 - 1.0);

            // Arêtes horizontales (base)
            values[4] = 4.0 * L0 * L1;
            values[5] = 4.0 * L1 * L2;
            values[6] = 4.0 * L2 * L0;

            // Arêtes verticales (vers L3)
            values[7] = 4.0 * L0 * L3;
            values[8] = 4.0 * L1 * L3;
            values[9] = 4.0 * L2 * L3;
        }

        void evaluateShapesDerivatives(const Vector3D& xi, 
                                 CfemReal* dN_dxi, 
                                 CfemReal* dN_deta, 
                                 CfemReal* dN_dzeta) const override 
        {
            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;
            const CfemReal L3 = xi.z;

            // Dérivées par rapport à xi (x)
            dN_dxi[0] = -3.0 + 4.0*L1 + 4.0*L2 + 4.0*L3;
            dN_dxi[1] =  4.0*L1 - 1.0;
            dN_dxi[2] =  0.0;
            dN_dxi[3] =  0.0;
            dN_dxi[4] =  4.0 - 8.0*L1 - 4.0*L2 - 4.0*L3;
            dN_dxi[5] =  4.0*L2;
            dN_dxi[6] = -4.0*L2;
            dN_dxi[7] = -4.0*L3;
            dN_dxi[8] =  4.0*L3;
            dN_dxi[9] =  0.0;

            // Dérivées par rapport à eta (y)
            dN_deta[0] = -3.0 + 4.0*L1 + 4.0*L2 + 4.0*L3;
            dN_deta[1] =  0.0;
            dN_deta[2] =  4.0*L2 - 1.0;
            dN_deta[3] =  0.0;
            dN_deta[4] = -4.0*L1;
            dN_deta[5] =  4.0*L1;
            dN_deta[6] =  4.0 - 4.0*L1 - 8.0*L2 - 4.0*L3;
            dN_deta[7] = -4.0*L3;
            dN_deta[8] =  0.0;
            dN_deta[9] =  4.0*L3;

            // Dérivées par rapport à zeta (z)
            dN_dzeta[0] = -3.0 + 4.0*L1 + 4.0*L2 + 4.0*L3;
            dN_dzeta[1] =  0.0;
            dN_dzeta[2] =  0.0;
            dN_dzeta[3] =  4.0*L3 - 1.0;
            dN_dzeta[4] = -4.0*L1;
            dN_dzeta[5] =  0.0;
            dN_dzeta[6] = -4.0*L2;
            dN_dzeta[7] =  4.0 - 4.0*L1 - 4.0*L2 - 8.0*L3;
            dN_dzeta[8] =  4.0*L1;
            dN_dzeta[9] =  4.0*L2;
        }

        void evaluateShapesSecondDerivatives(const Vector3D&,
                                       CfemReal* d2N_dxi2, 
                                       CfemReal* d2N_deta2, 
                                       CfemReal* d2N_dzeta2,
                                       CfemReal* d2N_dxideta, 
                                       CfemReal* d2N_dxidzeta, 
                                       CfemReal* d2N_detadzeta) const override 
        {
            // Nœud 0
            d2N_dxi2[0]    = 4.0; d2N_deta2[0]    = 4.0; d2N_dzeta2[0]    = 4.0;
            d2N_dxideta[0] = 4.0; d2N_dxidzeta[0] = 4.0; d2N_detadzeta[0] = 4.0;

            // Nœud 1
            d2N_dxi2[1]    = 4.0; d2N_deta2[1] = 0.0;    d2N_dzeta2[1]    = 0.0;
            d2N_dxideta[1] = 0.0; d2N_dxidzeta[1] = 0.0; d2N_detadzeta[1] = 0.0;

            // Nœud 2
            d2N_dxi2[2] = 0.0; d2N_deta2[2] = 4.0; d2N_dzeta2[2] = 0.0;
            d2N_dxideta[2] = 0.0; d2N_dxidzeta[2] = 0.0; d2N_detadzeta[2] = 0.0;

            // Nœud 3
            d2N_dxi2[3] = 0.0; d2N_deta2[3] = 0.0; d2N_dzeta2[3] = 4.0;
            d2N_dxideta[3] = 0.0; d2N_dxidzeta[3] = 0.0; d2N_detadzeta[3] = 0.0;

            // Nœud 4
            d2N_dxi2[4] =-8.0; d2N_deta2[4] = 0.0; d2N_dzeta2[4] = 0.0;
            d2N_dxideta[4] =-4.0; d2N_dxidzeta[4] =-4.0; d2N_detadzeta[4] = 0.0;

            // Nœud 5
            d2N_dxi2[5] = 0.0; d2N_deta2[5] = 0.0; d2N_dzeta2[5] = 0.0;
            d2N_dxideta[5] = 4.0; d2N_dxidzeta[5] = 0.0; d2N_detadzeta[5] = 0.0;

            // Nœud 6
            d2N_dxi2[6] = 0.0; d2N_deta2[6] =-8.0; d2N_dzeta2[6] = 0.0;
            d2N_dxideta[6] =-4.0; d2N_dxidzeta[6] = 0.0; d2N_detadzeta[6] =-4.0;

            // Nœud 7
            d2N_dxi2[7] = 0.0; d2N_deta2[7] = 0.0; d2N_dzeta2[7] =-8.0;
            d2N_dxideta[7] = 0.0; d2N_dxidzeta[7] =-4.0; d2N_detadzeta[7] =-4.0;

            // Nœud 8
            d2N_dxi2[8] = 0.0; d2N_deta2[8] = 0.0; d2N_dzeta2[8] = 0.0;
            d2N_dxideta[8] = 0.0; d2N_dxidzeta[8] = 4.0; d2N_detadzeta[8] = 0.0;

            // Nœud 9
            d2N_dxi2[9] = 0.0; d2N_deta2[9] = 0.0; d2N_dzeta2[9] = 0.0;
            d2N_dxideta[9] = 0.0; d2N_dxidzeta[9] = 0.0; d2N_detadzeta[9] = 4.0;
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, DynamicMatrix<CfemReal>& values) const override
        {
            // Constant piecewise field value
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                const Vector3D& xi =  xiPoints[q];
                const CfemReal L1 = xi.x;
                const CfemReal L2 = xi.y;
                const CfemReal L3 = xi.z;
                const CfemReal L0 = 1.0 - L1 - L2 - L3;

                // Sommets
                values[q][0] = L0 * (2.0 * L0 - 1.0);
                values[q][1] = L1 * (2.0 * L1 - 1.0);
                values[q][2] = L2 * (2.0 * L2 - 1.0);
                values[q][3] = L3 * (2.0 * L3 - 1.0);

                // Arêtes horizontales (base)
                values[q][4] = 4.0 * L0 * L1;
                values[q][5] = 4.0 * L1 * L2;
                values[q][6] = 4.0 * L2 * L0;

                // Arêtes verticales (vers L3)
                values[q][7] = 4.0 * L0 * L3;
                values[q][8] = 4.0 * L1 * L3;
                values[q][9] = 4.0 * L2 * L3;
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> & xiPoints, 
                                      DynamicMatrix<CfemReal> &dN_dxi, 
                                      DynamicMatrix<CfemReal> &dN_deta, 
                                      DynamicMatrix<CfemReal> &dN_dzeta) const override
        {
            // Constant function results in flat zero gradients everywhere
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                const Vector3D& xi =  xiPoints[q];
                const CfemReal L1 = xi.x;
                const CfemReal L2 = xi.y;
                const CfemReal L3 = xi.z;

                // Dérivées par rapport à xi (x)
                dN_dxi[q][0] = -3.0 + 4.0*L1 + 4.0*L2 + 4.0*L3;
                dN_dxi[q][1] =  4.0*L1 - 1.0;
                dN_dxi[q][2] =  0.0;
                dN_dxi[q][3] =  0.0;
                dN_dxi[q][4] =  4.0 - 8.0*L1 - 4.0*L2 - 4.0*L3;
                dN_dxi[q][5] =  4.0*L2;
                dN_dxi[q][6] = -4.0*L2;
                dN_dxi[q][7] = -4.0*L3;
                dN_dxi[q][8] =  4.0*L3;
                dN_dxi[q][9] =  0.0;

                // Dérivées par rapport à eta (y)
                dN_deta[q][0] = -3.0 + 4.0*L1 + 4.0*L2 + 4.0*L3;
                dN_deta[q][1] =  0.0;
                dN_deta[q][2] =  4.0*L2 - 1.0;
                dN_deta[q][3] =  0.0;
                dN_deta[q][4] = -4.0*L1;
                dN_deta[q][5] =  4.0*L1;
                dN_deta[q][6] =  4.0 - 4.0*L1 - 8.0*L2 - 4.0*L3;
                dN_deta[q][7] = -4.0*L3;
                dN_deta[q][8] =  0.0;
                dN_deta[q][9] =  4.0*L3;

                // Dérivées par rapport à zeta (z)
                dN_dzeta[q][0] = -3.0 + 4.0*L1 + 4.0*L2 + 4.0*L3;
                dN_dzeta[q][1] =  0.0;
                dN_dzeta[q][2] =  0.0;
                dN_dzeta[q][3] =  4.0*L3 - 1.0;
                dN_dzeta[q][4] = -4.0*L1;
                dN_dzeta[q][5] =  0.0;
                dN_dzeta[q][6] = -4.0*L2;
                dN_dzeta[q][7] =  4.0 - 4.0*L1 - 4.0*L2 - 8.0*L3;
                dN_dzeta[q][8] =  4.0*L1;
                dN_dzeta[q][9] =  4.0*L2;
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
            for (CfemInt q = 0; q < nQp; ++q) {
                // Nœud 0
                d2N_dxi2[q][0]    = 4.0; d2N_deta2[q][0]    = 4.0; d2N_dzeta2[q][0]    = 4.0;
                d2N_dxideta[q][0] = 4.0; d2N_dxidzeta[q][0] = 4.0; d2N_detadzeta[q][0] = 4.0;

                // Nœud 1
                d2N_dxi2[q][1]    = 4.0; d2N_deta2[q][1]    = 0.0;    d2N_dzeta2[q][1] = 0.0;
                d2N_dxideta[q][1] = 0.0; d2N_dxidzeta[q][1] = 0.0; d2N_detadzeta[q][1] = 0.0;

                // Nœud 2
                d2N_dxi2[q][2]    = 0.0; d2N_deta2[q][2]    = 4.0; d2N_dzeta2[q][2]    = 0.0;
                d2N_dxideta[q][2] = 0.0; d2N_dxidzeta[q][2] = 0.0; d2N_detadzeta[q][2] = 0.0;

                // Nœud 3
                d2N_dxi2[q][3]    = 0.0; d2N_deta2[q][3]    = 0.0; d2N_dzeta2[q][3]    = 4.0;
                d2N_dxideta[q][3] = 0.0; d2N_dxidzeta[q][3] = 0.0; d2N_detadzeta[q][3] = 0.0;

                // Nœud 4
                d2N_dxi2[q][4]    =-8.0; d2N_deta2[q][4]    = 0.0; d2N_dzeta2[q][4]    = 0.0;
                d2N_dxideta[q][4] =-4.0; d2N_dxidzeta[q][4] =-4.0; d2N_detadzeta[q][4] = 0.0;

                // Nœud 5
                d2N_dxi2[q][5]    = 0.0; d2N_deta2[q][5]    = 0.0; d2N_dzeta2[q][5]    = 0.0;
                d2N_dxideta[q][5] = 4.0; d2N_dxidzeta[q][5] = 0.0; d2N_detadzeta[q][5] = 0.0;

                // Nœud 6
                d2N_dxi2[q][6]    = 0.0; d2N_deta2[q][6]    =-8.0; d2N_dzeta2[q][6]    = 0.0;
                d2N_dxideta[q][6] =-4.0; d2N_dxidzeta[q][6] = 0.0; d2N_detadzeta[q][6] =-4.0;

                // Nœud 7
                d2N_dxi2[q][7]    = 0.0; d2N_deta2[q][7]    = 0.0; d2N_dzeta2[q][7]    =-8.0;
                d2N_dxideta[q][7] = 0.0; d2N_dxidzeta[q][7] =-4.0; d2N_detadzeta[q][7] =-4.0;

                // Nœud 8
                d2N_dxi2[q][8]    = 0.0; d2N_deta2[q][8]    = 0.0; d2N_dzeta2[q][8]    = 0.0;
                d2N_dxideta[q][8] = 0.0; d2N_dxidzeta[q][8] = 4.0; d2N_detadzeta[q][8] = 0.0;

                // Nœud 9
                d2N_dxi2[q][9]    = 0.0;  d2N_deta2[q][9]   = 0.0; d2N_dzeta2[q][9]    = 0.0;
                d2N_dxideta[q][9] = 0.0; d2N_dxidzeta[q][9] = 0.0; d2N_detadzeta[q][9] = 4.0;
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToTetraRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tet10_face_nodes[iFace];
        }
    };

    class MiniReferenceTetrahedron : public ReferenceElement {
    public:
        MiniReferenceTetrahedron() {
            m_numNodes = 5; // 4 sommets + 1 bulle
            m_dimension = 3;
            
            m_nodes = { 
                Vector3D{0.0, 0.0 ,0.0},
                Vector3D{1.0, 0.0 ,0.0},
                Vector3D{0.0, 1.0 ,0.0},
                Vector3D{0.0, 0.0 ,1.0},
                Vector3D{0.25, 0.25, 0.25}
            };
        }

        void evaluateShapesValues(const Vector3D& xi, CfemReal* values) const override 
        {
            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;
            const CfemReal L3 = xi.z;
            const CfemReal L0 = 1.0 - L1 - L2 - L3;

            // Bubble function: 256 * L0 * L1 * L2 * L3
            const CfemReal b = 256.0 * L0 * L1 * L2 * L3;

            values[4] = b;
            values[0] = L0 - b * 0.25;
            values[1] = L1 - b * 0.25;
            values[2] = L2 - b * 0.25;
            values[3] = L3 - b * 0.25;
        }

        void evaluateShapesDerivatives(const Vector3D& xi, 
                                 CfemReal* dN_dxi, 
                                 CfemReal* dN_deta, 
                                 CfemReal* dN_dzeta) const override 
        {
            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;
            const CfemReal L3 = xi.z;
            const CfemReal L0 = 1.0 - L1 - L2 - L3;

            // Bubble gradients
            const CfemReal db_dxi  = 256.0 * L2 * L3 * (L0 - L1);
            const CfemReal db_deta = 256.0 * L1 * L3 * (L0 - L2);
            const CfemReal db_dzeta = 256.0 * L1 * L2 * (L0 - L3);

            // Node 4 (Bubble)
            dN_dxi[4]   = db_dxi;
            dN_deta[4]  = db_deta;
            dN_dzeta[4] = db_dzeta;

            // Node 0
            dN_dxi[0]   = -1.0 - db_dxi * 0.25;
            dN_deta[0]  = -1.0 - db_deta * 0.25;
            dN_dzeta[0] = -1.0 - db_dzeta * 0.25;

            // Node 1
            dN_dxi[1]   =  1.0 - db_dxi * 0.25;
            dN_deta[1]  =  0.0 - db_deta * 0.25;
            dN_dzeta[1] =  0.0 - db_dzeta * 0.25;

            // Node 2
            dN_dxi[2]   =  0.0 - db_dxi * 0.25;
            dN_deta[2]  =  1.0 - db_deta * 0.25;
            dN_dzeta[2] =  0.0 - db_dzeta * 0.25;

            // Node 3
            dN_dxi[3]   =  0.0 - db_dxi * 0.25;
            dN_deta[3]  =  0.0 - db_deta * 0.25;
            dN_dzeta[3] =  1.0 - db_dzeta * 0.25;
        }

        void evaluateShapesSecondDerivatives(const Vector3D& xi,
                                       CfemReal* d2N_dxi2, 
                                       CfemReal* d2N_deta2, 
                                       CfemReal* d2N_dzeta2,
                                       CfemReal* d2N_dxideta, 
                                       CfemReal* d2N_dxidzeta, 
                                       CfemReal* d2N_detadzeta) const override 
        {
            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;
            const CfemReal L3 = xi.z;
            const CfemReal L0 = 1.0 - L1 - L2 - L3;

            // Second derivatives of the bubble function
            const CfemReal d2b_dxi2     = -512.0 * L2 * L3;
            const CfemReal d2b_deta2    = -512.0 * L1 * L3;
            const CfemReal d2b_dzeta2   = -512.0 * L1 * L2;
            const CfemReal d2b_dxideta  =  256.0 * L3 * (L0 - L1 - L2);
            const CfemReal d2b_dxidzeta =  256.0 * L2 * (L0 - L1 - L3);
            const CfemReal d2b_detadzeta=  256.0 * L1 * (L0 - L2 - L3);

            // Node 4 (Bubble)
            d2N_dxi2[4]      = d2b_dxi2; 
            d2N_deta2[4]     = d2b_deta2; 
            d2N_dzeta2[4]    = d2b_dzeta2;
            d2N_dxideta[4]   = d2b_dxideta; 
            d2N_dxidzeta[4]  = d2b_dxidzeta; 
            d2N_detadzeta[4] = d2b_detadzeta;

            // Nodes 0 to 3 inherit the bubble's correction divided by -4.0
            for (CfemInt i = 0; i < 4; ++i) {
                d2N_dxi2[i]      = -d2b_dxi2 * 0.25;
                d2N_deta2[i]     = -d2b_deta2 * 0.25;
                d2N_dzeta2[i]    = -d2b_dzeta2 * 0.25;
                d2N_dxideta[i]   = -d2b_dxideta * 0.25;
                d2N_dxidzeta[i]  = -d2b_dxidzeta * 0.25;
                d2N_detadzeta[i] = -d2b_detadzeta * 0.25;
            }
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, DynamicMatrix<CfemReal>& values) const override
        {
            // Constant piecewise field value
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                const CfemReal L1 = xiPoints[q].x;
                const CfemReal L2 = xiPoints[q].y;
                const CfemReal L3 = xiPoints[q].z;
                const CfemReal L0 = 1.0 - L1 - L2 - L3;

                // Bubble function: 256 * L0 * L1 * L2 * L3
                const CfemReal b = 256.0 * L0 * L1 * L2 * L3;

                values[q][4] = b;
                values[q][0] = L0 - b * 0.25;
                values[q][1] = L1 - b * 0.25;
                values[q][2] = L2 - b * 0.25;
                values[q][3] = L3 - b * 0.25;
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> & xiPoints, 
                                      DynamicMatrix<CfemReal> &dN_dxi, 
                                      DynamicMatrix<CfemReal> &dN_deta, 
                                      DynamicMatrix<CfemReal> &dN_dzeta) const override
        {
            // Constant function results in flat zero gradients everywhere
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                const Vector3D& xi =  xiPoints[q];
                const CfemReal L1 = xi.x;
                const CfemReal L2 = xi.y;
                const CfemReal L3 = xi.z;
                const CfemReal L0 = 1.0 - L1 - L2 - L3;

                // Bubble gradients
                const CfemReal db_dxi  = 256.0 * L2 * L3 * (L0 - L1);
                const CfemReal db_deta = 256.0 * L1 * L3 * (L0 - L2);
                const CfemReal db_dzeta = 256.0 * L1 * L2 * (L0 - L3);

                // Node 4 (Bubble)
                dN_dxi[q][4]   = db_dxi;
                dN_deta[q][4]  = db_deta;
                dN_dzeta[q][4] = db_dzeta;

                // Node 0
                dN_dxi[q][0]   = -1.0 - db_dxi * 0.25;
                dN_deta[q][0]  = -1.0 - db_deta * 0.25;
                dN_dzeta[q][0] = -1.0 - db_dzeta * 0.25;

                // Node 1
                dN_dxi[q][1]   =  1.0 - db_dxi * 0.25;
                dN_deta[q][1]  =  0.0 - db_deta * 0.25;
                dN_dzeta[q][1] =  0.0 - db_dzeta * 0.25;

                // Node 2
                dN_dxi[q][2]   =  0.0 - db_dxi * 0.25;
                dN_deta[q][2]  =  1.0 - db_deta * 0.25;
                dN_dzeta[q][2] =  0.0 - db_dzeta * 0.25;

                // Node 3
                dN_dxi[q][3]   =  0.0 - db_dxi * 0.25;
                dN_deta[q][3]  =  0.0 - db_deta * 0.25;
                dN_dzeta[q][3] =  1.0 - db_dzeta * 0.25;
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
            for (CfemInt q = 0; q < nQp; ++q) {
                d2N_dxi2[q][0]      = 0.0;
                d2N_deta2[q][0]     = 0.0;
                d2N_dzeta2[q][0]    = 0.0;
                d2N_dxideta[q][0]   = 0.0;
                d2N_dxidzeta[q][0]  = 0.0;
                d2N_detadzeta[q][0] = 0.0;
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToTetraRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tet4_face_nodes[iFace];
        }
    };

class P0ReferenceTetrahedron : public ReferenceElement {
    public:
        P0ReferenceTetrahedron() {
            m_numNodes = 1;
            m_dimension = 3;
        
            m_nodes = { Vector3D{0.25, 0.25, 0.25} };
        }

        void evaluateShapesValues(const Vector3D&, CfemReal* values) const override 
        {
            values[0] = 1.0;
        }

        void evaluateShapesDerivatives(const Vector3D&, 
                                 CfemReal* dN_dxi, 
                                 CfemReal* dN_deta, 
                                 CfemReal* dN_dzeta) const override 
        {
            dN_dxi[0]   = 0.0;
            dN_deta[0]  = 0.0;
            dN_dzeta[0] = 0.0;
        }

        void evaluateShapesSecondDerivatives(const Vector3D&,
                                       CfemReal* d2N_dxi2, 
                                       CfemReal* d2N_deta2, 
                                       CfemReal* d2N_dzeta2,
                                       CfemReal* d2N_dxideta, 
                                       CfemReal* d2N_dxidzeta, 
                                       CfemReal* d2N_detadzeta) const override 
        {
            d2N_dxi2[0]      = 0.0;
            d2N_deta2[0]     = 0.0;
            d2N_dzeta2[0]    = 0.0;
            d2N_dxideta[0]   = 0.0;
            d2N_dxidzeta[0]  = 0.0;
            d2N_detadzeta[0] = 0.0;
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, DynamicMatrix<CfemReal>& values) const override
        {
            // Constant piecewise field value
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                values[q][0] = 1.0;
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> & xiPoints, 
                                      DynamicMatrix<CfemReal> &dN_dxi, 
                                      DynamicMatrix<CfemReal> &dN_deta, 
                                      DynamicMatrix<CfemReal> &dN_dzeta) const override
        {
            // Constant function results in flat zero gradients everywhere
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                dN_dxi[q][0]   = 0.0;
                dN_deta[q][0]  = 0.0;
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
            for (CfemInt q = 0; q < nQp; ++q) {
                d2N_dxi2[q][0]      = 0.0;
                d2N_deta2[q][0]     = 0.0;
                d2N_dzeta2[q][0]    = 0.0;
                d2N_dxideta[q][0]   = 0.0;
                d2N_dxidzeta[q][0]  = 0.0;
                d2N_detadzeta[q][0] = 0.0;
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToTetraRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tet1_face_nodes[iFace];
        }
    };
}
#endif