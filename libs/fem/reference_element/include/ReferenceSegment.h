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

#ifndef CFEM_REFERENCE_SEGMENT
#define CFEM_REFERENCE_SEGMENT

#include "ReferenceElement.h"

namespace cfem::reference_elem::utils
{

    inline Vector3D projectFacetRefCoordToSegmentRefCoord(CfemInt localFacetId,
                                                          const Vector3D &xiFacet)
    {
        // Une facette d'une ligne est un point. xiFacet est ignoré.
        // Selon ReferenceSegment.h (P2), le segment de référence est [-1, 1]
        // Face 0 = Noeud 0 (-1.0), Face 1 = Noeud 1 (1.0)
        return Vector3D((localFacetId == 0) ? -1.0 : 1.0, 0.0, 0.0);
    }

    static const std::vector<std::vector<int>> line_face_nodes = {
        {0}, // Face 0
        {1}  // Face 1
    };

    static const std::vector<std::vector<int>> line1_face_nodes = {
        {0}, // Face 0 (x = -1) -> Noeud 0
        {0}  // Face 1 (x =  1) -> Noeud 0
    };
}

namespace cfem
{

    class P1ReferenceSegment : public ReferenceElement
    {
    public:
        P1ReferenceSegment()
        {
            m_numNodes = 2;
            m_dimension = 1;
            m_nodes.resize(m_numNodes);
            // Domaine de référence [-1, 1]
            m_nodes[0] = Vector3D{-1.0, 0.0, 0.0};
            m_nodes[1] = Vector3D{1.0, 0.0, 0.0};
        }

        void evaluateShapesValues(const Vector3D &xi,
                            CfemReal *values) const override
        {
            const CfemReal x = xi.x;
            values[0] = 0.5 * (1.0 - x);
            values[1] = 0.5 * (1.0 + x);
        }

        void evaluateShapesDerivatives(const Vector3D &xi,
                                 CfemReal *dN_dxi,
                                 CfemReal *dN_deta,
                                 CfemReal *dN_dzeta) const override
        {
            // Dérivées / x
            dN_dxi[0] = -0.5;
            dN_dxi[1] = 0.5;

            // L'élément est 1D : dérivées y et z sont nulles
            dN_deta[0] = 0.0;
            dN_deta[1] = 0.0;
            dN_dzeta[0] = 0.0;
            dN_dzeta[1] = 0.0;
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xi,
                                       CfemReal *d2N_dxi2,
                                       CfemReal *d2N_deta2,
                                       CfemReal *d2N_dzeta2,
                                       CfemReal *d2N_dxideta,
                                       CfemReal *d2N_dxidzeta,
                                       CfemReal *d2N_detadzeta) const override
        {
            for (CfemInt i = 0; i < 2; ++i)
            {
                d2N_dxi2[i] = 0.0;
                d2N_deta2[i] = 0.0;
                d2N_dzeta2[i] = 0.0;
                d2N_dxideta[i] = 0.0;
                d2N_dxidzeta[i] = 0.0;
                d2N_detadzeta[i] = 0.0;
            }
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints,
                                 DynamicMatrix<CfemReal> &valuesBatch) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                const CfemReal x = xiPoints[q].x;
                valuesBatch[q][0] = 0.5 * (1.0 - x);
                valuesBatch[q][1] = 0.5 * (1.0 + x);
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
                // Dérivées / x
                dN_dxi[q][0] = -0.5;
                dN_dxi[q][1] = 0.5;

                // L'élément est 1D : dérivées y et z sont nulles
                dN_deta[q][0] = 0.0;
                dN_deta[q][1] = 0.0;
                dN_dzeta[q][0] = 0.0;
                dN_dzeta[q][1] = 0.0;
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
                for (CfemInt i = 0; i < 2; ++i)
                {
                    d2N_dxi2[q][i] = 0.0;
                    d2N_deta2[q][i] = 0.0;
                    d2N_dzeta2[q][i] = 0.0;
                    d2N_dxideta[q][i] = 0.0;
                    d2N_dxidzeta[q][i] = 0.0;
                    d2N_detadzeta[q][i] = 0.0;
                }
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override
        {
            return cfem::reference_elem::utils::projectFacetRefCoordToSegmentRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override
        {
            return cfem::reference_elem::utils::line_face_nodes[iFace];
        }
    };

    class P2ReferenceSegment : public ReferenceElement
    {
    public:
        P2ReferenceSegment()
        {
            m_numNodes = 3;
            m_dimension = 1;
            m_nodes.resize(m_numNodes);
            // Domaine de référence [-1, 1], ordre VTK : extrémités puis milieu
            m_nodes[0] = Vector3D{-1.0, 0.0, 0.0};
            m_nodes[1] = Vector3D{1.0, 0.0, 0.0};
            m_nodes[2] = Vector3D{0.0, 0.0, 0.0}; // Nœud central en x=0
        }

        void evaluateShapesValues(const Vector3D &xi,
                            CfemReal *values) const override
        {
            const CfemReal x = xi.x;
            const CfemReal x2 = x * x;

            values[0] = 0.5 * (x2 - x);
            values[1] = 0.5 * (x2 + x);
            values[2] = 1.0 - x2;
        }

        void evaluateShapesDerivatives(const Vector3D &xi,
                                 CfemReal *dN_dxi,
                                 CfemReal *dN_deta,
                                 CfemReal *dN_dzeta) const override
        {
            const CfemReal x = xi.x;

            dN_dxi[0] = x - 0.5;
            dN_dxi[1] = x + 0.5;
            dN_dxi[2] = -2.0 * x;

            for (CfemInt i = 0; i < 3; ++i)
            {
                dN_deta[i] = 0.0;
                dN_dzeta[i] = 0.0;
            }
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xiPoints,
                                       CfemReal *d2N_dxi2,
                                       CfemReal *d2N_deta2,
                                       CfemReal *d2N_dzeta2,
                                       CfemReal *d2N_dxideta,
                                       CfemReal *d2N_dxidzeta,
                                       CfemReal *d2N_detadzeta) const override
        {

            d2N_dxi2[0] = 1.0;
            d2N_dxi2[1] = 1.0;
            d2N_dxi2[2] = -2.0;

            for (CfemInt i = 0; i < 3; ++i)
            {
                d2N_deta2[i] = 0.0;
                d2N_dzeta2[i] = 0.0;
                d2N_dxideta[i] = 0.0;
                d2N_dxidzeta[i] = 0.0;
                d2N_detadzeta[i] = 0.0;
            }
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints,
                                 DynamicMatrix<CfemReal> &valuesBatch) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {
                const CfemReal x = xiPoints[q].x;
                const CfemReal x2 = x * x;

                valuesBatch[q][0] = 0.5 * (x2 - x);
                valuesBatch[q][1] = 0.5 * (x2 + x);
                valuesBatch[q][2] = 1.0 - x2;
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
                const CfemReal x = xiPoints[q].x;

                dN_dxi[q][0] = x - 0.5;
                dN_dxi[q][1] = x + 0.5;
                dN_dxi[q][2] = -2.0 * x;

                for (CfemInt i = 0; i < 3; ++i)
                {
                    dN_deta[q][i] = 0.0;
                    dN_dzeta[q][i] = 0.0;
                }
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints,
                                            DynamicMatrix<CfemReal> &d2N_dxi2, DynamicMatrix<CfemReal> &d2N_deta2, DynamicMatrix<CfemReal> &d2N_dzeta2,
                                            DynamicMatrix<CfemReal> &d2N_dxideta, DynamicMatrix<CfemReal> &d2N_dxidzeta, DynamicMatrix<CfemReal> &d2N_detadzeta) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q)
            {

                d2N_dxi2[q][0] = 1.0;
                d2N_dxi2[q][1] = 1.0;
                d2N_dxi2[q][2] = -2.0;

                for (CfemInt i = 0; i < 3; ++i)
                {
                    d2N_deta2[q][i] = 0.0;
                    d2N_dzeta2[q][i] = 0.0;
                    d2N_dxideta[q][i] = 0.0;
                    d2N_dxidzeta[q][i] = 0.0;
                    d2N_detadzeta[q][i] = 0.0;
                }
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D &xiFacet) const override
        {
            return cfem::reference_elem::utils::projectFacetRefCoordToSegmentRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override
        {
            return cfem::reference_elem::utils::line_face_nodes[iFace];
        }
    };

    class P0ReferenceSegment : public ReferenceElement
    {
    public:
        P0ReferenceSegment()
        {
            m_numNodes = 1;
            m_dimension = 1;
            m_nodes.resize(m_numNodes);
            // Un seul nœud au centre du segment de référence [-1, 1]
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
            return cfem::reference_elem::utils::projectFacetRefCoordToSegmentRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override
        {
            return cfem::reference_elem::utils::line1_face_nodes[iFace];
        }
    };
}

#endif