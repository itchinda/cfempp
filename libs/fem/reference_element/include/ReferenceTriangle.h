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
   
    static const std::vector<std::vector<int>> tri1_face_nodes = {
        {0}, {0}, {0}
    };

    static const std::vector<std::vector<int>> tri3_face_nodes = {
        {0, 1}, {1, 2}, {2, 0}
    };
    
    static const std::vector<std::vector<int>> tri6_face_nodes = {
        {0, 1, 3}, {1, 2, 4}, {2, 0, 5}
    };

}

namespace cfem {

    /**
     * @class P1ReferenceTriangle
     * @brief Linear 3-node reference triangle.
     */
    class P1ReferenceTriangle : public ReferenceElement {

    public:
        
        P1ReferenceTriangle() {
            m_numNodes = 3;
            m_dimension = 2;

            m_nodes = { Vector3D{0.0, 0.0, 0.0}, 
                        Vector3D{1.0, 0.0, 0.0}, 
                        Vector3D{0.0, 1.0, 0.0}  };
        }

        void evaluateShapesValues(const Vector3D &xi, CfemReal* values) const override
        {
            values[0] = 1.0 - xi.x - xi.y;
            values[1] = xi.x;
            values[2] = xi.y;
        }

        void evaluateShapesDerivatives(const Vector3D &, 
                                 CfemReal* dN_dxi, 
                                 CfemReal* dN_deta, 
                                 CfemReal* dN_dzeta) const override
        {
            dN_dxi[0]  = -1.0; dN_deta[0]  = -1.0; dN_dzeta[0]  = 0.0;
            dN_dxi[1]  =  1.0; dN_deta[1]  =  0.0; dN_dzeta[1]  = 0.0;
            dN_dxi[2]  =  0.0; dN_deta[2]  =  1.0; dN_dzeta[2]  = 0.0;
        }

        void evaluateShapesSecondDerivatives(const Vector3D &,
                                       CfemReal* d2N_dxi2, 
                                       CfemReal* d2N_deta2, 
                                       CfemReal* d2N_dzeta2,
                                       CfemReal* d2N_dxideta, 
                                       CfemReal* d2N_dxidzeta, 
                                       CfemReal* d2N_detadzeta) const override
        {
            for (CfemInt i = 0; i < m_numNodes; ++i) {
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
                const Vector3D& xi = xiPoints[q];
                values[q][0] = 1.0 - xi.x - xi.y;
                values[q][1] = xi.x;
                values[q][2] = xi.y;
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
                const Vector3D& xi = xiPoints[q];
                dN_dxi[q][0] = -1.0; dN_deta[q][0] = -1.0; dN_dzeta[q][0] = 0.0;
                dN_dxi[q][1] =  1.0; dN_deta[q][1] =  0.0; dN_dzeta[q][1] = 0.0;
                dN_dxi[q][2] =  0.0; dN_deta[q][2] =  1.0; dN_dzeta[q][2] = 0.0;
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
                for (CfemInt i = 0; i < m_numNodes; ++i) {
                    d2N_dxi2[q][i]      = 0.0;
                    d2N_deta2[q][i]     = 0.0;
                    d2N_dzeta2[q][i]    = 0.0;
                    d2N_dxideta[q][i]   = 0.0;
                    d2N_dxidzeta[q][i]  = 0.0;
                    d2N_detadzeta[q][i] = 0.0;
                }
            }
        }
        
        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tri3_face_nodes[iFace];
        }

    };


    /**
     * @class P2ReferenceTriangle
     * @brief Quadratic 6-node reference triangle.
     */
    class P2ReferenceTriangle : public ReferenceElement {
    public:
        P2ReferenceTriangle() {
            m_numNodes = 6;
            m_dimension = 2;

            m_nodes = {
                Vector3D{0.0, 0.0, 0.0},
                Vector3D{1.0, 0.0, 0.0},
                Vector3D{0.0, 1.0, 0.0},
                Vector3D{0.5, 0.0, 0.0},
                Vector3D{0.5, 0.5, 0.0},
                Vector3D{0.0, 0.5, 0.0}
            };
        }

        void evaluateShapesValues(const Vector3D &xi, CfemReal* values) const override 
        {
            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;
            const CfemReal L0 = 1.0 - L1 - L2;

            values[0] = L0 * (2.0 * L0 - 1.0);
            values[1] = L1 * (2.0 * L1 - 1.0);
            values[2] = L2 * (2.0 * L2 - 1.0);
            values[3] = 4.0 * L0 * L1;
            values[4] = 4.0 * L1 * L2;
            values[5] = 4.0 * L2 * L0;
        }

        void evaluateShapesDerivatives(const Vector3D &xi, 
                                 CfemReal* dN_dxi, 
                                 CfemReal* dN_deta, 
                                 CfemReal* dN_dzeta) const override
        {
            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;

            // Derivatives with respect to xi (Local x direction)
            dN_dxi[0] =-3.0 + 4.0 * L1 + 4.0 * L2;
            dN_dxi[1] = 4.0 * L1 - 1.0;
            dN_dxi[2] = 0.0;
            dN_dxi[3] = 4.0 - 8.0 * L1 - 4.0 * L2;
            dN_dxi[4] = 4.0 * L2;
            dN_dxi[5] =-4.0 * L2;

            // Derivatives with respect to eta (Local y direction)
            dN_deta[0] = -3.0 + 4.0 * L1 + 4.0 * L2;
            dN_deta[1] = 0.0;
            dN_deta[2] = 4.0 * L2 - 1.0;
            dN_deta[3] = -4.0 * L1;
            dN_deta[4] = 4.0 * L1;
            dN_deta[5] = 4.0 - 4.0 * L1 - 8.0 * L2;

            // Derivatives with respect to zeta (Z-direction is zero for a 2D planar element)
            for (CfemInt i = 0; i < 6; ++i) {
                dN_dzeta[i] = 0.0;
            }
        }

        void evaluateShapesSecondDerivatives(const Vector3D &,
                                       CfemReal* d2N_dxi2, 
                                       CfemReal* d2N_deta2, 
                                       CfemReal* d2N_dzeta2,
                                       CfemReal* d2N_dxideta, 
                                       CfemReal* d2N_dxidzeta, 
                                       CfemReal* d2N_detadzeta) const override 
        {
            // Quadratic shape function second derivatives are constants.
            d2N_dxi2[0]    = 4.0;  d2N_deta2[0]    = 4.0;  d2N_dxideta[0]   = 4.0;
            d2N_dxi2[1]    = 4.0;  d2N_deta2[1]    = 0.0;  d2N_dxideta[1]   = 0.0;
            d2N_dxi2[2]    = 0.0;  d2N_deta2[2]    = 4.0;  d2N_dxideta[2]   = 0.0;
            d2N_dxi2[3]    =-8.0;  d2N_deta2[3]    = 0.0;  d2N_dxideta[3]   =-4.0;
            d2N_dxi2[4]    = 0.0;  d2N_deta2[4]    = 0.0;  d2N_dxideta[4]   = 4.0;
            d2N_dxi2[5]    = 0.0;  d2N_deta2[5]    =-8.0;  d2N_dxideta[5]   =-4.0;

            // All derivatives involving the out-of-plane reference coordinate zeta are strictly zero
            for (CfemInt i = 0; i < 6; ++i) {
                d2N_dzeta2[i]    = 0.0;
                d2N_dxidzeta[i]  = 0.0;
                d2N_detadzeta[i] = 0.0;
            }
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, DynamicMatrix<CfemReal>& values) const override
        {
            // Constant piecewise field value
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                const Vector3D& xi = xiPoints[q];
                const CfemReal L1 = xi.x;
                const CfemReal L2 = xi.y;
                const CfemReal L0 = 1.0 - L1 - L2;

                values[q][0] = L0 * (2.0 * L0 - 1.0);
                values[q][1] = L1 * (2.0 * L1 - 1.0);
                values[q][2] = L2 * (2.0 * L2 - 1.0);
                values[q][3] = 4.0 * L0 * L1;
                values[q][4] = 4.0 * L1 * L2;
                values[q][5] = 4.0 * L2 * L0;
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
                const Vector3D& xi = xiPoints[q];
                const CfemReal L1 = xi.x;
                const CfemReal L2 = xi.y;

                // Derivatives with respect to xi (Local x direction)
                dN_dxi[q][0] =-3.0 + 4.0 * L1 + 4.0 * L2;
                dN_dxi[q][1] = 4.0 * L1 - 1.0;
                dN_dxi[q][2] = 0.0;
                dN_dxi[q][3] = 4.0 - 8.0 * L1 - 4.0 * L2;
                dN_dxi[q][4] = 4.0 * L2;
                dN_dxi[q][5] =-4.0 * L2;

                // Derivatives with respect to eta (Local y direction)
                dN_deta[q][0] = -3.0 + 4.0 * L1 + 4.0 * L2;
                dN_deta[q][1] = 0.0;
                dN_deta[q][2] = 4.0 * L2 - 1.0;
                dN_deta[q][3] = -4.0 * L1;
                dN_deta[q][4] = 4.0 * L1;
                dN_deta[q][5] = 4.0 - 4.0 * L1 - 8.0 * L2;

                // Derivatives with respect to zeta (Z-direction is zero for a 2D planar element)
                for (CfemInt i = 0; i < 6; ++i) {
                    dN_dzeta[q][i] = 0.0;
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
            // Flat null second derivatives
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                // Quadratic shape function second derivatives are constants.
                d2N_dxi2[q][0] = 4.0;  d2N_deta2[q][0] = 4.0;  d2N_dxideta[q][0] = 4.0;
                d2N_dxi2[q][1] = 4.0;  d2N_deta2[q][1] = 0.0;  d2N_dxideta[q][1] = 0.0;
                d2N_dxi2[q][2] = 0.0;  d2N_deta2[q][2] = 4.0;  d2N_dxideta[q][2] = 0.0;
                d2N_dxi2[q][3] =-8.0;  d2N_deta2[q][3] = 0.0;  d2N_dxideta[q][3] =-4.0;
                d2N_dxi2[q][4] = 0.0;  d2N_deta2[q][4] = 0.0;  d2N_dxideta[q][4] = 4.0;
                d2N_dxi2[q][5] = 0.0;  d2N_deta2[q][5] =-8.0;  d2N_dxideta[q][5] =-4.0;

                // All derivatives involving the out-of-plane reference coordinate zeta are strictly zero
                for (CfemInt i = 0; i < 6; ++i) {
                    d2N_dzeta2[q][i]    = 0.0;
                    d2N_dxidzeta[q][i]  = 0.0;
                    d2N_detadzeta[q][i] = 0.0;
                }
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tri6_face_nodes[iFace];
        }
    };

    class MiniReferenceTriangle : public ReferenceElement {
    private:
        static constexpr inline CfemReal m_1_3 {1.0 / 3.0};
    public:
        MiniReferenceTriangle() {
            m_numNodes = 4; // 3 sommets + 1 bulle
            m_dimension = 2;

            m_nodes = {
                Vector3D{0.0, 0.0, 0.0},
                Vector3D{1.0, 0.0, 0.0},
                Vector3D{0.0, 1.0, 0.0},
                Vector3D{m_1_3, m_1_3, 0.0}
            };
        }

        void evaluateShapesValues(const Vector3D &xi, CfemReal* values) const override
        {
            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;
            const CfemReal L0 = 1.0 - L1 - L2;

            const CfemReal b = 27.0 * L0 * L1 * L2;
            const CfemReal b_x_m_1_3 = b * m_1_3;

            values[3] = b;
            values[0] = L0 - b_x_m_1_3;
            values[1] = L1 - b_x_m_1_3;
            values[2] = L2 - b_x_m_1_3;
        }

        void evaluateShapesDerivatives(const Vector3D &xi, 
                                 CfemReal* dN_dxi, 
                                 CfemReal* dN_deta, 
                                 CfemReal* dN_dzeta) const override
        {
            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;

            // Evaluate bubble gradients
            const CfemReal db_dxi  = 27.0 * L2 * (1.0 - 2.0 * L1 - L2);
            const CfemReal db_deta = 27.0 * L1 * (1.0 - L1 - 2.0 * L2);

            // Derivatives with respect to xi
            dN_dxi[3] = db_dxi;
            dN_dxi[0] = -1.0 - db_dxi * m_1_3;
            dN_dxi[1] =  1.0 - db_dxi * m_1_3;
            dN_dxi[2] =  0.0 - db_dxi * m_1_3;

            // Derivatives with respect to eta
            dN_deta[3] = db_deta;
            dN_deta[0] = -1.0 - db_deta * m_1_3;
            dN_deta[1] =  0.0 - db_deta * m_1_3;
            dN_deta[2] =  1.0 - db_deta * m_1_3;

            // Out-of-plane derivatives are zero
            for (CfemInt i = 0; i < 4; ++i) {
                dN_dzeta[i] = 0.0;
            }
        }

        void evaluateShapesSecondDerivatives(const Vector3D &xi,
                                       CfemReal* d2N_dxi2, 
                                       CfemReal* d2N_deta2, 
                                       CfemReal* d2N_dzeta2,
                                       CfemReal* d2N_dxideta, 
                                       CfemReal* d2N_dxidzeta, 
                                       CfemReal* d2N_detadzeta) const override
        {
            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;

            // Second derivatives of the bubble function
            const CfemReal d2b_dxi2  = -54.0 * L2;
            const CfemReal d2b_deta2 = -54.0 * L1;
            const CfemReal d2b_dxideta = 27.0 * (1.0 - 2.0 * L1 - 2.0 * L2);

            // Node 3 (Bubble)
            d2N_dxi2[3]    = d2b_dxi2;
            d2N_deta2[3]   = d2b_deta2;
            d2N_dxideta[3] = d2b_dxideta;

            // Nodes 0, 1, 2 inherit the bubble's linear correction term
            for (CfemInt i = 0; i < 3; ++i) {
                d2N_dxi2[i]    = -d2b_dxi2 * m_1_3;
                d2N_deta2[i]   = -d2b_deta2 * m_1_3;
                d2N_dxideta[i] = -d2b_dxideta * m_1_3;
            }

            // Initialize all remaining out-of-plane and zeta cross-derivatives to zero
            for (CfemInt i = 0; i < 4; ++i) {
                d2N_dzeta2[i]    = 0.0;
                d2N_dxidzeta[i]  = 0.0;
                d2N_detadzeta[i] = 0.0;
            }
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, DynamicMatrix<CfemReal>& values) const override
        {
            // Constant piecewise field value
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                const Vector3D& xi =  xiPoints[q];
                const CfemReal L1 = xi.x;
                const CfemReal L2 = xi.y;
                const CfemReal L0 = 1.0 - L1 - L2;

                const CfemReal b = 27.0 * L0 * L1 * L2;

                values[q][3] = b;
                values[q][0] = L0 - b * m_1_3;
                values[q][1] = L1 - b * m_1_3;
                values[q][2] = L2 - b * m_1_3;
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

                // Evaluate bubble gradients
                const CfemReal db_dxi  = 27.0 * L2 * (1.0 - 2.0 * L1 - L2);
                const CfemReal db_deta = 27.0 * L1 * (1.0 - L1 - 2.0 * L2);

                // Derivatives with respect to xi
                dN_dxi[q][3] = db_dxi;
                dN_dxi[q][0] = -1.0 - db_dxi * m_1_3;
                dN_dxi[q][1] =  1.0 - db_dxi * m_1_3;
                dN_dxi[q][2] =  0.0 - db_dxi * m_1_3;

                // Derivatives with respect to eta
                dN_deta[q][3] = db_deta;
                dN_deta[q][0] = -1.0 - db_deta * m_1_3;
                dN_deta[q][1] =  0.0 - db_deta * m_1_3;
                dN_deta[q][2] =  1.0 - db_deta * m_1_3;

                // Out-of-plane derivatives are zero
                for (CfemInt i = 0; i < 4; ++i) {
                    dN_dzeta[q][i] = 0.0;
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
            // Flat null second derivatives
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                const Vector3D& xi =  xiPoints[q];
                const CfemReal L1 = xi.x;
                const CfemReal L2 = xi.y;

                // Second derivatives of the bubble function
                const CfemReal d2b_dxi2  = -54.0 * L2;
                const CfemReal d2b_deta2 = -54.0 * L1;
                const CfemReal d2b_dxideta = 27.0 * (1.0 - 2.0 * L1 - 2.0 * L2);

                // Node 3 (Bubble)
                d2N_dxi2[q][3]    = d2b_dxi2;
                d2N_deta2[q][3]   = d2b_deta2;
                d2N_dxideta[q][3] = d2b_dxideta;

                // Nodes 0, 1, 2 inherit the bubble's linear correction term
                for (CfemInt i = 0; i < 3; ++i) {
                    d2N_dxi2[q][i]    = -d2b_dxi2 * m_1_3;
                    d2N_deta2[q][i]   = -d2b_deta2 * m_1_3;
                    d2N_dxideta[q][i] = -d2b_dxideta * m_1_3;
                }

                // Initialize all remaining out-of-plane and zeta cross-derivatives to zero
                for (CfemInt i = 0; i < 4; ++i) {
                    d2N_dzeta2[q][i]    = 0.0;
                    d2N_dxidzeta[q][i]  = 0.0;
                    d2N_detadzeta[q][i] = 0.0;
                }
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet);
        }

        // On réutilise tri3_face_nodes car la bulle n'apparaît pas sur les faces !
        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tri3_face_nodes[iFace];
        }
    };

    class P0ReferenceTriangle : public ReferenceElement {
    public:
        P0ReferenceTriangle() {
            m_numNodes = 1;
            m_dimension = 2;

            m_nodes = { Vector3D{1.0/3.0, 1.0/3.0, 0.0} };
        }

        void evaluateShapesValues(const Vector3D &, CfemReal* values) const override
        {
            // Constant piecewise field value
            values[0] = 1.0;
        }

        void evaluateShapesDerivatives(const Vector3D &, 
                                 CfemReal* dN_dxi, 
                                 CfemReal* dN_deta, 
                                 CfemReal* dN_dzeta) const override
        {
            // Constant function results in flat zero gradients everywhere
            dN_dxi[0]   = 0.0;
            dN_deta[0]  = 0.0;
            dN_dzeta[0] = 0.0;
        }

        void evaluateShapesSecondDerivatives(const Vector3D &,
                                       CfemReal* d2N_dxi2, 
                                       CfemReal* d2N_deta2, 
                                       CfemReal* d2N_dzeta2,
                                       CfemReal* d2N_dxideta, 
                                       CfemReal* d2N_dxidzeta, 
                                       CfemReal* d2N_detadzeta) const override
        {
            // Flat null second derivatives
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
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tri1_face_nodes[iFace];
        }
    };


}  // namespace


#endif