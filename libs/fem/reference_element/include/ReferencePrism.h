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

    static const std::vector<std::vector<int>> prism6_face_nodes = {
        {0, 2, 1},       // Face 0 (Tri)
        {3, 4, 5},       // Face 1 (Tri)
        {0, 1, 4, 3},    // Face 2 (Quad)
        {1, 2, 5, 4},    // Face 3 (Quad)
        {2, 0, 3, 5}     // Face 4 (Quad)
    };

    static const std::vector<std::vector<int>> prism18_face_nodes = {
        {0, 2, 1, 8, 7, 6},                         // Face 0 (Tri6)
        {3, 4, 5, 9, 10, 11},                       // Face 1 (Tri6)
        {0, 1, 4, 3, 6, 13, 9, 12, 15},             // Face 2 (Quad9)
        {1, 2, 5, 4, 7, 14, 10, 13, 16},            // Face 3 (Quad9)
        {2, 0, 3, 5, 8, 12, 11, 14, 17}             // Face 4 (Quad9)
    };

    static const std::vector<std::vector<int>> prism1_face_nodes = {
        {0}, {0}, {0}, {0}, {0}
    };
}

namespace cfem{

    class P1ReferencePrism : public ReferenceElement {
    private:
        static constexpr CfemInt node_map[3][2] = {
            {0, 3}, // Nœuds du Triangle 0
            {1, 4}, // Nœuds du Triangle 1
            {2, 5}  // Nœuds du Triangle 2
        };
    public:
        P1ReferencePrism() {
            m_numNodes = 6;
            m_dimension = 3;
            m_nodes.resize(m_numNodes);
            // Ordre VTK : Triangle inférieur Z=-1 (0,1,2), puis supérieur Z=1 (3,4,5)
            m_nodes[0] = {0,0,-1}; m_nodes[1] = {1,0,-1}; m_nodes[2] = {0,1,-1};
            m_nodes[3] = {0,0, 1}; m_nodes[4] = {1,0, 1}; m_nodes[5] = {0,1, 1};
        }

        void evaluateShapesValues(const Vector3D& xi, 
                            CfemReal* values) const override 
        {
            constexpr CfemReal o2 = 0.5;
            const CfemReal x = xi.x;
            const CfemReal y = xi.y;
            const CfemReal z = xi.z;

            // Fonctions du Triangle P1
            const CfemReal Nt[3] = {1.0 - x - y, x, y};
            // Fonctions de la Ligne P1
            const CfemReal Nl[2] = {o2 * (1.0 - z), o2 * (1.0 + z)};

            for (CfemInt t = 0; t < 3; ++t) {
                for (CfemInt l = 0; l < 2; ++l) {
                    values[node_map[t][l]] = Nt[t] * Nl[l];
                }
            }
        }

        void evaluateShapesDerivatives(const Vector3D& xi, 
                                 CfemReal* dN_dxi, 
                                 CfemReal* dN_deta, 
                                 CfemReal* dN_dzeta) const override 
        {

            constexpr CfemReal o2 = 0.5;

            const CfemReal x = xi.x;
            const CfemReal y = xi.y;
            const CfemReal z = xi.z;

            const CfemReal Nt[3] = {1.0 - x - y, x, y};
            const CfemReal dNt_dx[3] = {-1.0, 1.0, 0.0};
            const CfemReal dNt_dy[3] = {-1.0, 0.0, 1.0};

            const CfemReal Nl[2] = {o2 * (1.0 - z), o2 * (1.0 + z)};
            const CfemReal dNl_dz[2] = {-o2, o2};

            for (CfemInt t = 0; t < 3; ++t) {
                for (CfemInt l = 0; l < 2; ++l) {
                    const CfemInt node = node_map[t][l];
                    dN_dxi[node]   = dNt_dx[t] * Nl[l];
                    dN_deta[node]  = dNt_dy[t] * Nl[l];
                    dN_dzeta[node] = Nt[t]     * dNl_dz[l];
                }
            }
        }

        void evaluateShapesSecondDerivatives(const Vector3D& xiPoints,
                                       CfemReal* d2N_dxi2,
                                       CfemReal* d2N_deta2,
                                       CfemReal* d2N_dzeta2,
                                       CfemReal* d2N_dxideta,
                                       CfemReal* d2N_dxidzeta,
                                       CfemReal* d2N_detadzeta) const override 
        {

            constexpr CfemReal o2 = 0.5;

            const CfemReal dNt_dx[3] = {-1.0, 1.0, 0.0};
            const CfemReal dNt_dy[3] = {-1.0, 0.0, 1.0};
            const CfemReal dNl_dz[2] = {-o2, o2};

            for (CfemInt t = 0; t < 3; ++t) {
                for (CfemInt l = 0; l < 2; ++l) {
                    const CfemInt node = node_map[t][l];
                        
                    // Les dérivées secondes pures sont nulles pour un champ linéaire (P1xP1)
                    d2N_dxi2[node]    = 0.0;
                    d2N_deta2[node]   = 0.0;
                    d2N_dzeta2[node]  = 0.0;
                    d2N_dxideta[node] = 0.0;

                    // Les dérivées croisées impliquant 'z' survivent au produit tensoriel !
                    d2N_dxidzeta[node]  = dNt_dx[t] * dNl_dz[l];
                    d2N_detadzeta[node] = dNt_dy[t] * dNl_dz[l];
                }
            }
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, 
                                 DynamicMatrix<CfemReal>& valuesBatch) const override 
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            constexpr CfemReal o2 = 0.5;

            for (CfemInt q = 0; q < nQp; ++q) {
                const CfemReal x = xiPoints[q].x;
                const CfemReal y = xiPoints[q].y;
                const CfemReal z = xiPoints[q].z;

                // Fonctions du Triangle P1
                const CfemReal Nt[3] = {1.0 - x - y, x, y};
                // Fonctions de la Ligne P1
                const CfemReal Nl[2] = {o2 * (1.0 - z), o2 * (1.0 + z)};

                for (CfemInt t = 0; t < 3; ++t) {
                    for (CfemInt l = 0; l < 2; ++l) {
                        valuesBatch[q][node_map[t][l]] = Nt[t] * Nl[l];
                    }
                }
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, 
                                      DynamicMatrix<CfemReal>& dN_dxi, 
                                      DynamicMatrix<CfemReal>& dN_deta, 
                                      DynamicMatrix<CfemReal>& dN_dzeta) const override 
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            constexpr CfemReal o2 = 0.5;

            for (CfemInt q = 0; q < nQp; ++q) {
                const CfemReal x = xiPoints[q].x;
                const CfemReal y = xiPoints[q].y;
                const CfemReal z = xiPoints[q].z;

                const CfemReal Nt[3] = {1.0 - x - y, x, y};
                const CfemReal dNt_dx[3] = {-1.0, 1.0, 0.0};
                const CfemReal dNt_dy[3] = {-1.0, 0.0, 1.0};

                const CfemReal Nl[2] = {o2 * (1.0 - z), o2 * (1.0 + z)};
                const CfemReal dNl_dz[2] = {-o2, o2};

                for (CfemInt t = 0; t < 3; ++t) {
                    for (CfemInt l = 0; l < 2; ++l) {
                        const CfemInt node = node_map[t][l];
                        dN_dxi[q][node]   = dNt_dx[t] * Nl[l];
                        dN_deta[q][node]  = dNt_dy[t] * Nl[l];
                        dN_dzeta[q][node] = Nt[t]     * dNl_dz[l];
                    }
                }
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
            constexpr CfemReal o2 = 0.5;

            for (CfemInt q = 0; q < nQp; ++q) {
                const CfemReal dNt_dx[3] = {-1.0, 1.0, 0.0};
                const CfemReal dNt_dy[3] = {-1.0, 0.0, 1.0};
                const CfemReal dNl_dz[2] = {-o2, o2};

                for (CfemInt t = 0; t < 3; ++t) {
                    for (CfemInt l = 0; l < 2; ++l) {
                        const CfemInt node = node_map[t][l];
                        
                        // Les dérivées secondes pures sont nulles pour un champ linéaire (P1xP1)
                        d2N_dxi2[q][node]    = 0.0;
                        d2N_deta2[q][node]   = 0.0;
                        d2N_dzeta2[q][node]  = 0.0;
                        d2N_dxideta[q][node] = 0.0;

                        // Les dérivées croisées impliquant 'z' survivent au produit tensoriel !
                        d2N_dxidzeta[q][node]  = dNt_dx[t] * dNl_dz[l];
                        d2N_detadzeta[q][node] = dNt_dy[t] * dNl_dz[l];
                    }
                }
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToPrismRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::prism6_face_nodes[iFace];
        }
    };


    class P2ReferencePrism : public ReferenceElement {
    protected:
        // Mapping strict selon la connectivité de P2ReferencePrism
        // t (Triangle P2) : 0..5
        // l (Ligne P2)    : 0 (-1), 1 (0), 2 (+1)
        static constexpr CfemInt node_map[6][3] = {
            {0, 12, 3}, // Sommet 0 (0,0)
            {1, 13, 4}, // Sommet 1 (1,0)
            {2, 14, 5}, // Sommet 2 (0,1)
            {6, 15, 9}, // Milieu  3 (0.5,0)
            {7, 16, 10},// Milieu  4 (0.5,0.5)
            {8, 17, 11} // Milieu  5 (0,0.5)
        };

        // Helper: Triangle P2
        static inline void eval2D_tri_all(CfemReal x, CfemReal y, 
                                          CfemReal* N, CfemReal* dNx, CfemReal* dNy, 
                                          CfemReal* d2Nx, CfemReal* d2Ny, CfemReal* d2Nxy) 
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

        // Helper: Ligne P2 (Ordre -1.0, 0.0, 1.0)
        static inline void eval1D_line_all(CfemReal z, CfemReal* N, CfemReal* dNz, CfemReal* d2Nz) 
        {
            const CfemReal z2 = z * z;
            N[0] = 0.5 * (z2 - z); dNz[0] = z - 0.5;  d2Nz[0] = 1.0;
            N[1] = 1.0 - z2;       dNz[1] = -2.0 * z; d2Nz[1] = -2.0;
            N[2] = 0.5 * (z2 + z); dNz[2] = z + 0.5;  d2Nz[2] = 1.0;
        }

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

        void evaluateShapesValues(const Vector3D& xi, 
                            CfemReal* values) const override 
        {
            CfemReal Nt[6], dNtx[6], dNty[6], d2Ntx[6], d2Nty[6], d2Ntxy[6];
            CfemReal Nl[3], dNlz[3], d2Nlz[3];

            eval2D_tri_all(xi.x, xi.y, Nt, dNtx, dNty, d2Ntx, d2Nty, d2Ntxy);
            eval1D_line_all(xi.z, Nl, dNlz, d2Nlz);

            for (CfemInt t = 0; t < 6; ++t) {
                for (CfemInt l = 0; l < 3; ++l) {
                    values[node_map[t][l]] = Nt[t] * Nl[l];
                }
            }
        }

        void evaluateShapesDerivatives(const Vector3D& xi, 
                                 CfemReal* dN_dxi, 
                                 CfemReal* dN_deta, 
                                 CfemReal* dN_dzeta) const override 
        {
            CfemReal Nt[6], dNtx[6], dNty[6], d2Ntx[6], d2Nty[6], d2Ntxy[6];
            CfemReal Nl[3], dNlz[3], d2Nlz[3];

            eval2D_tri_all(xi.x, xi.y, Nt, dNtx, dNty, d2Ntx, d2Nty, d2Ntxy);
            eval1D_line_all(xi.z, Nl, dNlz, d2Nlz);

            for (CfemInt t = 0; t < 6; ++t) {
                for (CfemInt l = 0; l < 3; ++l) {
                    const CfemInt node = node_map[t][l];
                    dN_dxi[node]   = dNtx[t] * Nl[l];
                    dN_deta[node]  = dNty[t] * Nl[l];
                    dN_dzeta[node] = Nt[t]   * dNlz[l];
                }
            }
        }

        void evaluateShapesSecondDerivatives(const Vector3D& xi,
                                       CfemReal* d2N_dxi2,
                                       CfemReal* d2N_deta2,
                                       CfemReal* d2N_dzeta2,
                                       CfemReal* d2N_dxideta,
                                       CfemReal* d2N_dxidzeta,
                                       CfemReal* d2N_detadzeta) const override 
        {
            CfemReal Nt[6], dNtx[6], dNty[6], d2Ntx[6], d2Nty[6], d2Ntxy[6];
            CfemReal Nl[3], dNlz[3], d2Nlz[3];

            eval2D_tri_all(xi.x, xi.y, Nt, dNtx, dNty, d2Ntx, d2Nty, d2Ntxy);
            eval1D_line_all(xi.z, Nl, dNlz, d2Nlz);

            for (CfemInt t = 0; t < 6; ++t) {
                for (CfemInt l = 0; l < 3; ++l) {
                    const CfemInt node = node_map[t][l];
                    
                    // Règle du produit tensoriel pour les dérivées secondes
                    d2N_dxi2[node]      = d2Ntx[t]  * Nl[l];
                    d2N_deta2[node]     = d2Nty[t]  * Nl[l];
                    d2N_dzeta2[node]    = Nt[t]     * d2Nlz[l];
                        
                    d2N_dxideta[node]   = d2Ntxy[t] * Nl[l];
                    d2N_dxidzeta[node]  = dNtx[t]   * dNlz[l];
                    d2N_detadzeta[node] = dNty[t]   * dNlz[l];
                }
            }
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, 
                                DynamicMatrix<CfemReal>& valuesBatch) const override 
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                CfemReal Nt[6], dNtx[6], dNty[6], d2Ntx[6], d2Nty[6], d2Ntxy[6];
                CfemReal Nl[3], dNlz[3], d2Nlz[3];

                eval2D_tri_all(xiPoints[q].x, xiPoints[q].y, Nt, dNtx, dNty, d2Ntx, d2Nty, d2Ntxy);
                eval1D_line_all(xiPoints[q].z, Nl, dNlz, d2Nlz);

                for (CfemInt t = 0; t < 6; ++t) {
                    for (CfemInt l = 0; l < 3; ++l) {
                        valuesBatch[q][node_map[t][l]] = Nt[t] * Nl[l];
                    }
                }
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, 
                                    DynamicMatrix<CfemReal>& dN_dxi, 
                                    DynamicMatrix<CfemReal>& dN_deta, 
                                    DynamicMatrix<CfemReal>& dN_dzeta) const override 
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                CfemReal Nt[6], dNtx[6], dNty[6], d2Ntx[6], d2Nty[6], d2Ntxy[6];
                CfemReal Nl[3], dNlz[3], d2Nlz[3];

                eval2D_tri_all(xiPoints[q].x, xiPoints[q].y, Nt, dNtx, dNty, d2Ntx, d2Nty, d2Ntxy);
                eval1D_line_all(xiPoints[q].z, Nl, dNlz, d2Nlz);

                for (CfemInt t = 0; t < 6; ++t) {
                    for (CfemInt l = 0; l < 3; ++l) {
                        const CfemInt node = node_map[t][l];
                        dN_dxi[q][node]   = dNtx[t] * Nl[l];
                        dN_deta[q][node]  = dNty[t] * Nl[l];
                        dN_dzeta[q][node] = Nt[t]   * dNlz[l];
                    }
                }
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
            for (CfemInt q = 0; q < nQp; ++q) {
                CfemReal Nt[6], dNtx[6], dNty[6], d2Ntx[6], d2Nty[6], d2Ntxy[6];
                CfemReal Nl[3], dNlz[3], d2Nlz[3];

                eval2D_tri_all(xiPoints[q].x, xiPoints[q].y, Nt, dNtx, dNty, d2Ntx, d2Nty, d2Ntxy);
                eval1D_line_all(xiPoints[q].z, Nl, dNlz, d2Nlz);

                for (CfemInt t = 0; t < 6; ++t) {
                    for (CfemInt l = 0; l < 3; ++l) {
                        const CfemInt node = node_map[t][l];
                        
                        // Règle du produit tensoriel pour les dérivées secondes
                        d2N_dxi2[q][node]      = d2Ntx[t]  * Nl[l];
                        d2N_deta2[q][node]     = d2Nty[t]  * Nl[l];
                        d2N_dzeta2[q][node]    = Nt[t]     * d2Nlz[l];
                        
                        d2N_dxideta[q][node]   = d2Ntxy[t] * Nl[l];
                        d2N_dxidzeta[q][node]  = dNtx[t]   * dNlz[l];
                        d2N_detadzeta[q][node] = dNty[t]   * dNlz[l];
                    }
                }
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToPrismRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::prism18_face_nodes[iFace];
        }
    };

    class MiniReferencePrism : public ReferenceElement {
    public:
        MiniReferencePrism() {
            m_numNodes = 7;
            m_dimension = 3;
            m_nodes.resize(7);
            m_nodes[0]={0,0,-1}; m_nodes[1]={1,0,-1}; m_nodes[2]={0,1,-1};
            m_nodes[3]={0,0, 1}; m_nodes[4]={1,0, 1}; m_nodes[5]={0,1, 1};
            m_nodes[6]={1.0/3.0, 1.0/3.0, 0.0}; // Bulle
        }

        void evaluateShapesValues(const Vector3D& xi, 
                            CfemReal* values) const override 
        {
            constexpr CfemReal o6 = 1.0 / 6.0;

            const CfemReal r = xi.x;
            const CfemReal s = xi.y;
            const CfemReal t = xi.z;

            const CfemReal L[3] = {1.0 - r - s, r, s};
            const CfemReal M[2] = {0.5 * (1.0 - t), 0.5 * (1.0 + t)};

            // Calcul de la Bulle (Nœud 6)
            const CfemReal triB = 27.0 * L[0] * L[1] * L[2];
            const CfemReal zB = 1.0 - t * t;
            const CfemReal B = triB * zB;

            values[6] = B;

            // 2. Fonctions des Sommets (0-5) avec correction
            const CfemReal correction = B * o6;
            for (CfemInt i = 0; i < 6; ++i) {
                const CfemInt iT = i % 3;
                const CfemInt iZ = i / 3;
                values[i] = L[iT] * M[iZ] - correction;
            }
        }

        void evaluateShapesDerivatives(const Vector3D& xi,
                                 CfemReal* dN_dxi, 
                                 CfemReal* dN_deta, 
                                 CfemReal* dN_dzeta) const override 
        {
            constexpr CfemReal o6 = 1.0 / 6.0;

            const CfemReal r = xi.x;
            const CfemReal s = xi.y;
            const CfemReal t = xi.z;

            const CfemReal L[3] = {1.0 - r - s, r, s};
            const CfemReal M[2] = {0.5 * (1.0 - t), 0.5 * (1.0 + t)};
            const CfemReal dM[2] = {-0.5, 0.5};

            // Gradients de la bulle
            const CfemReal triB = 27.0 * L[0] * L[1] * L[2];
            const CfemReal zB = 1.0 - t * t;

            const CfemReal gBx = 27.0 * s * (1.0 - 2.0*r - s) * zB;
            const CfemReal gBy = 27.0 * r * (1.0 - r - 2.0*s) * zB;
            const CfemReal gBz = triB * (-2.0 * t);

            dN_dxi[6]   = gBx;
            dN_deta[6]  = gBy;
            dN_dzeta[6] = gBz;

                // Gradients des sommets avec correction
            const CfemReal c_dx = gBx * o6;
            const CfemReal c_dy = gBy * o6;
            const CfemReal c_dz = gBz * o6;

            for (CfemInt i = 0; i < 6; ++i) {
                const CfemInt iT = i % 3;
                const CfemInt iZ = i / 3;
                const CfemReal drL = (iT == 0 ? -1.0 : (iT == 1 ? 1.0 : 0.0));
                const CfemReal dsL = (iT == 0 ? -1.0 : (iT == 2 ? 1.0 : 0.0));

                dN_dxi[i]   = drL * M[iZ] - c_dx;
                dN_deta[i]  = dsL * M[iZ] - c_dy;
                dN_dzeta[i] = L[iT] * dM[iZ] - c_dz;
            }
        }

        void evaluateShapesSecondDerivatives(const Vector3D& xi,
                                            CfemReal* d2N_dxi2,   
                                            CfemReal* d2N_deta2,
                                            CfemReal* d2N_dzeta2,
                                            CfemReal* d2N_dxideta,
                                            CfemReal* d2N_dxidzeta,
                                            CfemReal* d2N_detadzeta) const override 
        {
            constexpr CfemReal o6 = 1.0 / 6.0;

            const CfemReal r = xi.x;
            const CfemReal s = xi.y;
            const CfemReal t = xi.z;

            const CfemReal L[3] = {1.0 - r - s, r, s};
            const CfemReal dM[2] = {-0.5, 0.5};

            const CfemReal triB = 27.0 * L[0] * L[1] * L[2];
            const CfemReal zB = 1.0 - t * t;

            // Hessiens de la bulle (calculs simplifiés sans divisions par zéro)
            const CfemReal hBxx = -54.0 * s * zB;
            const CfemReal hByy = -54.0 * r * zB;
            const CfemReal hBzz = -2.0 * triB;
            const CfemReal hBxy =  27.0 * (1.0 - 2.0*r - 2.0*s) * zB;
            const CfemReal hBxz = -54.0 * s * t * (1.0 - 2.0*r - s); // CORRECTION: evite (gBx / zB) * (-2t)
            const CfemReal hByz = -54.0 * r * t * (1.0 - r - 2.0*s); // CORRECTION: evite (gBy / zB) * (-2t)

            d2N_dxi2[6]      = hBxx;
            d2N_deta2[6]     = hByy;
            d2N_dzeta2[6]    = hBzz;
            d2N_dxideta[6]   = hBxy;
            d2N_dxidzeta[6]  = hBxz;
            d2N_detadzeta[6] = hByz;

            // Hessiens des sommets avec correction
            const CfemReal c_xx = hBxx * o6;
            const CfemReal c_yy = hByy * o6;
            const CfemReal c_zz = hBzz * o6;
            const CfemReal c_xy = hBxy * o6;
            const CfemReal c_xz = hBxz * o6;
            const CfemReal c_yz = hByz * o6;

            for (CfemInt i = 0; i < 6; ++i) {
                const CfemInt iT = i % 3;
                const CfemInt iZ = i / 3;
                const CfemReal drL = (iT == 0 ? -1.0 : (iT == 1 ? 1.0 : 0.0));
                const CfemReal dsL = (iT == 0 ? -1.0 : (iT == 2 ? 1.0 : 0.0));

                d2N_dxi2[i]      = -c_xx;
                d2N_deta2[i]     = -c_yy;
                d2N_dzeta2[i]    = -c_zz;
                d2N_dxideta[i]   = -c_xy;
                // La partie croisée linéaire-linéaire n'est non nulle que pour les couples impliquant z
                d2N_dxidzeta[i]  = drL * dM[iZ] - c_xz;
                d2N_detadzeta[i] = dsL * dM[iZ] - c_yz;
            }
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, 
                                DynamicMatrix<CfemReal>& valuesBatch) const override 
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            constexpr CfemReal o6 = 1.0 / 6.0;

            for (CfemInt q = 0; q < nQp; ++q) {
                const CfemReal r = xiPoints[q].x;
                const CfemReal s = xiPoints[q].y;
                const CfemReal t = xiPoints[q].z;

                const CfemReal L[3] = {1.0 - r - s, r, s};
                const CfemReal M[2] = {0.5 * (1.0 - t), 0.5 * (1.0 + t)};

                // 1. Calcul de la Bulle (Nœud 6)
                const CfemReal triB = 27.0 * L[0] * L[1] * L[2];
                const CfemReal zB = 1.0 - t * t;
                const CfemReal B = triB * zB;

                valuesBatch[q][6] = B;

                // 2. Fonctions des Sommets (0-5) avec correction
                const CfemReal correction = B * o6;
                for (CfemInt i = 0; i < 6; ++i) {
                    const CfemInt iT = i % 3;
                    const CfemInt iZ = i / 3;
                    valuesBatch[q][i] = L[iT] * M[iZ] - correction;
                }
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, 
                                    DynamicMatrix<CfemReal>& dN_dxi, 
                                    DynamicMatrix<CfemReal>& dN_deta, 
                                    DynamicMatrix<CfemReal>& dN_dzeta) const override 
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            constexpr CfemReal o6 = 1.0 / 6.0;

            for (CfemInt q = 0; q < nQp; ++q) {
                const CfemReal r = xiPoints[q].x;
                const CfemReal s = xiPoints[q].y;
                const CfemReal t = xiPoints[q].z;

                const CfemReal L[3] = {1.0 - r - s, r, s};
                const CfemReal M[2] = {0.5 * (1.0 - t), 0.5 * (1.0 + t)};
                const CfemReal dM[2] = {-0.5, 0.5};

                // Gradients de la bulle
                const CfemReal triB = 27.0 * L[0] * L[1] * L[2];
                const CfemReal zB = 1.0 - t * t;

                const CfemReal gBx = 27.0 * s * (1.0 - 2.0*r - s) * zB;
                const CfemReal gBy = 27.0 * r * (1.0 - r - 2.0*s) * zB;
                const CfemReal gBz = triB * (-2.0 * t);

                dN_dxi[q][6]   = gBx;
                dN_deta[q][6]  = gBy;
                dN_dzeta[q][6] = gBz;

                // Gradients des sommets avec correction
                const CfemReal c_dx = gBx * o6;
                const CfemReal c_dy = gBy * o6;
                const CfemReal c_dz = gBz * o6;

                for (CfemInt i = 0; i < 6; ++i) {
                    const CfemInt iT = i % 3;
                    const CfemInt iZ = i / 3;
                    const CfemReal drL = (iT == 0 ? -1.0 : (iT == 1 ? 1.0 : 0.0));
                    const CfemReal dsL = (iT == 0 ? -1.0 : (iT == 2 ? 1.0 : 0.0));

                    dN_dxi[q][i]   = drL * M[iZ] - c_dx;
                    dN_deta[q][i]  = dsL * M[iZ] - c_dy;
                    dN_dzeta[q][i] = L[iT] * dM[iZ] - c_dz;
                }
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
            constexpr CfemReal o6 = 1.0 / 6.0;

            for (CfemInt q = 0; q < nQp; ++q) {
                const CfemReal r = xiPoints[q].x;
                const CfemReal s = xiPoints[q].y;
                const CfemReal t = xiPoints[q].z;

                const CfemReal L[3] = {1.0 - r - s, r, s};
                const CfemReal dM[2] = {-0.5, 0.5};

                const CfemReal triB = 27.0 * L[0] * L[1] * L[2];
                const CfemReal zB = 1.0 - t * t;

                // Hessiens de la bulle (calculs simplifiés sans divisions par zéro)
                const CfemReal hBxx = -54.0 * s * zB;
                const CfemReal hByy = -54.0 * r * zB;
                const CfemReal hBzz = -2.0 * triB;
                const CfemReal hBxy =  27.0 * (1.0 - 2.0*r - 2.0*s) * zB;
                const CfemReal hBxz = -54.0 * s * t * (1.0 - 2.0*r - s); // CORRECTION: evite (gBx / zB) * (-2t)
                const CfemReal hByz = -54.0 * r * t * (1.0 - r - 2.0*s); // CORRECTION: evite (gBy / zB) * (-2t)

                d2N_dxi2[q][6]      = hBxx;
                d2N_deta2[q][6]     = hByy;
                d2N_dzeta2[q][6]    = hBzz;
                d2N_dxideta[q][6]   = hBxy;
                d2N_dxidzeta[q][6]  = hBxz;
                d2N_detadzeta[q][6] = hByz;

                // Hessiens des sommets avec correction
                const CfemReal c_xx = hBxx * o6;
                const CfemReal c_yy = hByy * o6;
                const CfemReal c_zz = hBzz * o6;
                const CfemReal c_xy = hBxy * o6;
                const CfemReal c_xz = hBxz * o6;
                const CfemReal c_yz = hByz * o6;

                for (CfemInt i = 0; i < 6; ++i) {
                    const CfemInt iT = i % 3;
                    const CfemInt iZ = i / 3;
                    const CfemReal drL = (iT == 0 ? -1.0 : (iT == 1 ? 1.0 : 0.0));
                    const CfemReal dsL = (iT == 0 ? -1.0 : (iT == 2 ? 1.0 : 0.0));

                    d2N_dxi2[q][i]      = -c_xx;
                    d2N_deta2[q][i]     = -c_yy;
                    d2N_dzeta2[q][i]    = -c_zz;
                    d2N_dxideta[q][i]   = -c_xy;
                    // La partie croisée linéaire-linéaire n'est non nulle que pour les couples impliquant z
                    d2N_dxidzeta[q][i]  = drL * dM[iZ] - c_xz;
                    d2N_detadzeta[q][i] = dsL * dM[iZ] - c_yz;
                }
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return reference_elem::utils::projectFacetRefCoordToPrismRefCoord(localFacetId, xiFacet);
        }
        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return reference_elem::utils::prism6_face_nodes[iFace];
        }
    };

    class P0ReferencePrism : public ReferenceElement {
    public:
        P0ReferencePrism() {
            m_numNodes = 1;
            m_dimension = 3;

            // Un seul nœud au centre de gravité du prisme de référence
            // Triangle base (0,0), (1,0), (0,1) -> (1/3, 1/3)
            // Segment Z de -1 à 1 -> 0
            m_nodes = { Vector3D{1.0/3.0, 1.0/3.0, 0.0} };
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
            return cfem::reference_elem::utils::projectFacetRefCoordToPrismRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::prism1_face_nodes[iFace];
        }
    };

#endif

}