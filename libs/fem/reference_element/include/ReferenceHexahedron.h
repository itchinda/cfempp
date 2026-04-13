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

namespace cfem::reference_elem::utils{

    inline Vector3D projectFacetRefCoordToHexaRefCoord(CfemInt localFacetId, 
                                                      const Vector3D& xiFacet) 
    {
        const CfemReal u = xiFacet.x;
        const CfemReal v = xiFacet.y;

        switch (localFacetId) {
            case 0: return Vector3D(v, u, -1.0);               // Face {0,3,2,1} : z=-1
            case 1: return Vector3D(u, v, 1.0);                // Face {4,5,6,7} : z=1
            case 2: return Vector3D(u, -1.0, v);               // Face {0,1,5,4} : y=-1
            case 3: return Vector3D(1.0, u, v);                // Face {1,2,6,5} : x=1
            case 4: return Vector3D(-u, 1.0, v);               // Face {2,3,7,6} : y=1
            case 5: return Vector3D(-1.0, -u, v);              // Face {3,0,4,7} : x=-1
            default: CFEM_ERROR << "Invalid facet ID for Hexa";
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
        {0, 3, 2, 1, 11, 10, 9, 8, 24},          // Face 0
        {4, 5, 6, 7, 12, 13, 14, 15, 25},        // Face 1
        {0, 1, 5, 4, 8, 17, 12, 16, 22},         // Face 2
        {1, 2, 6, 5, 9, 18, 13, 17, 21},         // Face 3
        {2, 3, 7, 6, 10, 19, 14, 18, 23},        // Face 4
        {3, 0, 4, 7, 11, 16, 15, 19, 20}         // Face 5 
    };

    static const std::vector<std::vector<int>> hex1_face_nodes = {
        {0}, {0}, {0}, {0}, {0}, {0}
    };
}

namespace cfem{

    class Q1ReferenceHexahedron : public ReferenceElement {
    public:
        Q1ReferenceHexahedron() {
            m_numNodes = 8;
            m_dimension = 3;
            m_nodes.resize(m_numNodes);
            // Ordre VTK : Face Z=-1 (0,1,2,3), puis Face Z=1 (4,5,6,7)
            m_nodes[0] = {-1,-1,-1}; m_nodes[1] = { 1,-1,-1}; m_nodes[2] = { 1, 1,-1}; m_nodes[3] = {-1, 1,-1};
            m_nodes[4] = {-1,-1, 1}; m_nodes[5] = { 1,-1, 1}; m_nodes[6] = { 1, 1, 1}; m_nodes[7] = {-1, 1, 1};
        }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = m_numNodes;
            if (out.values.size() != n) out.setup(n);
            if ((hasFlag(requested, ShapeEvaluationFlags::Hessian)) && out.hessians.size() != n)
                out.hessians.assign(n, SymTensorO2());

            const CfemReal r = xi.x, s = xi.y, t = xi.z;
            const CfemReal xr[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
            const CfemReal xs[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
            const CfemReal xt[8] = {-1, -1, -1, -1, 1, 1, 1, 1};

            for (int i = 0; i < n; ++i) {
                CfemReal pr = 1.0 + xr[i]*r;
                CfemReal ps = 1.0 + xs[i]*s;
                CfemReal pt = 1.0 + xt[i]*t;

                if (hasFlag(requested, ShapeEvaluationFlags::Value)) 
                    out.values[i] = 0.125 * pr * ps * pt;

                if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                    out.gradients[i].x = 0.125 * xr[i] * ps * pt;
                    out.gradients[i].y = 0.125 * xs[i] * pr * pt;
                    out.gradients[i].z = 0.125 * xt[i] * pr * ps;
                }

                if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                    out.hessians[i].setZero();
                    out.hessians[i].xy() = 0.125 * xr[i] * xs[i] * pt;
                    out.hessians[i].xz() = 0.125 * xr[i] * ps * xt[i];
                    out.hessians[i].yz() = 0.125 * pr * xs[i] * xt[i];
                }
            }
            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToHexaRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::hex8_face_nodes[iFace];
        }
    };


    class Q2ReferenceHexahedron : public ReferenceElement {
    public:
        Q2ReferenceHexahedron() {
            m_numNodes = 27;
            m_dimension = 3;
            m_nodes.resize(m_numNodes);
            // On définit les positions 1D pour r, s, t : -1.0, 0.0, 1.0
            CfemReal pos[3] = {-1.0, 0.0, 1.0};

            // Mapping VTK pour HEX27 (i, j, k correspondent aux indices de pos[3])
            // Cet ordre est crucial pour passer le test de Kronecker
            int vtk_map[3][3][3];
            // Sommets
            vtk_map[0][0][0] = 0; vtk_map[2][0][0] = 1; vtk_map[2][2][0] = 2; vtk_map[0][2][0] = 3;
            vtk_map[0][0][2] = 4; vtk_map[2][0][2] = 5; vtk_map[2][2][2] = 6; vtk_map[0][2][2] = 7;
            // Milieux d'arêtes
            vtk_map[1][0][0] = 8;  vtk_map[2][1][0] = 9;  vtk_map[1][2][0] = 10; vtk_map[0][1][0] = 11;
            vtk_map[1][0][2] = 12; vtk_map[2][1][2] = 13; vtk_map[1][2][2] = 14; vtk_map[0][1][2] = 15;
            vtk_map[0][0][1] = 16; vtk_map[2][0][1] = 17; vtk_map[2][2][1] = 18; vtk_map[0][2][1] = 19;
            // Centres de faces
            vtk_map[0][1][1] = 20; vtk_map[2][1][1] = 21; vtk_map[1][0][1] = 22;
            vtk_map[1][2][1] = 23; vtk_map[1][1][0] = 24; vtk_map[1][1][2] = 25;
            // Centre volume
            vtk_map[1][1][1] = 26;

            for (int i = 0; i < 3; ++i)
                for (int j = 0; j < 3; ++j)
                    for (int k = 0; k < 3; ++k)
                        m_nodes[vtk_map[i][j][k]] = {pos[i], pos[j], pos[k]};
        }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = m_numNodes;
            if (out.values.size() != n) out.setup(n);
            if ((hasFlag(requested, ShapeEvaluationFlags::Hessian)) && out.hessians.size() != n)
                out.hessians.assign(n, SymTensorO2());

            const CfemReal r = xi.x, s = xi.y, t = xi.z;

            // Fonctions 1D P2 et leurs dérivées
            CfemReal fr[3] = {0.5*r*(r-1.0), 1.0-r*r, 0.5*r*(r+1.0)};
            CfemReal fs[3] = {0.5*s*(s-1.0), 1.0-s*s, 0.5*s*(s+1.0)};
            CfemReal ft[3] = {0.5*t*(t-1.0), 1.0-t*t, 0.5*t*(t+1.0)};

            CfemReal dfr[3] = {r-0.5, -2.0*r, r+0.5};
            CfemReal dfs[3] = {s-0.5, -2.0*s, s+0.5};
            CfemReal dft[3] = {t-0.5, -2.0*t, t+0.5};

            CfemReal ddfr[3] = {1.0, -2.0, 1.0};
            CfemReal ddfs[3] = {1.0, -2.0, 1.0};
            CfemReal ddft[3] = {1.0, -2.0, 1.0};

            // Re-déclaration locale du map pour la boucle (ou utilisation d'un membre static)
            static const int map[3][3][3] = {
                {{0, 16, 4}, {11, 20, 15}, {3, 19, 7}},
                {{8, 22, 12}, {24, 26, 25}, {10, 23, 14}},
                {{1, 17, 5}, {9, 21, 13}, {2, 18, 6}}
            };

            for (int i = 0; i < 3; ++i) { // r
                for (int j = 0; j < 3; ++j) { // s
                    for (int k = 0; k < 3; ++k) { // t
                        int idx = map[i][j][k];

                        // if (hasFlag(requested, ShapeEvaluationFlags::Value))
                        //     out.values[idx] = fr[i] * fs[j] * ft[idx < 0 ? 0 : k]; // Correction logique simple

                        // Plus proprement :
                        CfemReal val_i = fr[i], val_j = fs[j], val_k = ft[k];
                        
                        if (hasFlag(requested, ShapeEvaluationFlags::Value))
                            out.values[idx] = val_i * val_j * val_k;

                        if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                            out.gradients[idx].x = dfr[i] * val_j * val_k;
                            out.gradients[idx].y = val_i * dfs[j] * val_k;
                            out.gradients[idx].z = val_i * val_j * dft[k];
                        }

                        if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                            out.hessians[idx].setZero();
                            // Diagonales
                            out.hessians[idx].xx() = ddfr[i] * val_j * val_k;
                            out.hessians[idx].yy() = val_i * ddfs[j] * val_k;
                            out.hessians[idx].zz() = val_i * val_j * ddft[k];
                            // Croisées
                            out.hessians[idx].xy() = dfr[i] * dfs[j] * val_k;
                            out.hessians[idx].xz() = dfr[i] * val_j * dft[k];
                            out.hessians[idx].yz() = val_i * dfs[j] * dft[k];
                        }
                    }
                }
            }
            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToHexaRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::hex27_face_nodes[iFace];
        }
    };

    class MiniReferenceHexahedron : public ReferenceElement {
    public:
        MiniReferenceHexahedron() {
            m_numNodes = 9; // 8 sommets + 1 bulle
            m_dimension = 3;
            m_nodes.resize(m_numNodes);
            
            // Ordre VTK Q1
            m_nodes[0] = {-1,-1,-1}; m_nodes[1] = { 1,-1,-1}; m_nodes[2] = { 1, 1,-1}; m_nodes[3] = {-1, 1,-1};
            m_nodes[4] = {-1,-1, 1}; m_nodes[5] = { 1,-1, 1}; m_nodes[6] = { 1, 1, 1}; m_nodes[7] = {-1, 1, 1};
            // Nœud Bulle au centre
            m_nodes[8] = {0, 0, 0};
        }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = m_numNodes;
            if (out.values.size() != n) out.setup(n);
            if ((hasFlag(requested, ShapeEvaluationFlags::Hessian)) && out.hessians.size() != n)
                out.hessians.assign(n, SymTensorO2());

            const CfemReal x = xi.x, y = xi.y, z = xi.z;
            const CfemReal xr[8] = {-1, 1, 1, -1, -1, 1, 1, -1};
            const CfemReal ys[8] = {-1, -1, 1, 1, -1, -1, 1, 1};
            const CfemReal zt[8] = {-1, -1, -1, -1, 1, 1, 1, 1};

            // --- Calcul de la Bulle ---
            CfemReal b = (1.0 - x*x) * (1.0 - y*y) * (1.0 - z*z);
            CfemReal db_dx = -2.0 * x * (1.0 - y*y) * (1.0 - z*z);
            CfemReal db_dy = -2.0 * y * (1.0 - x*x) * (1.0 - z*z);
            CfemReal db_dz = -2.0 * z * (1.0 - x*x) * (1.0 - y*y);

            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[8] = b;
            }
            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                out.gradients[8] = Vector3D{db_dx, db_dy, db_dz};
            }
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                out.hessians[8].setZero();
                out.hessians[8].xx() = -2.0 * (1.0 - y*y) * (1.0 - z*z);
                out.hessians[8].yy() = -2.0 * (1.0 - x*x) * (1.0 - z*z);
                out.hessians[8].zz() = -2.0 * (1.0 - x*x) * (1.0 - y*y);
                out.hessians[8].xy() =  4.0 * x * y * (1.0 - z*z);
                out.hessians[8].xz() =  4.0 * x * z * (1.0 - y*y);
                out.hessians[8].yz() =  4.0 * y * z * (1.0 - x*x);
            }

            // --- Calcul des fonctions Q1 corrigées ---
            for (int i = 0; i < 8; ++i) {
                CfemReal px = 1.0 + xr[i]*x;
                CfemReal py = 1.0 + ys[i]*y;
                CfemReal pz = 1.0 + zt[i]*z;

                if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                    out.values[i] = (0.125 * px * py * pz) - (b / 8.0);
                }

                if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                    out.gradients[i].x = (0.125 * xr[i] * py * pz) - (db_dx / 8.0);
                    out.gradients[i].y = (0.125 * ys[i] * px * pz) - (db_dy / 8.0);
                    out.gradients[i].z = (0.125 * zt[i] * px * py) - (db_dz / 8.0);
                }

                if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                    out.hessians[i].setZero();
                    out.hessians[i].xx() = -out.hessians[8].xx() / 8.0;
                    out.hessians[i].yy() = -out.hessians[8].yy() / 8.0;
                    out.hessians[i].zz() = -out.hessians[8].zz() / 8.0;
                    // Les croisées ont une partie Q1 + la correction bulle
                    out.hessians[i].xy() = (0.125 * xr[i] * ys[i] * pz) - (out.hessians[8].xy() / 8.0);
                    out.hessians[i].xz() = (0.125 * xr[i] * py * zt[i]) - (out.hessians[8].xz() / 8.0);
                    out.hessians[i].yz() = (0.125 * px * ys[i] * zt[i]) - (out.hessians[8].yz() / 8.0);
                }
            }
            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToHexaRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::hex8_face_nodes[iFace];
        }
    };

    class Q0ReferenceHexahedron : public ReferenceElement {
    public:
        Q0ReferenceHexahedron() {
            m_numNodes = 1;
            m_dimension = 3;
            m_nodes.resize(m_numNodes);
            // Un seul nœud au centre du cube de référence [-1, 1]^3
            m_nodes = { {0.0, 0.0, 0.0} };
        }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = m_numNodes;
            if (out.values.size() != n) out.setup(n);

            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = 1.0;
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                out.gradients[0] = Vector3D{0.0, 0.0, 0.0};
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                if (out.hessians.size() != n) out.hessians.assign(n, SymTensorO2());
                out.hessians[0].setZero();
            }

            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToHexaRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::hex1_face_nodes[iFace];
        }
    };


}



#endif