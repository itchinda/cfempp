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

namespace cfem::reference_elem::utils{

    inline Vector3D projectFacetRefCoordToQuadRefCoord(CfemInt localFacetId, 
                                                       const Vector3D& xiFacet) 
    {
        const CfemReal u = xiFacet.x;

        switch (localFacetId) {
            case 0: return Vector3D(u, -1.0, 0.0);             // Face {0,1} : y=-1
            case 1: return Vector3D(1.0, u, 0.0);              // Face {1,2} : x=1
            case 2: return Vector3D(-u, 1.0, 0.0);             // Face {2,3} : y=1
            case 3: return Vector3D(-1.0, -u, 0.0);            // Face {3,0} : x=-1
            default: CFEM_ERROR << "Invalid facet ID for Quad cell";
        }
        return Vector3D{}; // Make the compiler happy
    }

    static const std::vector<std::vector<int>> quad4_face_nodes = {
    {0, 1}, {1, 2}, {2, 3}, {3, 0}
    };
    static const std::vector<std::vector<int>> quad9_face_nodes = {
        {0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7}
    };
    static const std::vector<std::vector<int>> quad1_face_nodes = {
        {0}, {0}, {0}, {0}
    };
}

namespace cfem {

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

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 4;
            if (out.values.size() != n) out.setup(n);
            
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            const CfemReal x = xi.x;
            const CfemReal y = xi.y;

            // Facteurs de signes pour les 4 nœuds
            const CfemReal sx[4] = {-1, 1, 1, -1};
            const CfemReal sy[4] = {-1, -1, 1, 1};

            for (int i = 0; i < n; ++i) {
                if ( hasFlag(requested, ShapeEvaluationFlags::Value)) {
                    out.values[i] = 0.25 * (1.0 + sx[i] * x) * (1.0 + sy[i] * y);
                }

                if ( hasFlag(requested, ShapeEvaluationFlags::Gradient) ) {
                    out.gradients[i].x = 0.25 * sx[i] * (1.0 + sy[i] * y);
                    out.gradients[i].y = 0.25 * sy[i] * (1.0 + sx[i] * x);
                    out.gradients[i].z = 0.0;
                }

                if ( hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                    out.hessians[i].setZero();
                    // Seule la dérivée croisée est non nulle : d2N / dx dy
                    CfemReal d2N_dxy = 0.25 * sx[i] * sy[i];
                    out.hessians[i].xy() = d2N_dxy;
                }
            }
            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToQuadRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::quad4_face_nodes[iFace];
        }
    };


class Q2ReferenceQuad : public ReferenceElement {
    public:
        Q2ReferenceQuad() {
            m_numNodes = 9;
            m_dimension = 2;
            m_nodes.resize(m_numNodes);
            // Coordonnées de référence sur [-1, 1]
            // Sommets
            m_nodes[0] = Vector3D{-1.0, -1.0, 0.0}; // i=0, j=0
            m_nodes[1] = Vector3D{ 1.0, -1.0, 0.0}; // i=2, j=0
            m_nodes[2] = Vector3D{ 1.0,  1.0, 0.0}; // i=2, j=2
            m_nodes[3] = Vector3D{-1.0,  1.0, 0.0}; // i=0, j=2

            // Milieux d'arêtes
            m_nodes[4] = Vector3D{ 0.0, -1.0, 0.0}; // i=1, j=0
            m_nodes[5] = Vector3D{ 1.0,  0.0, 0.0}; // i=2, j=1
            m_nodes[6] = Vector3D{ 0.0,  1.0, 0.0}; // i=1, j=2
            m_nodes[7] = Vector3D{-1.0,  0.0, 0.0}; // i=0, j=1

            // Centre
            m_nodes[8] = Vector3D{ 0.0,  0.0, 0.0}; // i=1, j=1
        }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 9;
            if (out.values.size() != n) out.setup(n);
            
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            const CfemReal r = xi.x;
            const CfemReal s = xi.y;

            // Fonctions de base Lagrange 1D (P2)
            CfemReal fr[3] = { 0.5 * r * (r - 1.0), 1.0 - r * r, 0.5 * r * (r + 1.0) };
            CfemReal fs[3] = { 0.5 * s * (s - 1.0), 1.0 - s * s, 0.5 * s * (s + 1.0) };
            
            // Dérivées premières 1D
            CfemReal dfr[3] = { r - 0.5, -2.0 * r, r + 0.5 };
            CfemReal dfs[3] = { s - 0.5, -2.0 * s, s + 0.5 };

            // Dérivées secondes 1D
            CfemReal ddfr[3] = { 1.0, -2.0, 1.0 };
            CfemReal ddfs[3] = { 1.0, -2.0, 1.0 };

            // Mapping (i,j) -> index global m_nodes
            // i=r index (0:-1, 1:0, 2:1), j=s index (0:-1, 1:0, 2:1)
            static const int map[3][3] = {
                {0, 7, 3}, // r = -1.0
                {4, 8, 6}, // r =  0.0
                {1, 5, 2}  // r =  1.0
            };

            for (int i = 0; i < 3; ++i) {
                for (int j = 0; j < 3; ++j) {
                    int idx = map[i][j];

                    if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                        out.values[idx] = fr[i] * fs[j];
                    }

                    if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                        out.gradients[idx].x = dfr[i] * fs[j];
                        out.gradients[idx].y = fr[i] * dfs[j];
                        out.gradients[idx].z = 0.0;
                    }

                    if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                        out.hessians[idx].setZero();
                        out.hessians[idx].xx() = ddfr[i] * fs[j];   // d2N/dr2
                        out.hessians[idx].yy() = fr[i] * ddfs[j];   // d2N/ds2
                        out.hessians[idx].xy() = dfr[i] * dfs[j];   // d2N/drds
                    }
                }
            }
            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToQuadRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::quad9_face_nodes[iFace];
        }
    };

    class MiniReferenceQuad : public ReferenceElement {
    public:
        MiniReferenceQuad() {
            m_numNodes = 5;
            m_dimension = 2;
            m_nodes.resize(5);
            m_nodes[0] = {-1,-1,0}; m_nodes[1] = {1,-1,0};
            m_nodes[2] = { 1, 1, 0}; m_nodes[3] = {-1, 1,0};
            m_nodes[4] = { 0, 0, 0}; // Bulle
        }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 5;
            if (out.values.size() != n) out.setup(n);
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n)
                out.hessians.assign(n, SymTensorO2());

            const CfemReal x = xi.x, y = xi.y;
            
            // 1. Calcul de la Bulle (Nœud 4)
            const CfemReal bx = 1.0 - x*x, by = 1.0 - y*y;
            const CfemReal B = bx * by;
            const Vector3D gradB = {-2.0*x*by, -2.0*y*bx, 0.0};
            SymTensorO2 hessB;
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                hessB.setZero();
                hessB.xx() = -2.0 * by;
                hessB.yy() = -2.0 * bx;
                hessB.xy() = 4.0 * x * y;
            }

            // 2. Fonctions des Sommets (0-3) corrigées par -B/4
            const CfemReal sx[4] = {-1, 1, 1, -1}, sy[4] = {-1, -1, 1, 1};
            for (int i = 0; i < 4; ++i) {
                if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                    CfemReal Q1 = 0.25 * (1.0 + sx[i]*x) * (1.0 + sy[i]*y);
                    out.values[i] = Q1 - 0.25 * B;
                }
                if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                    out.gradients[i].x = 0.25 * sx[i] * (1.0 + sy[i]*y) - 0.25 * gradB.x;
                    out.gradients[i].y = 0.25 * sy[i] * (1.0 + sx[i]*x) - 0.25 * gradB.y;
                }
                if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                    out.hessians[i].setZero();
                    out.hessians[i].xy() = 0.25 * sx[i] * sy[i] - 0.25 * hessB.xy();
                    out.hessians[i].xx() = -0.25 * hessB.xx();
                    out.hessians[i].yy() = -0.25 * hessB.yy();
                }
            }

            // 3. Stockage final de la Bulle
            if (hasFlag(requested, ShapeEvaluationFlags::Value))    out.values[4] = B;
            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) out.gradients[4] = gradB;
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian))  out.hessians[4] = hessB;

            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override {
            return reference_elem::utils::projectFacetRefCoordToQuadRefCoord(localFacetId, xiFacet);
        }
        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return reference_elem::utils::quad4_face_nodes[iFace];
        }
    };

class Q0ReferenceQuad : public ReferenceElement {
    public:
        Q0ReferenceQuad() {
            m_numNodes = 1;
            m_dimension = 2;
            m_nodes.resize(m_numNodes);
            // Un seul nœud au centre du carré de référence [-1, 1]x[-1, 1]
            m_nodes = { {0.0, 0.0, 0.0} };
        }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 1;
            if (out.values.size() != n) out.setup(n);

            // Valeur constante N0 = 1.0
            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = 1.0;
            }

            // Gradient nul
            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                out.gradients[0] = Vector3D{0.0, 0.0, 0.0};
            }

            // Hessien nul
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                if (out.hessians.size() != n) out.hessians.assign(n, SymTensorO2());
                out.hessians[0].setZero();
            }

            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToQuadRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::quad1_face_nodes[iFace];
        }
    };
}
#endif