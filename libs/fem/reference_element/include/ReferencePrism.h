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
    public:
        P1ReferencePrism() {
            m_numNodes = 6;
            m_dimension = 3;
            m_nodes.resize(m_numNodes);
            // Ordre VTK : Triangle inférieur Z=-1 (0,1,2), puis supérieur Z=1 (3,4,5)
            m_nodes[0] = {0,0,-1}; m_nodes[1] = {1,0,-1}; m_nodes[2] = {0,1,-1};
            m_nodes[3] = {0,0, 1}; m_nodes[4] = {1,0, 1}; m_nodes[5] = {0,1, 1};
        }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const CfemInt n = 6;
            if (out.values.size() != n) out.setup(n);
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n)
                out.hessians.assign(n, SymTensorO2());

            const CfemReal r = xi.x, s = xi.y, t = xi.z;
            
            // Fonctions Triangle P1
            CfemReal LT[3] = { 1.0 - r - s, r, s };
            CfemReal dLT_r[3] = { -1.0, 1.0, 0.0 };
            CfemReal dLT_s[3] = { -1.0, 0.0, 1.0 };

            // Fonctions Segment P1 (sur t de -1 à 1)
            CfemReal LZ[2] = { 0.5 * (1.0 - t), 0.5 * (1.0 + t) };
            CfemReal dLZ_t[2] = { -0.5, 0.5 };

            for (CfemInt i = 0; i < 6; ++i) {
                CfemInt triIdx = i % 3;      // 0, 1, 2
                CfemInt zIdx   = i / 3;      // 0 pour Z=-1, 1 pour Z=1

                if (hasFlag(requested, ShapeEvaluationFlags::Value))
                    out.values[i] = LT[triIdx] * LZ[zIdx];

                if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                    out.gradients[i].x = dLT_r[triIdx] * LZ[zIdx];
                    out.gradients[i].y = dLT_s[triIdx] * LZ[zIdx];
                    out.gradients[i].z = LT[triIdx]    * dLZ_t[zIdx];
                }

                if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                    out.hessians[i].setZero();
                    // Dérivées croisées (Plan x Hauteur)
                    out.hessians[i].xz() = dLT_r[triIdx] * dLZ_t[zIdx];
                    out.hessians[i].yz() = dLT_s[triIdx] * dLZ_t[zIdx];
                    // Les dérivées d2N/drds sont nulles ici car Triangle P1 est linéaire.
                }
            }
            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToPrismRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::prism6_face_nodes[iFace];
        }
    };


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

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const CfemInt n = 18;
            if (out.values.size() != n) out.setup(n);
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            const CfemReal r = xi.x;
            const CfemReal s = xi.y;
            const CfemReal t = xi.z;
            const CfemReal L0 = 1.0 - r - s;

            // Fonctions 2D Triangle P2
            CfemReal fT[6]   = { L0*(2*L0-1), r*(2*r-1), s*(2*s-1), 4*r*L0, 4*r*s, 4*s*L0 };
            CfemReal dfTr[6] = { -3.0+4*r+4*s, 4*r-1, 0, 4-8*r-4*s, 4*s, -4*s };
            CfemReal dfTs[6] = { -3.0+4*r+4*s, 0, 4*s-1, -4*r, 4*r, 4-4*r-8*s };
            
            // Hessiennes constantes du triangle P2 (d2N/dr2, d2N/ds2, d2N/drds)
            static const CfemReal hT[6][3] = {
                {4, 4, 4}, {4, 0, 0}, {0, 4, 0}, {-8, 0, -4}, {0, 0, 4}, {0, -8, -4}
            };

            // Fonctions 1D Segment P2 (sur t)
            CfemReal fZ[3]   = { 0.5*t*(t-1), 0.5*t*(t+1), 1-t*t }; // Bas, Haut, Milieu
            CfemReal dfZt[3] = { t-0.5, t+0.5, -2.0*t };
            CfemReal ddfZtt[3]= { 1.0, 1.0, -2.0 };

            // Mapping Topologique VTK (Index Triangle, Index Segment)
            static const CfemInt mapT[18] = {0,1,2, 0,1,2, 3,4,5, 3,4,5, 0,1,2, 3,4,5};
            static const CfemInt mapZ[18] = {0,0,0, 1,1,1, 0,0,0, 1,1,1, 2,2,2, 2,2,2};

            for (CfemInt i = 0; i < 18; ++i) {
                CfemInt iT = mapT[i];
                CfemInt iZ = mapZ[i];

                if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                    out.values[i] = fT[iT] * fZ[iZ];
                }

                if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                    out.gradients[i].x = dfTr[iT] * fZ[iZ];
                    out.gradients[i].y = dfTs[iT] * fZ[iZ];
                    out.gradients[i].z = fT[iT]   * dfZt[iZ];
                }

                if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                    out.hessians[i].setZero();
                    // d2N/dr2, d2N/ds2, d2N/drds
                    out.hessians[i].xx() = hT[iT][0] * fZ[iZ];
                    out.hessians[i].yy() = hT[iT][1] * fZ[iZ];
                    out.hessians[i].xy() = hT[iT][2] * fZ[iZ];
                    // d2N/dt2
                    out.hessians[i].zz() = fT[iT] * ddfZtt[iZ];
                    // Croisées (Plan x t)
                    out.hessians[i].xz() = dfTr[iT] * dfZt[iZ];
                    out.hessians[i].yz() = dfTs[iT] * dfZt[iZ];
                }
            }
            out.validFlags = requested;
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

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const CfemInt n = 7;
            if (out.values.size() != n) out.setup(n);
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n)
                out.hessians.assign(n, SymTensorO2());

            const CfemReal r = xi.x, s = xi.y, t = xi.z;
            const CfemReal L[3] = {1.0 - r - s, r, s};
            const CfemReal triB = 27.0 * L[0] * L[1] * L[2];
            const CfemReal zB = 1.0 - t*t;

            // 1. Calcul de la Bulle (Nœud 6)
            const CfemReal B = triB * zB;
            Vector3D gB = { 27.0 * s * (1.0 - 2.0*r - s) * zB, 
                            27.0 * r * (1.0 - r - 2.0*s) * zB, 
                            triB * (-2.0 * t) };
            
            SymTensorO2 hB; hB.setZero();
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                hB.xx() = -54.0 * s * zB;
                hB.yy() = -54.0 * r * zB;
                hB.zz() = -2.0 * triB;
                hB.xy() = 27.0 * (1.0 - 2.0*r - 2.0*s) * zB;
                hB.xz() = gB.x / zB * (-2.0 * t);
                hB.yz() = gB.y / zB * (-2.0 * t);
            }

            // 2. Fonctions des Sommets (0-5) corrigées par -B/6
            const CfemReal M[2] = {0.5*(1.0-t), 0.5*(1.0+t)};
            const CfemReal dM[2] = {-0.5, 0.5};

            for (CfemInt i = 0; i < 6; ++i) {
                CfemInt iT = i % 3, iZ = i / 3;
                CfemReal drL = (iT==0?-1:(iT==1?1:0)), dsL = (iT==0?-1:(iT==2?1:0));

                if (hasFlag(requested, ShapeEvaluationFlags::Value))
                    out.values[i] = L[iT] * M[iZ] - B / 6.0;

                if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                    out.gradients[i].x = drL * M[iZ] - gB.x / 6.0;
                    out.gradients[i].y = dsL * M[iZ] - gB.y / 6.0;
                    out.gradients[i].z = L[iT] * dM[iZ] - gB.z / 6.0;
                }

                if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                    out.hessians[i].setZero();
                    out.hessians[i].xz() = drL * dM[iZ] - hB.xz() / 6.0;
                    out.hessians[i].yz() = dsL * dM[iZ] - hB.yz() / 6.0;
                    out.hessians[i].xx() = -hB.xx() / 6.0;
                    out.hessians[i].yy() = -hB.yy() / 6.0;
                    out.hessians[i].zz() = -hB.zz() / 6.0;
                    out.hessians[i].xy() = -hB.xy() / 6.0;
                }
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Value))    out.values[6] = B;
            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) out.gradients[6] = gB;
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian))  out.hessians[6] = hB;

            out.validFlags = requested;
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
            m_nodes.resize(m_numNodes);
            // Un seul nœud au centre de gravité du prisme de référence
            // Triangle base (0,0), (1,0), (0,1) -> (1/3, 1/3)
            // Segment Z de -1 à 1 -> 0
            m_nodes = { {1.0/3.0, 1.0/3.0, 0.0} };
        }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const CfemInt n = 1;
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

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToPrismRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::prism1_face_nodes[iFace];
        }
    };

#endif

}