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
}

namespace cfem{

class P1ReferencePrism : public ReferenceElement {
public:
    P1ReferencePrism() {
        m_nodes.resize(6);
        // Ordre VTK : Triangle inférieur Z=-1 (0,1,2), puis supérieur Z=1 (3,4,5)
        m_nodes[0] = {0,0,-1}; m_nodes[1] = {1,0,-1}; m_nodes[2] = {0,1,-1};
        m_nodes[3] = {0,0, 1}; m_nodes[4] = {1,0, 1}; m_nodes[5] = {0,1, 1};
    }

    int dimension() const override { return 3; }
    int getNumNodes() const override { return 6; }

    void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
        const int n = 6;
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

        for (int i = 0; i < 6; ++i) {
            int triIdx = i % 3;      // 0, 1, 2
            int zIdx   = i / 3;      // 0 pour Z=-1, 1 pour Z=1

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

    Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
        return cfem::reference_elem::utils::projectFacetRefCoordToPrismRefCoord(localFacetId, xiFacet); // Make the compiler happy
    }
};


    class P2ReferencePrism : public ReferenceElement {
    public:
        P2ReferencePrism() {
            m_nodes.resize(18);
            
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

        int dimension() const override { return 3; }
        int getNumNodes() const override { return 18; }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 18;
            if (out.values.size() != n) out.setup(n);
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            const CfemReal r = xi.x;
            const CfemReal s = xi.y;
            const CfemReal t = xi.z;
            const CfemReal L0 = 1.0 - r - s;

            // --- Fonctions 2D Triangle P2 ---
            CfemReal fT[6]   = { L0*(2*L0-1), r*(2*r-1), s*(2*s-1), 4*r*L0, 4*r*s, 4*s*L0 };
            CfemReal dfTr[6] = { -3.0+4*r+4*s, 4*r-1, 0, 4-8*r-4*s, 4*s, -4*s };
            CfemReal dfTs[6] = { -3.0+4*r+4*s, 0, 4*s-1, -4*r, 4*r, 4-4*r-8*s };
            
            // Hessiennes constantes du triangle P2 (d2N/dr2, d2N/ds2, d2N/drds)
            static const CfemReal hT[6][3] = {
                {4, 4, 4}, {4, 0, 0}, {0, 4, 0}, {-8, 0, -4}, {0, 0, 4}, {0, -8, -4}
            };

            // --- Fonctions 1D Segment P2 (sur t) ---
            CfemReal fZ[3]   = { 0.5*t*(t-1), 0.5*t*(t+1), 1-t*t }; // Bas, Haut, Milieu
            CfemReal dfZt[3] = { t-0.5, t+0.5, -2.0*t };
            CfemReal ddfZtt[3]= { 1.0, 1.0, -2.0 };

            // --- Mapping Topologique VTK (Index Triangle, Index Segment) ---
            static const int mapT[18] = {0,1,2, 0,1,2, 3,4,5, 3,4,5, 0,1,2, 3,4,5};
            static const int mapZ[18] = {0,0,0, 1,1,1, 0,0,0, 1,1,1, 2,2,2, 2,2,2};

            for (int i = 0; i < 18; ++i) {
                int iT = mapT[i];
                int iZ = mapZ[i];

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

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToPrismRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }
    };

#endif

}