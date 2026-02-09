#ifndef CFEM_REFERENCE_QUAD
#define CFEM_REFERENCE_QUAD

#include "ReferenceElement.h"

namespace cfem {

    class Q1ReferenceQuad : public ReferenceElement {
    public:
        Q1ReferenceQuad() {
            m_nodes.resize(4);
            m_nodes[0] = Vector3D{-1.0, -1.0, 0.0};
            m_nodes[1] = Vector3D{ 1.0, -1.0, 0.0};
            m_nodes[2] = Vector3D{ 1.0,  1.0, 0.0};
            m_nodes[3] = Vector3D{-1.0,  1.0, 0.0};
        }

        int dimension() const override { return 2; }
        int numNodes() const override { return 4; }

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
    };


class Q2ReferenceQuad : public ReferenceElement {
    public:
        Q2ReferenceQuad() {
            m_nodes.resize(9);
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

        int dimension() const override { return 2; }
        int numNodes() const override { return 9; }

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
    };
}
#endif