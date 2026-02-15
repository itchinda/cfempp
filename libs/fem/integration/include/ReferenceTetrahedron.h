#ifndef CFEM_REFERENCE_TETRAHEDRON
#define CFEM_REFERENCE_TETRAHEDRON

#include "ReferenceElement.h"

namespace cfem {

    /**
     * @brief Tétraèdre Linéaire (4 nœuds)
     */
    class P1ReferenceTetrahedron : public ReferenceElement {
    public:
        P1ReferenceTetrahedron() {
            m_nodes = { {0,0,0}, {1,0,0}, {0,1,0}, {0,0,1} };
        }

        int dimension() const override { return 3; }
        int getNumNodes() const override { return 4; }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 4;
            if (out.values.size() != n) out.setup(n);
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;
            const CfemReal L3 = xi.z;
            const CfemReal L0 = 1.0 - L1 - L2 - L3;

            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = L0; out.values[1] = L1;
                out.values[2] = L2; out.values[3] = L3;
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                out.gradients[0] = {-1, -1, -1};
                out.gradients[1] = { 1,  0,  0};
                out.gradients[2] = { 0,  1,  0};
                out.gradients[3] = { 0,  0,  1};
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                for (auto& H : out.hessians) H.setZero();
            }
            out.validFlags = requested;
        }
    };

    /**
     * @brief Tétraèdre Quadratique (10 nœuds)
     */
    class P2ReferenceTetrahedron : public ReferenceElement {
    public:
        P2ReferenceTetrahedron() {
            m_nodes.resize(10);
            // Sommets
            m_nodes[0]={0,0,0}; m_nodes[1]={1,0,0}; m_nodes[2]={0,1,0}; m_nodes[3]={0,0,1};
            // Milieux d'arêtes (Standard convention: 01, 12, 20, 03, 13, 23)
            m_nodes[4]={0.5, 0.0, 0.0}; m_nodes[5]={0.5, 0.5, 0.0}; m_nodes[6]={0.0, 0.5, 0.0};
            m_nodes[7]={0.0, 0.0, 0.5}; m_nodes[8]={0.5, 0.0, 0.5}; m_nodes[9]={0.0, 0.5, 0.5};
        }

        int dimension() const override { return 3; }
        int getNumNodes() const override { return 10; }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 10;
            if (out.values.size() != n) out.setup(n);
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;
            const CfemReal L3 = xi.z;
            const CfemReal L0 = 1.0 - L1 - L2 - L3;

            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = L0*(2*L0-1); out.values[1] = L1*(2*L1-1);
                out.values[2] = L2*(2*L2-1); out.values[3] = L3*(2*L3-1);
                out.values[4] = 4*L0*L1;     out.values[5] = 4*L1*L2;
                out.values[6] = 4*L2*L0;     out.values[7] = 4*L0*L3;
                out.values[8] = 4*L1*L3;     out.values[9] = 4*L2*L3;
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                // Exemple pour N0: Grad(L0(2L0-1)) = (4L0-1)*Grad(L0) = (4L0-1)*(-1,-1,-1)
                auto gL0 = Vector3D{-1,-1,-1};
                out.gradients[0] = gL0 * (4*L0-1);
                out.gradients[1] = {4*L1-1, 0, 0};
                out.gradients[2] = {0, 4*L2-1, 0};
                out.gradients[3] = {0, 0, 4*L3-1};
                // Milieux (ex N4 = 4L0L1 -> Grad = 4(L1*GradL0 + L0*GradL1))
                out.gradients[4] = {4*(L0-L1), -4*L1, -4*L1};
                out.gradients[5] = {4*L2, 4*L1, 0};
                out.gradients[6] = {-4*L2, 4*(L0-L2), -4*L2};
                out.gradients[7] = {-4*L3, -4*L3, 4*(L0-L3)};
                out.gradients[8] = {4*L3, 0, 4*L1};
                out.gradients[9] = {0, 4*L3, 4*L2};
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                for (int i = 0; i < n; ++i) out.hessians[i].setZero();

                // --- Sommets ---
                // N0 = L0(2L0-1) avec L0 = 1-x-y-z. 
                // d2N0/dx2 = d2N0/dy2 = d2N0/dz2 = 4 | d2N0/dxdy = d2N0/dxdz = d2N0/dydz = 4
                out.hessians[0].setTo(4.0);
                out.hessians[1].xx() = 4.0; // Ni = Li(2Li-1)
                out.hessians[2].yy() = 4.0;
                out.hessians[3].zz() = 4.0;

                // --- Milieux d'arêtes ---
                // N4 = 4*L0*L1 = 4*(1-x-y-z)*x = 4x - 4x^2 - 4xy - 4xz
                out.hessians[4].xx() = -8.0; 
                out.hessians[4].xy() = -4.0;
                out.hessians[4].xz() = -4.0;

                // N5 = 4*L1*L2 = 4*x*y
                out.hessians[5].xy() =  4.0;

                // N6 = 4*L2*L0 = 4*y*(1-x-y-z) = 4y - 4xy - 4y^2 - 4yz
                out.hessians[6].yy() = -8.0;
                out.hessians[6].xy() = -4.0;
                out.hessians[6].yz() = -4.0;

                // N7 = 4*L0*L3 = 4*(1-x-y-z)*z = 4z - 4xz - 4yz - 4z^2
                out.hessians[7].zz() = -8.0;
                out.hessians[7].xz() = -4.0;
                out.hessians[7].yz() = -4.0;

                // N8 = 4*L1*L3 = 4*x*z
                out.hessians[8].xz() =  4.0;

                // N9 = 4*L2*L3 = 4*y*z
                out.hessians[9].yz() =  4.0;
            }
            out.validFlags = requested;
        }
    };
}
#endif