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
        {0, 1, 2}, 
        {0, 1, 3}, 
        {1, 2, 3}, 
        {2, 0, 3}
    };

    static const std::vector<std::vector<int>> tet10_face_nodes = {
        {0, 1, 2, 4, 5, 6}, 
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

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToTetraRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::tet4_face_nodes[iFace];
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

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToTetraRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::tet10_face_nodes[iFace];
        }
    };

class P0ReferenceTetrahedron : public ReferenceElement {
    public:
        P0ReferenceTetrahedron() {
            // Un seul nœud au centre de gravité du tétraèdre de référence
            // Sommets: (0,0,0), (1,0,0), (0,1,0), (0,0,1) -> Barycentre: (1/4, 1/4, 1/4)
            m_nodes = { {0.25, 0.25, 0.25} };
        }

        int dimension() const override { return 3; }
        int getNumNodes() const override { return 1; }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 1;
            if (out.values.size() != n) out.setup(n);

            // Valeur constante N0 = 1.0
            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = 1.0;
            }

            // Gradient nul dans les 3 directions
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
            return cfem::reference_elem::utils::projectFacetRefCoordToTetraRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::tet1_face_nodes[iFace];
        }
    };
}
#endif