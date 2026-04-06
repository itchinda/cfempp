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
}

namespace cfem {
    class P1ReferenceTriangle : public ReferenceElement {

    public:
        
        P1ReferenceTriangle() {
            m_nodes = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} };
        }

        int dimension() const override { return 2; }
        int getNumNodes() const override { return 3; }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = getNumNodes();
    
            if (out.values.size() != n) {
                out.setup(n); 
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            // N0 = 1 - xi - eta, N1 = xi, N2 = eta
            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = 1.0 - xi.x - xi.y;
                out.values[1] = xi.x;
                out.values[2] = xi.y;
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                out.gradients[0] = Vector3D{-1.0, -1.0, 0.0};
                out.gradients[1] = Vector3D{ 1.0,  0.0, 0.0};
                out.gradients[2] = Vector3D{ 0.0,  1.0, 0.0};
            }

            // Hessiens par rapport à (xi, eta)
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                for (auto& H : out.hessians) H.setZero();
            }

            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

    };


    class P2ReferenceTriangle : public ReferenceElement {
    public:
        P2ReferenceTriangle() {
            m_nodes.resize(6);
            // Main nodes
            m_nodes[0] = Vector3D{0.0, 0.0, 0.0};
            m_nodes[1] = Vector3D{1.0, 0.0, 0.0};
            m_nodes[2] = Vector3D{0.0, 1.0, 0.0};
            // Edges Middles
            m_nodes[3] = Vector3D{0.5, 0.0, 0.0};
            m_nodes[4] = Vector3D{0.5, 0.5, 0.0};
            m_nodes[5] = Vector3D{0.0, 0.5, 0.0};
        }

        int dimension() const override { return 2; }
        int getNumNodes() const override { return 6; }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 6;
            if (out.values.size() != n) out.setup(n);
            
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;
            const CfemReal L0 = 1.0 - L1 - L2;

            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = L0 * (2.0 * L0 - 1.0);
                out.values[1] = L1 * (2.0 * L1 - 1.0);
                out.values[2] = L2 * (2.0 * L2 - 1.0);
                out.values[3] = 4.0 * L0 * L1;
                out.values[4] = 4.0 * L1 * L2;
                out.values[5] = 4.0 * L2 * L0;
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                out.gradients[0] = Vector3D{-3.0 + 4.0*L1 + 4.0*L2, -3.0 + 4.0*L1 + 4.0*L2, 0.0};
                out.gradients[1] = Vector3D{4.0*L1 - 1.0, 0.0, 0.0};
                out.gradients[2] = Vector3D{0.0, 4.0*L2 - 1.0, 0.0};
                out.gradients[3] = Vector3D{4.0 - 8.0*L1 - 4.0*L2, -4.0*L1, 0.0};
                out.gradients[4] = Vector3D{4.0*L2, 4.0*L1, 0.0};
                out.gradients[5] = Vector3D{-4.0*L2, 4.0 - 4.0*L1 - 8.0*L2, 0.0};
            }

            // --- Hessiennes (Dérivées secondes constantes) ---
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                for (int i = 0; i < n; ++i) {
                    out.hessians[i].setZero(); 
                }

                out.hessians[0].xx() = 4.0; out.hessians[0].yy() = 4.0; out.hessians[0].xy()=4.0;
                out.hessians[1].xx() = 4.0;
                out.hessians[2].yy() = 4.0;
                out.hessians[3].xx() =-8.0; out.hessians[3].xy() = -4.0;
                out.hessians[4].xy() = 4.0; out.hessians[4].xy() =  4.0;
                out.hessians[5].yy() =-8.0; out.hessians[5].xy() = -4.0;
            }

            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }
    };


}  // namespace


#endif