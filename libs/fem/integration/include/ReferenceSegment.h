#ifndef CFEM_REFERENCE_SEGMENT
#define CFEM_REFERENCE_SEGMENT

#include "ReferenceElement.h"

namespace cfem::reference_elem::utils{

    inline Vector3D projectFacetRefCoordToSegmentRefCoord(CfemInt localFacetId, 
                                                       const Vector3D& xiFacet) 
    {
        // Une facette d'une ligne est un point. xiFacet est ignoré.
        // Selon ReferenceSegment.h (P2), le segment de référence est [-1, 1]
        // Face 0 = Noeud 0 (-1.0), Face 1 = Noeud 1 (1.0)
        return Vector3D( (localFacetId == 0) ? -1.0 : 1.0, 0.0, 0.0 );
    }
}


namespace cfem {

    class P1ReferenceSegment : public ReferenceElement {
    public:
        P1ReferenceSegment() {
            m_nodes.resize(2);
            // Domaine de référence [-1, 1]
            m_nodes[0] = Vector3D{-1.0, 0.0, 0.0};
            m_nodes[1] = Vector3D{ 1.0, 0.0, 0.0};
        }

        int dimension() const override { return 1; }
        int getNumNodes() const override { return 2; }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 2;
            if (out.values.size() != n) out.setup(n);

            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            const CfemReal x = xi.x;

            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = 0.5 * (1.0 - x);
                out.values[1] = 0.5 * (1.0 + x);
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                out.gradients[0] = Vector3D{-0.5, 0.0, 0.0};
                out.gradients[1] = Vector3D{ 0.5, 0.0, 0.0};
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                for (auto& H : out.hessians) H.setZero();
            }
            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToSegmentRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }
    };

    class P2ReferenceSegment : public ReferenceElement {
    public:
        P2ReferenceSegment() {
            m_nodes.resize(3);
            // Domaine de référence [-1, 1], ordre VTK : extrémités puis milieu
            m_nodes[0] = Vector3D{-1.0, 0.0, 0.0};
            m_nodes[1] = Vector3D{ 1.0, 0.0, 0.0};
            m_nodes[2] = Vector3D{ 0.0, 0.0, 0.0}; // Nœud central en x=0
        }

        int dimension() const override { return 1; }
        int getNumNodes() const override { return 3; }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 3;
            if (out.values.size() != n) out.setup(n);
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            const CfemReal x = xi.x;

            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = 0.5 * x * (x - 1.0);
                out.values[1] = 0.5 * x * (x + 1.0);
                out.values[2] = 1.0 - x * x;
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                out.gradients[0] = Vector3D{x - 0.5, 0.0, 0.0};
                out.gradients[1] = Vector3D{x + 0.5, 0.0, 0.0};
                out.gradients[2] = Vector3D{-2.0 * x, 0.0, 0.0};
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                for (auto& H : out.hessians) H.setZero();
                out.hessians[0].xx() =  1.0;
                out.hessians[1].xx() =  1.0;
                out.hessians[2].xx() = -2.0;
            }
            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToSegmentRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }
    };
}

#endif