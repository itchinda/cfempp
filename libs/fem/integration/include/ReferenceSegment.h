#ifndef CFEM_REFERENCE_SEGMENT
#define CFEM_REFERENCE_SEGMENT

#include "ReferenceElement.h"

namespace cfem {

    class P1ReferenceSegment : public ReferenceElement {
    public:
        P1ReferenceSegment() {
            m_nodes.resize(2);
            m_nodes[0] = Vector3D{0.0, 0.0, 0.0};
            m_nodes[1] = Vector3D{1.0, 0.0, 0.0};
        }

        int dimension() const override { return 1; }
        int numNodes() const override { return 2; }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 2;
            if (out.values.size() != n) out.setup(n);

            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            const CfemReal x = xi.x;

            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = 1.0 - x;
                out.values[1] = x;
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                out.gradients[0] = Vector3D{-1.0, 0.0, 0.0};
                out.gradients[1] = Vector3D{ 1.0, 0.0, 0.0};
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                for (auto& H : out.hessians) H.setZero();
            }
            out.validFlags = requested;
        }
    };

    class P2ReferenceSegment : public ReferenceElement {
    public:
        P2ReferenceSegment() {
            m_nodes.resize(3);
            m_nodes[0] = Vector3D{0.0, 0.0, 0.0};
            m_nodes[1] = Vector3D{1.0, 0.0, 0.0};
            m_nodes[2] = Vector3D{0.5, 0.0, 0.0}; // Nœud central
        }

        int dimension() const override { return 1; }
        int numNodes() const override { return 3; }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 3;
            if (out.values.size() != n) out.setup(n);
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            const CfemReal x = xi.x;

            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = (1.0 - x) * (1.0 - 2.0 * x);
                out.values[1] = x * (2.0 * x - 1.0);
                out.values[2] = 4.0 * x * (1.0 - x);
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                out.gradients[0] = Vector3D{4.0 * x - 3.0, 0.0, 0.0};
                out.gradients[1] = Vector3D{4.0 * x - 1.0, 0.0, 0.0};
                out.gradients[2] = Vector3D{4.0 - 8.0 * x, 0.0, 0.0};
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                for (auto& H : out.hessians) H.setZero();
                out.hessians[0].xx() = 4.0;
                out.hessians[1].xx() = 4.0;
                out.hessians[2].xx() = -8.0;
            }
            out.validFlags = requested;
        }
    };
}

#endif