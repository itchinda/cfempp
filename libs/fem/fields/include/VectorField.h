#ifndef VECTOR_FIELD_H
#define VECTOR_FIELD_H

#include "TypedField.h"

namespace cfem {

class VectorField : public TypedField<Vector3D> {
public:
    using TypedField<Vector3D>::TypedField;

protected:
    Vector3D getComponentValue(const std::vector<CfemInt>& dofs, int localNodeIdx) const override {
        CfemInt start = dofs[localNodeIdx * 3]; // simplified indexing logic
        return {m_data[start], m_data[start+1], m_data[start+2]};
    }
};

}
#endif