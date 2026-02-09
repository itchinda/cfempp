#ifndef CFEM_SCALAR_FIELD_H
#define CFEM_SCALAR_FIELD_H

#include "TypedField.h"


namespace cfem {
class ScalarField : public TypedField<CfemReal> {
public:
    using TypedField<CfemReal>::TypedField;

protected:
    CfemReal getComponentValue(const std::vector<CfemInt>& dofs, int localNodeIdx) const override {
        return m_data[dofs[localNodeIdx]];
    }
};

}

#endif
