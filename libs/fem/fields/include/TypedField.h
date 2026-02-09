#ifndef CFEM_TYPEDFIELD_H
#define CFEM_TYPEDFIELD_H

#include "Field.h"


namespace cfem {
    
/**
 * @class TypedField
 * @brief Intermediate template to handle data storage and interpolation logic.
 * @tparam T The physical type (CfemReal, Vector3D, etc.)
 */
template <typename T>
class TypedField : public Field {
public:
    TypedField(std::shared_ptr<FESpace> space, const std::string& name)
        : Field(space, name) {
        m_data.resize(m_space->getGlobalNumDofs(), 0.0);
    }

    std::vector<CfemReal>& data() override { return m_data; }
    const std::vector<CfemReal>& data() const override { return m_data; }

    /** * @brief The core FEM interpolation: u(xi) = sum( Ni(xi) * ui )
     */
    virtual T evaluateAt(CfemInt cellId, const Vector3D& xi) const {
        std::vector<CfemInt> dofs;
        m_space->getElementDofs(cellId, dofs);

        std::vector<CfemReal> basis;
        m_space->evaluateBasis(m_mesh.getCell(cellId).type, xi, basis);

        T result = T{}; // Assumes T has a default constructor and supports += / * operators
        int nNodes = basis.size();
        int dim = m_space->getDimension();

        for (int i = 0; i < nNodes; ++i) {
            // This logic maps the flat m_data to the type T
            result += getComponentValue(dofs, i) * basis[i];
        }
        return result;
    }

protected:
    std::vector<CfemReal> m_data;

    // Internal helper to extract the T-sized chunk of data for node i
    virtual T getComponentValue(const std::vector<CfemInt>& dofs, int localNodeIdx) const = 0;
};

}
#endif