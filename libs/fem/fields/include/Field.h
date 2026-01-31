#ifndef FIELD_H
#define FIELD_H

#include "FESpace.h"
/**
 * @class Field
 * @brief Abstract base for all physical fields (Temperature, Velocity, etc).
 */
class Field {
public:
    Field(std::shared_ptr<FESpace> space, const std::string& name)
        : m_space(space), m_name(name) {}

    virtual ~Field() = default;

    std::string getName() const { return m_name; }
    std::shared_ptr<FESpace> getSpace() const { return m_space; }

    /** @brief Returns the raw data vector for linear algebra solvers. */
    virtual std::vector<CfemReal>& data() = 0;
    virtual const std::vector<CfemReal>& data() const = 0;

    /** @brief Interpolates the field value at a physical point using a cell's local coords. */
    // This will be implemented in the derived typed fields.
    
protected:
    std::shared_ptr<FESpace> m_space;
    std::string m_name;
};

#endif