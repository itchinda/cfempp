#ifndef CFEM_FESPACE_H
#define CFEM_FESPACE_H

#include "Mesh.h"

namespace cfem {

class FESpace {
public:
    enum class Order { P0 = 0, P1 = 1, P2 = 2 };
    enum class Continuity { CG, DG };

    FESpace(const Mesh& mesh, Order order, Continuity cont, int dimension = 1)
        : m_mesh(mesh), m_order(order), m_continuity(cont), m_dimension(dimension) {}

    virtual ~FESpace() = default;

    CfemInt getDimension() const { return m_dimension; } 
    
    virtual CfemInt getGlobalNumDofs() const = 0;
    virtual int getDofsPerCell(CellType type) const = 0;
    virtual void getElementDofs(CfemInt cellId, std::vector<CfemInt>& globalIndices) const = 0;
    virtual void evaluateBasis(CellType type, const cfem::Vector3D& xi, std::vector<CfemReal>& values) const = 0;

protected:
    const Mesh& m_mesh;
    Order m_order;
    Continuity m_continuity;
    int m_dimension; 
};
}

#endif