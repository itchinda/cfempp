#ifndef CFEM_LAGRANGE_SPACE_H
#define CFEM_LAGRANGE_SPACE_H

#include "FESpace.h"

class LagrangeSpace : public FESpace {
public:
    using FESpace::FESpace; // Use the base constructor

    CfemInt getGlobalNumDofs() const override {
        if (m_continuity == Continuity::CG) {
            // In CG, DoFs are attached to vertices/nodes. 
            // If P2, this includes mid-edge/mid-face nodes already in the Mesh.
            return m_mesh.getNumVertices() * m_dimension;
        } else {
            // In DG, every cell has its own private set of DoFs.
            CfemInt total = 0;
            for (CfemInt i = 0; i < m_mesh.getNumCells(); ++i) {
                total += getNodesPerCell(m_mesh.getCell(i).type) * m_dimension;
            }
            return total;
        }
    }

    int getNodesPerCell(CellType type) const {
        // If the user requested a P1 Space, we only use the corner nodes
        // regardless of whether the underlying mesh cell has 4, 8, or 9 nodes.
        if (m_order == Order::P1) {
            if (type == CellType::LINE2 || type == CellType::LINE3) return 2;
            if (type == CellType::TRI3 || type == CellType::TRI6) return 3;
            if (type == CellType::QUAD4 || type == CellType::QUAD9) return 4;
            if (type == CellType::TET4 || type == CellType::TET10) return 4;
            if (type == CellType::PRISM6 || type == CellType::PRISM18) return 6;
            if (type == CellType::HEX8 || type == CellType::HEX27) return 8;
        }
        
        // If P2, we use the full set
        if (m_order == Order::P2) {
            if ( type == CellType::LINE3) return 3;
            if ( type == CellType::TRI6) return 6;
            if ( type == CellType::QUAD9) return 9;
            if ( type == CellType::TET10) return 10;
            if ( type == CellType::PRISM18) return 18;
            if ( type == CellType::HEX27) return 27;
        }

        CFEM_ERROR("UNKNOWN CELL TYPE");
        return 0;
    }

    void getElementDofs(CfemInt cellId, std::vector<CfemInt>& globalIndices) const override {
        const auto& cell = m_mesh.getCell(cellId);
        int nNodes = getDofsPerCell(cell.type) / m_dimension; // Local nodes
        globalIndices.clear();

        if (m_continuity == Continuity::CG) {
            for (int i = 0; i < nNodes; ++i) {
                CfemInt meshNodeId = cell.vertex_ids[i];
                for (int d = 0; d < m_dimension; ++d) {
                    globalIndices.push_back(meshNodeId * m_dimension + d);
                }
            }
        } else {
            // DG: Offset by previous cells' total DoFs
            // For a uniform mesh, this is cellId * nodesPerCell * m_dimension
            CfemInt offset = 0; 
            for(CfemInt i=0; i < cellId; ++i) offset += getDofsPerCell(m_mesh.getCell(i).type);
            
            for (int i = 0; i < getDofsPerCell(cell.type); ++i) {
                globalIndices.push_back(offset + i);
            }
        }
    }

    void evaluateBasis(CellType type, const Vector3D& xi, std::vector<CfemReal>& values) const override {
        // This will call your Shape Function library
        // e.g., ShapeFunctions::evaluate(type, m_order, xi, values);
    }
};

#endif