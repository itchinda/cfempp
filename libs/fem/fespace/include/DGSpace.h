//************************************************************************
// --- CFEM++ Library
// ---
// --- Copyright 2024-2026 Ismaël Tchinda Ngueyong et al.
// --- ALL RIGHTS RESERVED
// ---
// --- This software is protected by international copyright laws.
// --- Use, distribution, or modification of this software in any form,
// --- source or binary, for personal, academic, or non-commercial 
// --- purposes is permitted free of charge, provided that this 
// --- copyright notice and this permission notice appear in all 
// --- copies and supporting documentation.
// ---
// --- The authors and contributors provide this code "as is" without 
// --- any express or implied warranty. In no event shall they be 
// --- held liable for any damages arising from the use of this software.
//************************************************************************

#ifndef CFEM_DGSPACE_H
#define CFEM_DGSPACE_H

#include "FESpace.h"
#include <map>

namespace cfem {

/**
 * @class DGSpace
 * @brief Discontinuous Galerkin (DG) Space. 
 * DOFs are local to each cell; no nodes are shared between neighbors.
 */
class DGSpace : public FESpace {
private:
    struct ConstructorToken {};
    
    CfemInt m_totalDofs = 0;

    std::vector<CfemInt> m_cellOffsets;
    void enumerateDofs() override;

public:
    DGSpace(ConstructorToken, std::shared_ptr<Mesh> mesh, Order order, int dimension = 1);
    
    static std::shared_ptr<DGSpace> create(std::shared_ptr<Mesh> mesh, Order order, FieldType type = FieldType::Scalar);

    // --- FESpace Overrides ---
    CfemInt getGlobalNumDofs() const override { return m_totalDofs; }
    int getDofsPerCell(CfemInt cellId) const override;
    void getElementDofs(CfemInt cellId, std::vector<CfemInt>& globalIndices) const override;
    
    // --- DG Specific Abstraction ---
    CfemInt getNumUniquePositions() const override { return m_totalDofs / m_dimension; }
    
    // In DG, the position ID is unique per (Cell, LocalNode)
    CfemInt getPositionId(CfemInt cellId, int localNodeIdx) const override;
    
    // Global index mapping
    CfemInt getGlobalDofIdx(CfemInt cellId, int localNodeIdx, int component) const override;
    
    // Legacy/Vertex mapping: In DG, a mesh node does not map to a single DOF!
    // We return -1 or the first occurrence depending on policy.
    CfemInt getGlobalDofFromMeshNode(CfemInt meshNodeId, int component) const override { 
        CFEM_ERROR <<  "DGSpace::getGlobalDofFromMeshNode() called. In Discontinuous Galerkin spaces, "
                   << "nodes are not uniquely mapped to global DOFs because they are duplicated across "
                   << "element interfaces. Use getGlobalDofIdx(cellId, localNodeId, component) instead.";
        return -1;
    }

    void computeSparsityPattern(std::vector<CfemInt>& nnz) const override;
    void getEntityDofs(int internalEntityId, std::vector<CfemInt>& boundaryDofs) const override;

    void evaluateBasis(CfemInt cellId, const Vector3D& xi, std::vector<CfemReal>& values) const override;
    void evaluateGradients(CfemInt cellId, const Vector3D& xi, std::vector<Vector3D>& grads) const override;
};

} // namespace cfem
#endif