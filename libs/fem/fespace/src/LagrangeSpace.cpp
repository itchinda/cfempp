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

#include "LagrangeSpace.h"
#include "ReferenceElements.h"
#include "Mesh.h"

#include <set>
#include <algorithm>

namespace cfem {

    /**
     * @brief Constructor for the Continuous Lagrange Finite Element Space.
     * @param mesh Shared pointer to the underlying mesh.
     * @param order Polynomial order (P1 or P2).
     * @param field_type The field type (Scalar, Vector2d/3d Tensor2d/3d etc.).
     */
    LagrangeSpace::LagrangeSpace(ConstructorToken, std::shared_ptr<Mesh> mesh, Order order, FieldType field_type)
        : FESpace(mesh, order, Continuity::CG, field_type) 
    {
        if (m_order == Order::P0){
            CFEM_ERROR << "Invalid Order: LagrangeSpace is designed for C0-continuity and requires at least Order::P1. "
                    << "For piecewise constants, use a DGSpace.";
        }
    }

    void LagrangeSpace::validate() const {
        FESpace::validate();
        MeshOrder meshOrder = m_mesh->getOrder();
        if (m_order == Order::P2 && (meshOrder == MeshOrder::None || meshOrder == MeshOrder::Linear)) {
            CFEM_ERROR << "LagrangeSpace Order Mismatch: P2 Space requires a Quadratic Mesh.";
        }
    }

    /**
     * @brief Factory method: The recommended way to instantiate a LagrangeSpace.
     * @param mesh Shared pointer to the geometry.
     * @param order Polynomial order (P1 or P2). P0 is not supported here.
     * @param field_type Scalar or Vectorial nature of the physics.
     * @return A shared pointer to the initialized and enumerated LagrangeSpace.
     * @throw CFEM_ERROR if order is P0 or if mesh cannot support requested order.
     */
    std::shared_ptr<LagrangeSpace> LagrangeSpace::create(std::shared_ptr<Mesh> mesh, 
                                                 Order order, 
                                                 FieldType field_type) {

        if (order == Order::P0) {
            CFEM_ERROR << "LagrangeSpace Error: P0 is a discontinuous space and is not supported by LagrangeSpace. "
                       "Please use DGSpace or upgrade the order to P1.";
        }

        if (order == Order::P2 && mesh->getOrder() == MeshOrder::Linear) {
            CFEM_INFO << "LagrangeSpace: Auto-upgrading mesh to quadratic to support P2 space.";
            mesh->upgradeToQuadratic(); 
        }

        // Instantiation & Initialization
        auto space = std::make_shared<LagrangeSpace>(ConstructorToken{}, mesh, order, field_type);
        
        // Connect mesh updates to this space
        mesh->addObserver(space);
        
        space->validate();
        space->enumerateDofs();

        return space;
    }


/**
 * @brief Determines the number of active nodes per cell based on the space order.
 * This filters higher-order nodes if the space is only P1.
 */
CfemInt LagrangeSpace::getNodesPerCell(CellType type) const {
    if (m_order == Order::P1) {
        switch (type) {
            case CellType::LINE2:  case CellType::LINE3:   return 2;
            case CellType::TRI3:   case CellType::TRI6:    return 3;
            case CellType::QUAD4:  case CellType::QUAD9:   return 4;
            case CellType::TET4:   case CellType::TET10:   return 4;
            case CellType::PRISM6: case CellType::PRISM18: return 6;
            case CellType::HEX8:   case CellType::HEX27:   return 8;
            default: return 0;
        }
    }
    if (m_order == Order::P2) {
        switch (type) {
            case CellType::LINE3:    return 3;
            case CellType::TRI6:     return 6;
            case CellType::QUAD9:    return 9;
            case CellType::TET10:    return 10;
            case CellType::PRISM18:  return 18;
            case CellType::HEX27:    return 27;
            default: return 0;
        }
    }
    if (m_order == Order::Mini) {
        switch (type) {
            case CellType::LINE2:  case CellType::LINE3:   return 3;
            case CellType::TRI3:   case CellType::TRI6:    return 4;
            case CellType::QUAD4:  case CellType::QUAD9:   return 5;
            case CellType::TET4:   case CellType::TET10:   return 5;
            case CellType::PRISM6: case CellType::PRISM18: return 7;
            case CellType::HEX8:   case CellType::HEX27:   return 9;
            default: return 0;
        }
    }
    return 0;
}

/**
 * @brief Assigns local and global IDs to active mesh vertices.
 * Filters nodes based on polynomial order and builds bi-directional mapping tables.
 */
void LagrangeSpace::enumerateDofs() {
    const CfemInt numMeshNodes = static_cast<CfemInt>(m_mesh->getNumVertices());
    const CfemInt numCells = static_cast<CfemInt>(m_mesh->getNumCells());
    
    // Initialisation des tables de correspondance
    m_nodeToLocalIdx.assign(numMeshNodes, -1);
    
    CfemInt counter = 0;

    // PREMIER PASSAGE : Énumération des nœuds du maillage (Partagés)
    // On ne boucle que sur les sommets qui appartiennent physiquement au maillage
    // Pour l'élément Mini, ce sont les 8 premiers nœuds (HEX8) ou 4 (TET4).
    for (CfemInt c = 0; c < numCells; ++c) {
        const auto& cell = m_mesh->getCell(c);
        const auto& cell_ids = m_mesh->getCellConnectivity(c);
        
        // On détermine combien de nœuds de la connectivité du maillage on doit utiliser
        // Pour Order::Mini sur HEX8, nodesInMesh = 8, même si getNodesPerCell() renvoie 9.
        CfemInt nodesInMesh = 0;
        if (m_order == Order::P1) nodesInMesh = getNodesPerCell(cell.type);
        else if (m_order == Order::Mini) {
            // Le Mini utilise les sommets linéaires du maillage + 1 bulle interne
            switch (cell.type) {
                case CellType::TRI3:  case CellType::TRI6:  nodesInMesh = 3; break;
                case CellType::QUAD4: case CellType::QUAD9: nodesInMesh = 4; break;
                case CellType::TET4:  case CellType::TET10: nodesInMesh = 4; break;
                case CellType::HEX8:  case CellType::HEX27: nodesInMesh = 8; break;
                default: nodesInMesh = 0;
            }
        }
        else if (m_order == Order::P2) nodesInMesh = getNodesPerCell(cell.type);

        for (CfemInt n = 0; n < nodesInMesh; ++n) {
            CfemInt vId = cell_ids[n];
            if (m_nodeToLocalIdx[vId] == -1) {
                m_nodeToLocalIdx[vId] = counter++;
            }
        }
    }

    // SECOND PASSAGE : Énumération des Bulles (DOFs internes, non partagés)
    // On stocke l'ID de départ des bulles pour faciliter le mapping dans getGlobalDofIdx
    m_firstBubbleLocalId = counter; 
    
    if (m_order == Order::Mini) {
        m_cellToBubbleIdx.assign(numCells, -1);
        for (CfemInt c = 0; c < numCells; ++c) {
            // Chaque cellule reçoit un index unique pour sa bulle
            m_cellToBubbleIdx[c] = counter++;
        }
    }

    m_numActiveNodes = counter; // Sommets uniques + Bulles uniques
    m_totalDofs = m_numActiveNodes * m_numComponents;

    // 4. Construction de la table inverse (Local -> Mesh)
    // Note : Les bulles n'ont pas de correspondance directe avec un sommet du maillage
    m_localToMeshNodeIdx.assign(m_numActiveNodes, -1);
    for (CfemInt meshId = 0; meshId < (CfemInt)m_nodeToLocalIdx.size(); ++meshId) {
        CfemInt local = m_nodeToLocalIdx[meshId];
        if (local != -1) {
            m_localToMeshNodeIdx[local] = meshId;
        }
    }
}

CfemInt LagrangeSpace::getPositionId(CfemInt cellId, CfemInt localNodeIdx) const {
    const auto& cell = m_mesh->getCell(cellId);
    const auto& connectivity = m_mesh->getCellConnectivity(cellId);

    // Si c'est un nœud de sommet (0-7 pour Hexa)
    if (localNodeIdx < (int)connectivity.size()) {
        return m_nodeToLocalIdx[connectivity[localNodeIdx]];
    }
    
    // Si c'est le nœud bulle (index 8 pour Hexa-Mini)
    if (m_order == Order::Mini && localNodeIdx == getNodesPerCell(cell.type) - 1) {
        return m_cellToBubbleIdx[cellId];
    }

    return -1;
}

CfemInt LagrangeSpace::getGlobalDofIdx(CfemInt cellId, CfemInt localNodeIdx, CfemInt component) const {
    CfemInt posId = getPositionId(cellId, localNodeIdx);
    if (posId == -1) return -1;
    // BLOCK-WISE ORDERING: [Comp 0 Nodes][Comp 1 Nodes]...[Comp N Nodes]
    return posId + (component * m_numActiveNodes);
}

CfemInt LagrangeSpace::getGlobalDofFromMeshNode(CfemInt meshNodeId, CfemInt component) const {
    CfemInt posId = m_nodeToLocalIdx[meshNodeId];
    if (posId == -1) return -1;
    return posId + (component * m_numActiveNodes);
}

// --- ELEMENT INTERFACE ---

CfemInt LagrangeSpace::getDofsPerCell(CfemInt cellId) const {
    return getNodesPerCell(m_mesh->getCell(cellId).type) * m_numComponents;
}

void LagrangeSpace::getElementDofs(CfemInt cellId, std::vector<CfemInt>& globalIndices) const {
    const auto& cell = m_mesh->getCell(cellId);
    CfemInt nodesInCell = getNodesPerCell(cell.type);
    const auto& cell_ids = m_mesh->getCellConnectivity(cellId);
    
    globalIndices.resize(nodesInCell * m_numComponents);

    for (CfemInt n = 0; n < nodesInCell; ++n) {
        CfemInt nodeIdx = m_nodeToLocalIdx[cell_ids[n]];
        for (CfemInt d = 0; d < m_numComponents; ++d) {
            // Fill with block-wise indices: node first, then component blocks
            globalIndices[n + d * nodesInCell] = nodeIdx + (d * m_numActiveNodes);
        }
    }
}

// --- MATH & EVALUATION ---

void LagrangeSpace::evaluateBasis(CfemInt cellId, const Vector3D& xi, std::vector<CfemReal>& values) const {
    const Cell& cell = m_mesh->getCell(cellId);
    const ReferenceElement& refEl = getRefElement(cell.type);

    ShapeInfo shape;
    refEl.evaluate(xi, ShapeEvaluationFlags::Value, shape);

    values.assign(shape.values.begin(), shape.values.end());
}

void LagrangeSpace::evaluateGradients(CfemInt cellId, const Vector3D& xi, std::vector<Vector3D>& grads) const {
    const Cell& cell = m_mesh->getCell(cellId);
    const ReferenceElement& refEl = getRefElement(cell.type);
    
    // Thread-local scratchpad to ensure OpenMP safety and avoid allocations
    static thread_local DynamicVector<Vector3D> t_physicalCoords;
    m_mesh->getCellVerticesCoords(cellId, t_physicalCoords);

    ShapeInfo shape;
    MappingInfo map;
    refEl.computeMapping(xi, t_physicalCoords, ShapeEvaluationFlags::Gradient, shape, map);

    grads.assign(map.physicalShapeGradients.begin(), map.physicalShapeGradients.end());
}

// --- PETSc & SYSTEM SUPPORT ---

/**
 * @brief Computes the number of non-zeros (nnz) per row for PETSc matrix pre-allocation.
 * Respects the block-wise numbering by replicating connectivity counts across components.
 */
void LagrangeSpace::computeSparsityPattern(std::vector<CfemInt>& nnz) const {
    nnz.assign(m_totalDofs, 0);

    for (CfemInt i = 0; i < (CfemInt)m_mesh->getNumVertices(); ++i) {
        CfemInt nodeIdx = m_nodeToLocalIdx[i];
        if (nodeIdx == -1) continue;

        const auto& vertex = m_mesh->getVertex(i);
        std::set<CfemInt> neighborNodes;

        // Traverse connected cells to find all adjacent active nodes
        for (CfemInt cellId : vertex.getConnectedCells()){
            const auto& cell = m_mesh->getCell(cellId);
            const auto& cell_ids = m_mesh->getCellConnectivity(cellId);

            CfemInt nInCell = getNodesPerCell(cell.type);
            for (CfemInt n = 0; n < nInCell; ++n) {
                neighborNodes.insert(m_nodeToLocalIdx[cell_ids[n]]);
            }
        }

        // Each component in the block has the same number of neighbors * dimensions
        CfemInt rowNNZ = static_cast<CfemInt>(neighborNodes.size()) * m_numComponents;

        for (CfemInt d = 0; d < m_numComponents; ++d) {
            nnz[nodeIdx + (d * m_numActiveNodes)] = rowNNZ;
        }
    }
}

void LagrangeSpace::getEntityDofs(CfemInt internalEntityId, std::vector<CfemInt>& boundaryDofs) const {
    boundaryDofs.clear();
    const auto& boundaryNodes = m_mesh->getEntityConst(internalEntityId).getVertexIds();

    for (CfemInt vId : boundaryNodes) {
        CfemInt nodeIdx = m_nodeToLocalIdx[vId];
        if (nodeIdx != -1) {
            for (CfemInt d = 0; d < m_numComponents; ++d) {
                boundaryDofs.push_back(nodeIdx + (d * m_numActiveNodes));
            }
        }
    }
}

} // namespace cfem