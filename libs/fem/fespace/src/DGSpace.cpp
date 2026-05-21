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

#include "DGSpace.h"
#include "Mesh.h"
#include "ReferenceElementFactory.h"

namespace cfem {

DGSpace::DGSpace(ConstructorToken, std::shared_ptr<Mesh> mesh, Order order, FieldType field_type)
    : FESpace(mesh, order, Continuity::DG, field_type) 
{
    if (order == Order::Mini){
        CFEM_THROW("Mini elements are not supported in DGSpace. "
                   << "Mini elements (P1+bubble) are designed for H1-stability; for DG, use P1 or P2 instead.");
    }

}

std::shared_ptr<DGSpace> DGSpace::create(std::shared_ptr<Mesh> mesh, 
                                         Order order, 
                                         FieldType field_type) {

    if (order == Order::Mini){
        CFEM_THROW("Mini elements are not supported in DGSpace. "
                   << "Mini elements (P1+bubble) are designed for H1-stability; for DG, use P1 or P2 instead.");
    }

    if (order == Order::P2 && mesh->getOrder() == MeshOrder::Linear) {
        CFEM_INFO << "DGSpace: Auto-upgrading mesh to quadratic to support P2 space.";
        mesh->upgradeToQuadratic(); 
    }

    // Instantiation & Initialization
    auto space = std::make_shared<DGSpace>(ConstructorToken{}, mesh, order, field_type);
    
    // Connect mesh updates to this space
    mesh->addObserver(space);
    
    space->validate();
    space->enumerateDofs();
    return space;
}


/**
 * @brief In DG, numbering is trivial. 
 * DOFs are simply stacked cell by cell.
 */
void DGSpace::enumerateDofs() {
    const CfemInt numCells = m_mesh->getNumCells();
    m_cellOffsets.assign(numCells + 1, 0);
    
    CfemInt currentOffset = 0;
    for (CfemInt i = 0; i < numCells; ++i) {
        m_cellOffsets[i] = currentOffset;
        const Cell& cell = m_mesh->getCell(i);
        const ReferenceElement& refEl = ReferenceElementFactory::getReferenceElement( cell.type, getOrderInt() );
        currentOffset += refEl.getNumNodes();
    }
    m_cellOffsets[numCells] = currentOffset; // Offset total final
    
    m_totalDofs = currentOffset * m_numComponents;
}

CfemInt DGSpace::getPositionId(CfemInt cellId, CfemInt localNodeIdx) const {
    return m_cellOffsets[cellId] + localNodeIdx;
}

CfemInt DGSpace::getGlobalDofIdx(CfemInt cellId, CfemInt localNodeIdx, CfemInt component) const {
    CfemInt posId = getPositionId(cellId, localNodeIdx);
    CfemInt numUniqueNodes = m_totalDofs / m_numComponents;
    
    // Keep the same BLOCK-wise logic: [All Comp 0][All Comp 1]...
    return posId + (component * numUniqueNodes);
}

CfemInt DGSpace::getDofsPerCell(CfemInt cellId) const {
    const Cell& cell = m_mesh->getCell(cellId);
    const ReferenceElement& refEl = ReferenceElementFactory::getReferenceElement( cell.type, getOrderInt() );
    return refEl.getNumNodes() * m_numComponents;
}

void DGSpace::getElementDofs(CfemInt cellId, std::vector<CfemInt>& globalIndices) const {
    const Cell& cell = m_mesh->getCell(cellId);
    const ReferenceElement& refEl = ReferenceElementFactory::getReferenceElement( cell.type, getOrderInt() );
    CfemInt nNodes = refEl.getNumNodes();
    globalIndices.resize(nNodes * m_numComponents);
    
    for (CfemInt c = 0; c < m_numComponents; ++c) {
        for (CfemInt i = 0; i < nNodes; ++i) {
            globalIndices[i + c * nNodes] = getGlobalDofIdx(cellId, i, c);
        }
    }
}

/**
 * @brief Calcule le pattern de sparsité exact pour le DG en utilisant l'adjacence du maillage.
 * * En DG, une ligne i de la matrice (associée à un DOF de la cellule C) possède des entrées non-nulles 
 * uniquement pour :
 * 1. Tous les DOFs de la cellule C elle-même (bloc diagonal).
 * 2. Tous les DOFs des cellules voisines partageant une face avec C (blocs extra-diagonaux).
 */
void DGSpace::computeSparsityPattern(std::vector<CfemInt>& nnz) const {
    const CfemInt numCells = m_mesh->getNumCells();
    nnz.assign(m_totalDofs, 0);

    for (CfemInt i = 0; i < numCells; ++i) {
        const Cell& cell = m_mesh->getCell(i);
        
        // 1. DOFs de la cellule elle-même
        CfemInt dofsSelf = getDofsPerCell(i);
        
        // 2. DOFs des voisins (en parcourant l'adjacence construite par buildAdjacency)
        CfemInt dofsNeighbors = 0;
        
        // On utilise directement les pointeurs initialisés par Mesh::buildAdjacency()
        // cell.num_faces correspond au nombre d'entrées dans le bloc CSR pour cette cellule
        for (CfemInt f = 0; f < cell.num_faces; ++f) {
            CfemInt neighborId = cell.neighbor_cells_ptr[f];
            
            // Si neighborId == -1, c'est une condition aux limites (pas de cellule voisine)
            if (neighborId != -1) {
                dofsNeighbors += getDofsPerCell(neighborId);
            }
        }

        // Le nombre total de non-zéros pour chaque ligne (DOF) de cette cellule
        CfemInt totalRowNNZ = dofsSelf + dofsNeighbors;

        // On remplit le tableau nnz pour chaque composante et chaque noeud local de la cellule
        // Rappel : Notre numérotation est Block-wise [Comp 0][Comp 1]...
        for (CfemInt comp = 0; comp < m_numComponents; ++comp) {
            const ReferenceElement& refEl = ReferenceElementFactory::getReferenceElement( cell.type, getOrderInt() );
            CfemInt nNodes = refEl.getNumNodes();
            for (CfemInt n = 0; n < nNodes; ++n) {
                CfemInt globalIdx = getGlobalDofIdx(i, n, comp);
                nnz[globalIdx] = totalRowNNZ;
            }
        }
    }
}


    /**
     * @brief Evaluates the reference basis function values at a local coordinate xi within a specific cell.
     * @param cellId Local index of the computational cell.
     * @param xi Reference/Local coordinates inside the cell's geometric space.
     * @param values Destination vector container to hold the evaluated scalar basis functions.
     */
    void DGSpace::evaluateBasis(CfemInt cellId, const Vector3D& xi, std::vector<CfemReal>& values) const 
    {
        ShapeInfo shape;
        const Cell& cell = m_mesh->getCell(cellId);
        const ReferenceElement& refEl = ReferenceElementFactory::getReferenceElement(cell.type, getOrderInt());
        
        // Evaluate reference basis values directly into the ShapeInfo container
        shape.values.resize(refEl.getNumNodes());
        refEl.evaluateShapesValues(xi, shape.values.data());
        
        // Direct contiguous copy remains valid for scalar values
        values.assign(shape.values.begin(), shape.values.end());
    }

    /**
     * @brief Evaluates the physical gradients of the basis functions at a local coordinate xi.
     * @param cellId Local index of the computational cell.
     * @param xi Reference/Local coordinates inside the cell's geometric space.
     * @param grads Destination vector container to pack the calculated physical Vector3D gradients.
     */
    void DGSpace::evaluateBasisGradients(CfemInt cellId, const Vector3D& xi, std::vector<Vector3D>& grads) const 
    {
        // Thread-local scratchpad to guarantee OpenMP safety and bypass dynamic allocations
        static thread_local DynamicVector<Vector3D> t_physCoords;
        m_mesh->getCellVerticesCoords(cellId, t_physCoords);

        const Cell& cell = m_mesh->getCell(cellId);
        const ReferenceElement& refEl = ReferenceElementFactory::getReferenceElement(cell.type, getOrderInt());
        
        ShapeInfo shape;
        MappingInfo map;
        
        refEl.computeMapping(
            xi, 
            t_physCoords, 
            ShapeEvaluationFlags::FirstDerivatives, 
            MappingEvaluationFlags::PhysicalFirstDerivatives,
            shape, 
            map
        );

        const CfemInt nNodes = static_cast<CfemInt>(map.dN_dx.size());
        const CfemInt spatialDimension = m_mesh->getSpatialDimension();

        // Resize the destination vector to accommodate the layout without memory spikes
        grads.resize(nNodes);

        for (CfemInt i = 0; i < nNodes; ++i)
        {
            grads[i].x = map.dN_dx[i];
            grads[i].y = (spatialDimension >= 2) ? map.dN_dy[i] : 0.0;
            grads[i].z = (spatialDimension == 3) ? map.dN_dz[i] : 0.0;
        }
    }

/**
 * @brief Récupère les DOFs associés à une entité géométrique (bord ou groupe physique).
 * @param internalEntityId L'identifiant (Gmsh tag) de l'entité.
 * @param entityDofs Vecteur à remplir avec les indices globaux.
 * * En DG, pour une entité de dimension (dim-1), on identifie les facettes appartenant 
 * à cette entité, puis on récupère les DOFs de la cellule adjacente (left_cell) 
 * qui sont localisés sur cette face.
 */
void DGSpace::getEntityDofs(CfemInt internalEntityId, std::vector<CfemInt>& entityDofs) const {
    entityDofs.clear();

    // 1. Récupérer toutes les facettes appartenant à cette entité physique
    const std::vector<CfemInt>& facetIds = m_mesh->getEntityConst(internalEntityId).getFacetIds();
    
    // Pour éviter les doublons si plusieurs facettes touchent la même cellule 
    // (rare mais possible selon la topologie)
    std::set<CfemInt> dofSet;

    for (CfemInt fId : facetIds) {
        const Facet& facet = m_mesh->getFacet(fId);
        
        // En DG, on s'intéresse à la cellule qui "possède" cette facette de bord
        // Par convention, les facettes de bord ont leur cellule à gauche (left_cell)
        CfemInt cellId = facet.left_cell;
        if (cellId == -1) continue; 

        const auto& cell = m_mesh->getCell(cellId);
        const ReferenceElement& refEl = ReferenceElementFactory::getReferenceElement( cell.type, getOrderInt() );
        
        // 2. Trouver quels nœuds locaux de la cellule sont sur cette facette précise
        // On récupère la connectivité de la facette pour comparer
        auto facetNodes = m_mesh->getFacetConnectivity(fId);
        
        // On parcourt les faces locales de la cellule pour identifier celle qui correspond à fId
        for (CfemInt iFace = 0; iFace < cell.num_faces; ++iFace) {
            if (cell.facet_ids_ptr[iFace] == fId) {
                // On a trouvé la face locale correspondante !
                // On récupère les indices des nœuds locaux de la cellule appartenant à cette face
                auto localNodesOnFace = refEl.getFaceNodes(iFace);

                // 3. Ajouter les DOFs globaux pour chaque composante
                for (CfemInt comp = 0; comp < m_numComponents; ++comp) {
                    for (CfemInt localIdx : localNodesOnFace) {
                        CfemInt globalIdx = getGlobalDofIdx(cellId, localIdx, comp);
                        dofSet.insert(globalIdx);
                    }
                }
                break; // Facette traitée pour cette cellule
            }
        }
    }

    // Transfert du set vers le vecteur de sortie
    entityDofs.assign(dofSet.begin(), dofSet.end());
}

} // namespace cfem