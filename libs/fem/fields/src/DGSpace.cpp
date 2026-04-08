#include "DGSpace.h"
#include "ReferenceElements.h" // Assuming a factory or manual init

namespace cfem {

DGSpace::DGSpace(ConstructorToken, std::shared_ptr<Mesh> mesh, Order order, int dimension)
    : FESpace(mesh, order, Continuity::DG, dimension) 
{
    initReferenceElements();
}

// std::shared_ptr<DGSpace> DGSpace::create(std::shared_ptr<Mesh> mesh, Order order, FieldType type) {
//     auto space = std::make_shared<DGSpace>(ConstructorToken{}, mesh, order, static_cast<int>(type));
//     space->enumerateDofs();
//     return space;
// }


std::shared_ptr<DGSpace> DGSpace::create(std::shared_ptr<Mesh> mesh, 
                                         Order order, 
                                         FieldType field_type) {

    if (order == Order::P2 && mesh->getOrder() == MeshOrder::Linear) {
        CFEM_INFO << "DGSpace: Auto-upgrading mesh to quadratic to support P2 space.";
        mesh->upgradeToQuadratic(); 
    }
    // Instantiation & Initialization
    auto space = std::make_shared<DGSpace>(ConstructorToken{}, mesh, order, static_cast<int>(field_type));
    // space->initReferenceElements();
    
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
        const auto& refEl = getRefElement(m_mesh->getCell(i).type);
        currentOffset += refEl.getNumNodes();
    }
    m_cellOffsets[numCells] = currentOffset; // Offset total final
    
    m_totalDofs = currentOffset * m_dimension;
}

CfemInt DGSpace::getPositionId(CfemInt cellId, int localNodeIdx) const {
    return m_cellOffsets[cellId] + localNodeIdx;
}

CfemInt DGSpace::getGlobalDofIdx(CfemInt cellId, int localNodeIdx, int component) const {
    CfemInt posId = getPositionId(cellId, localNodeIdx);
    CfemInt numUniqueNodes = m_totalDofs / m_dimension;
    
    // Keep the same BLOCK-wise logic: [All Comp 0][All Comp 1]...
    return posId + (component * numUniqueNodes);
}

int DGSpace::getDofsPerCell(CfemInt cellId) const {
    return getRefElement(m_mesh->getCell(cellId).type).getNumNodes() * m_dimension;
}

void DGSpace::getElementDofs(CfemInt cellId, std::vector<CfemInt>& globalIndices) const {
    int nNodes = getRefElement(m_mesh->getCell(cellId).type).getNumNodes();
    globalIndices.resize(nNodes * m_dimension);
    
    for (int c = 0; c < m_dimension; ++c) {
        for (int i = 0; i < nNodes; ++i) {
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
        int dofsSelf = getDofsPerCell(i);
        
        // 2. DOFs des voisins (en parcourant l'adjacence construite par buildAdjacency)
        int dofsNeighbors = 0;
        
        // On utilise directement les pointeurs initialisés par Mesh::buildAdjacency()
        // cell.num_faces correspond au nombre d'entrées dans le bloc CSR pour cette cellule
        for (int f = 0; f < cell.num_faces; ++f) {
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
        for (int comp = 0; comp < m_dimension; ++comp) {
            int nNodes = getRefElement(cell.type).getNumNodes();
            for (int n = 0; n < nNodes; ++n) {
                CfemInt globalIdx = getGlobalDofIdx(i, n, comp);
                nnz[globalIdx] = totalRowNNZ;
            }
        }
    }
}


// Math evaluations (evaluateBasis/evaluateGradients) are identical to Lagrange 
// because they both use the ReferenceElement interface!
void DGSpace::evaluateBasis(CfemInt cellId, const Vector3D& xi, std::vector<CfemReal>& values) const {
    ShapeInfo shape;
    getRefElement(m_mesh->getCell(cellId).type).evaluate(xi, ShapeEvaluationFlags::Value, shape);
    values.assign(shape.values.begin(), shape.values.end());
}

void DGSpace::evaluateGradients(CfemInt cellId, const Vector3D& xi, std::vector<Vector3D>& grads) const {
    static thread_local DynamicVector<Vector3D> physCoords;
    m_mesh->getCellVerticesCoords(cellId, physCoords);
    
    ShapeInfo shape;
    MappingInfo map;
    getRefElement(m_mesh->getCell(cellId).type).computeMapping(xi, physCoords, ShapeEvaluationFlags::Gradient, shape, map);
    grads.assign(map.physicalShapeGradients.begin(), map.physicalShapeGradients.end());
}

/**
 * @brief Récupère les DOFs associés à une entité géométrique (bord ou groupe physique).
 * @param internalEntityId L'identifiant (Gmsh tag) de l'entité.
 * @param entityDofs Vecteur à remplir avec les indices globaux.
 * * En DG, pour une entité de dimension (dim-1), on identifie les facettes appartenant 
 * à cette entité, puis on récupère les DOFs de la cellule adjacente (left_cell) 
 * qui sont localisés sur cette face.
 */
void DGSpace::getEntityDofs(int internalEntityId, std::vector<CfemInt>& entityDofs) const {
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
        const auto& refEl = getRefElement(cell.type);
        
        // 2. Trouver quels nœuds locaux de la cellule sont sur cette facette précise
        // On récupère la connectivité de la facette pour comparer
        auto facetNodes = m_mesh->getFacetConnectivity(fId);
        
        // On parcourt les faces locales de la cellule pour identifier celle qui correspond à fId
        for (int iFace = 0; iFace < cell.num_faces; ++iFace) {
            if (cell.facet_ids_ptr[iFace] == fId) {
                // On a trouvé la face locale correspondante !
                // On récupère les indices des nœuds locaux de la cellule appartenant à cette face
                auto localNodesOnFace = refEl.getFaceNodes(iFace);

                // 3. Ajouter les DOFs globaux pour chaque composante
                for (int comp = 0; comp < m_dimension; ++comp) {
                    for (int localIdx : localNodesOnFace) {
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