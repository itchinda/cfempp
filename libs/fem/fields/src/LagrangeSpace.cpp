#include "LagrangeSpace.h"
#include "ReferenceSegment.h"
#include "ReferenceTriangle.h"
#include "ReferenceQuad.h"
#include "ReferenceTetrahedron.h"
#include "ReferencePrism.h"
#include "ReferenceHexahedron.h"

#include "FEField.h"

#include <set>
#include <algorithm>

namespace cfem {

LagrangeSpace::LagrangeSpace(ConstructorToken, std::shared_ptr<Mesh> mesh, Order order, int dimension)
    : FESpace(mesh, order, Continuity::CG, dimension) 
{
}

void LagrangeSpace::initReferenceElements() {
    m_refElements.clear();

    if (m_mesh->getDimension() == 1){
        m_refElements[CellType::LINE2] = std::make_unique<P1ReferenceSegment>();
        m_refElements[CellType::LINE3] = std::make_unique<P2ReferenceSegment>();
    }
    else if (m_mesh->getDimension() == 2) {
        m_refElements[CellType::TRI3]  = std::make_unique<P1ReferenceTriangle>();
        m_refElements[CellType::QUAD4] = std::make_unique<Q1ReferenceQuad>();
        m_refElements[CellType::TRI6]  = std::make_unique<P2ReferenceTriangle>();
        m_refElements[CellType::QUAD9] = std::make_unique<Q2ReferenceQuad>();
    }
    else if (m_mesh->getDimension() == 3) {
        m_refElements[CellType::TET4]    = std::make_unique<P1ReferenceTetrahedron>();
        m_refElements[CellType::PRISM6]  = std::make_unique<P1ReferencePrism>();
        m_refElements[CellType::HEX8]    = std::make_unique<Q1ReferenceHexahedron>();
        m_refElements[CellType::TET10]   = std::make_unique<P2ReferenceTetrahedron>();
        m_refElements[CellType::PRISM18] = std::make_unique<P2ReferencePrism>();
        m_refElements[CellType::HEX27]   = std::make_unique<Q2ReferenceHexahedron>();
    }

}

/**
 * @brief Implémentation de la logique de filtrage P1/P2.
 */
int LagrangeSpace::getNodesPerCell(CellType type) const {
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
    return 0;
}

void LagrangeSpace::enumerateDofs() {
    const CfemInt numMeshNodes = static_cast<CfemInt>(m_mesh->getNumVertices());
    m_nodeToLocalIdx.assign(numMeshNodes, -1);
    
    CfemInt counter = 0;

    // On boucle sur les cellules pour identifier les nœuds actifs selon l'ordre
    for (CfemInt c = 0; c < (CfemInt)m_mesh->getNumCells(); ++c) {
        const auto& cell = m_mesh->getCell(c);
        int nodesToUse = getNodesPerCell(cell.type);
        auto cell_ids = m_mesh->getCellConnectivity(cell.id);
        
        for (int n = 0; n < nodesToUse; ++n) {
            CfemInt vId = cell_ids[n];
            if (m_nodeToLocalIdx[vId] == -1) {
                m_nodeToLocalIdx[vId] = counter++;
            }
        }
    }

    m_numActiveNodes = counter;
    m_totalDofs = m_numActiveNodes * m_dimension;

    m_localToMeshNodeIdx.assign(m_numActiveNodes, -1);
    for (CfemInt meshId = 0; meshId < (CfemInt)m_nodeToLocalIdx.size(); ++meshId) {
        CfemInt local = m_nodeToLocalIdx[meshId];
        if (local != -1) {
            m_localToMeshNodeIdx[local] = meshId;
        }
    }
}

int LagrangeSpace::getDofsPerCell(CfemInt cellId) const {
    return getNodesPerCell(m_mesh->getCell(cellId).type) * m_dimension;
}

/**
 * @brief Remplit les DOFs avec la numérotation par bloc :
 * GlobalID = localNodeIdx + (compIdx * m_numActiveNodes)
 */
void LagrangeSpace::getElementDofs(CfemInt cellId, std::vector<CfemInt>& globalIndices) const {
    const auto& cell = m_mesh->getCell(cellId);
    int nodesInCell = getNodesPerCell(cell.type);
    auto cell_ids = m_mesh->getCellConnectivity(cell.id);
    
    globalIndices.resize(nodesInCell * m_dimension);

    for (int n = 0; n < nodesInCell; ++n) {
        CfemInt nodeIdx = m_nodeToLocalIdx[cell_ids[n]];
        
        for (int d = 0; d < m_dimension; ++d) {
            // Bloc 0 : Composante 0, Bloc 1 : Composante 1...
            globalIndices[n + d * nodesInCell] = nodeIdx + (d * m_numActiveNodes);
        }
    }
}

/**
 * @brief Calcul du NNZ (Sparsity Pattern) par blocs.
 */
void LagrangeSpace::computeSparsityPattern(std::vector<CfemInt>& nnz) const {
    nnz.assign(m_totalDofs, 0);

    for (CfemInt i = 0; i < (CfemInt)m_mesh->getNumVertices(); ++i) {
        CfemInt nodeIdx = m_nodeToLocalIdx[i];
        if (nodeIdx == -1) continue;

        const auto& vertex = m_mesh->getVertex(i);
        std::set<CfemInt> neighborNodes;

        // On cherche les nœuds connectés
        for (CfemInt cellId : vertex.getConnectedCells()){
            const auto& cell = m_mesh->getCell(cellId);
            auto cell_ids = m_mesh->getCellConnectivity(cell.id);

            int nInCell = getNodesPerCell(cell.type);
            for (int n = 0; n < nInCell; ++n) {
                neighborNodes.insert(m_nodeToLocalIdx[cell_ids[n]]);
            }
        }

        // Nombre de connexions : Nb de nœuds voisins * Nb de composantes
        CfemInt rowNNZ = static_cast<CfemInt>(neighborNodes.size()) * m_dimension;

        // On reporte cette valeur pour chaque bloc (chaque composante du nœud)
        for (int d = 0; d < m_dimension; ++d) {
            nnz[nodeIdx + (d * m_numActiveNodes)] = rowNNZ;
        }
    }
}

void LagrangeSpace::getBoundaryDofs(int internalEntityId, std::vector<CfemInt>& boundaryDofs) const {
    boundaryDofs.clear();
    const auto& entities = m_mesh->getEntities();
    
    for (const auto& ent : entities) {
        if (ent->getId() == internalEntityId) {
            for (CfemInt vId : ent->getVertexIds()) {
                CfemInt nodeIdx = m_nodeToLocalIdx[vId];
                if (nodeIdx != -1) {
                    for (int d = 0; d < m_dimension; ++d) {
                        boundaryDofs.push_back(nodeIdx + (d * m_numActiveNodes));
                    }
                }
            }
            break;
        }
    }
}

// -----------------------------------------------------------------------------
// Helper pour récupérer l'élément de référence (cache)
// -----------------------------------------------------------------------------
const ReferenceElement& LagrangeSpace::getRefElement(CellType type) const {
    CFEM_ASSERT(!m_refElements.empty());
    auto it = m_refElements.find(type);
    if (it == m_refElements.end()) {
        CFEM_ERROR << "LagrangeSpace: Element type not supported or not initialized.";
    }
    return *it->second;
}

// -----------------------------------------------------------------------------
// Évaluation des Fonctions de Forme (Basis Functions)
// -----------------------------------------------------------------------------
void LagrangeSpace::evaluateBasis(CfemInt cellId, const Vector3D& xi, std::vector<CfemReal>& values) const {
    // 1. Récupération de la cellule et de son élément de référence
    const Cell& cell = m_mesh->getCell(cellId);
    const ReferenceElement& refEl = getRefElement(cell.type);

    // 2. Préparation du buffer léger (ShapeInfo)
    // On le déclare en local pour garantir la thread-safety
    ShapeInfo shape;
    
    // 3. Calcul des fonctions de forme N(xi)
    // On demande uniquement les VALEURS pour l'instant
    refEl.evaluate(xi, ShapeEvaluationFlags::Value, shape);

    // 4. Copie vers le vecteur de sortie std::vector
    // On redimensionne si nécessaire (optimisation possible: réserver avant)
    if (values.size() != shape.values.size()) {
        values.resize(shape.values.size());
    }
    
    for (size_t i = 0; i < shape.values.size(); ++i) {
        values[i] = shape.values[i];
    }
}

// -----------------------------------------------------------------------------
// Évaluation des Gradients Physiques (Physical Gradients)
// -----------------------------------------------------------------------------
void LagrangeSpace::evaluateGradients(CfemInt cellId, const Vector3D& xi, std::vector<Vector3D>& grads) const {
    // 1. Récupération de la cellule et de l'élément de référence
    const Cell& cell = m_mesh->getCell(cellId);
    const ReferenceElement& refEl = getRefElement(cell.type);

    // 2. Extraction des coordonnées PHYSIQUES des nœuds
    // Comme ton maillage stocke tous les nœuds (P1 ou P2), on prend tout.
    const auto& nodeIndices = m_mesh->getCellConnectivity( cellId);
    size_t numNodes = nodeIndices.size();

    // On utilise DynamicVector car ReferenceElement l'attend
    DynamicVector<Vector3D> physicalNodes(numNodes);

    for (size_t i = 0; i < numNodes; ++i) {
        // On récupère la coordonnée réelle (Vertex -> Vector3D)
        physicalNodes[i] = m_mesh->getVertex(nodeIndices[i]).coords;
    }

    // 3. Préparation des structures de calcul
    ShapeInfo shape;
    MappingInfo map;

    // 4. Calcul du Mapping Isoparamétrique complet
    // computeMapping fait tout : evaluate(Grad), Jacobian, invJacobian, push-forward
    refEl.computeMapping(xi, 
                         physicalNodes, 
                         ShapeEvaluationFlags::Gradient, // On veut dN/dx
                         shape, 
                         map);

    // 5. Transfert des résultats
    if (grads.size() != map.physicalShapeGradients.size()) {
        grads.resize(map.physicalShapeGradients.size());
    }

    for (size_t i = 0; i < map.physicalShapeGradients.size(); ++i) {
        grads[i] = map.physicalShapeGradients[i];
    }
}

} // namespace cfem