#include "MeshEntity.h"
#include "Mesh.h"

#include <algorithm>

MeshEntity::MeshEntity(const Mesh* mesh, 
                       CfemInt id, 
                       const std::string& name, 
                       CfemInt dim):
    m_mesh(mesh), m_id(id), m_name(name), m_dim(dim)
{

}

std::string MeshEntity::getEntityFormat() const {
    std::stringstream ss;
    ss << "------------------------------------------\n";
    ss << "ENTITY INFO: " << m_name << " (ID: " << m_id << ")\n";
    ss << "Dimension:   " << m_dim << "\n";
    ss << "------------------------------------------\n";

    ss << "$Nodes\n";
    ss << m_vertex_ids.size() << " # total number of nodes\n";
    for (size_t i = 0; i < m_vertex_ids.size(); ++i) {
        ss << m_vertex_ids[i] << (i == m_vertex_ids.size() - 1 ? "" : ", ");
    }
    ss << "\n$EndNodes\n";

    if (!m_cell_ids.empty()) {
        ss << "$Cells\n";
        ss << m_cell_ids.size() << " # total number of cells\n";
        for (size_t i = 0; i < m_cell_ids.size(); ++i) {
            ss << m_cell_ids[i] << (i == m_cell_ids.size() - 1 ? "" : ", ");
        }
        ss << "\n$EndCells\n";
    }

    if (!m_facet_ids.size() == 0) {
        ss << "$Facets\n";
        ss << m_facet_ids.size() << " # total number of facets\n";
        for (size_t i = 0; i < m_facet_ids.size(); ++i) {
            ss << m_facet_ids[i] << (i == m_facet_ids.size() - 1 ? "" : ", ");
        }
        ss << "\n$EndFacets\n";
    }
    
    ss << "------------------------------------------\n";
    return ss.str();
}


/**
 * @brief Synchronizes the vertex list based on current cell and facet IDs.
 * Call this after the parent mesh has been upgraded to quadratic.
 */
void MeshEntity::update() {
    std::unordered_set<CfemInt> updated_nodes;

    for (CfemInt cell_id : m_cell_ids) {
        const auto& cell = m_mesh->getCell(cell_id);
        for (CfemInt v_id : cell.vertex_ids) {
            updated_nodes.insert(v_id);
        }
    }

    for (CfemInt facet_id : m_facet_ids) {
        const auto& facet = m_mesh->getFacet(facet_id);
        for (int i = 0; i < facet.number_of_vertices; ++i) {
            updated_nodes.insert(facet.vertex_ids[i]);
        }
    }

    m_vertex_ids.assign(updated_nodes.begin(), updated_nodes.end());
    this->finalize(); 
}

void MeshEntity::updateUniverse() {
    m_vertex_ids.resize(m_mesh->getNumVertices());
    for (CfemInt i = 0; i < (CfemInt)m_mesh->getNumVertices(); ++i) m_vertex_ids[i] = i;

    m_cell_ids.resize(m_mesh->getNumCells());
    for (CfemInt i = 0; i < (CfemInt)m_mesh->getNumCells(); ++i) m_cell_ids[i] = i;

    m_facet_ids.resize(m_mesh->getNumFacets());
    for (CfemInt i = 0; i < (CfemInt)m_mesh->getNumFacets(); ++i) m_facet_ids[i] = i;

    // No need to finalize (sort or unique) since we inserted indices in order
}

void MeshEntity::updateBoundary() {
    std::unordered_set<CfemInt> boundary_nodes;
    m_facet_ids.clear();

    // Loop through all facets in the mesh
    for (CfemInt i = 0; i < (CfemInt)m_mesh->getNumFacets(); ++i) {
        const auto& facet = m_mesh->getFacet(i);
        
        if (facet.isBoundary()) {
            m_facet_ids.push_back(i);
            for (int n = 0; n < facet.number_of_vertices; ++n) {
                boundary_nodes.insert(facet.vertex_ids[n]);
            }
        }
    }

    m_vertex_ids.assign(boundary_nodes.begin(), boundary_nodes.end());
    this->finalize();
}

/**
 * @brief Ensures all ID vectors are sorted and contain no duplicates.
 * Should be called after initial construction or significant changes.
 */
void MeshEntity::finalize() {
    auto clean = [](std::vector<CfemInt>& vec) {
        std::sort(vec.begin(), vec.end());
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    };
    clean(m_vertex_ids);
    clean(m_facet_ids);
    clean(m_cell_ids);
}