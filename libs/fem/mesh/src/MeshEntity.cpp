#include "MeshEntity.h"
#include "Mesh.h"

#include <algorithm>

namespace cfem
{

    MeshEntity::MeshEntity(const Mesh *mesh,
                           CfemInt id,
                           const std::string &name,
                           CfemInt dim) : m_mesh(mesh), m_id(id), m_name(name), m_dim(dim)
    {}

    std::string MeshEntity::getEntityFormat(bool include_Ids) const
    {
        std::stringstream ss;
        ss << "------------------------------------------\n";
        ss << "ENTITY INFO: " << m_name << " (ID: " << m_id << ")\n";
        ss << "Dimension:   " << m_dim << "\n";
        ss << "------------------------------------------\n";

        ss << "$Nodes\n";
        ss << m_vertex_ids.size() << " # total number of nodes\n";
        if (include_Ids){
            for (size_t i = 0; i < m_vertex_ids.size(); ++i){
                ss << m_vertex_ids[i] << (i == m_vertex_ids.size() - 1 ? "" : ", ");
            }
        }
        ss << "\n$EndNodes\n";

        if (!m_cell_ids.empty()){
            ss << "$Cells\n";
            ss << m_cell_ids.size() << " # total number of cells\n";
            if (include_Ids){
                for (size_t i = 0; i < m_cell_ids.size(); ++i){
                    ss << m_cell_ids[i] << (i == m_cell_ids.size() - 1 ? "" : ", ");
                }
            }
            ss << "\n$EndCells\n";
        }

        if (!m_facet_ids.size() == 0){
            ss << "$Facets\n";
            ss << m_facet_ids.size() << " # total number of facets\n";
            if (include_Ids){
                for (size_t i = 0; i < m_facet_ids.size(); ++i){
                    ss << m_facet_ids[i] << (i == m_facet_ids.size() - 1 ? "" : ", ");
                }
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
    void MeshEntity::update()
    {
        m_vertex_ids.clear();

        // We use a temporary vector of flags to track which vertices are used by this entity
        std::vector<uint8_t> node_marker(m_mesh->getNumVertices(), 0);

        for (CfemInt cell_id : m_cell_ids){
            // Use the accessor to reach into the Big Block
            auto c_vids = m_mesh->getCellConnectivity(cell_id);
            for (CfemInt v_id : c_vids){
                node_marker[v_id] = 1;
            }
        }

        for (CfemInt facet_id : m_facet_ids){
            auto f_vids = m_mesh->getFacetConnectivity(facet_id);
            for (CfemInt v_id : f_vids){
                node_marker[v_id] = 1;
            }
        }

        // This naturally keeps node IDs sorted, making finalize() much cheaper
        for (CfemInt v_id = 0; v_id < (CfemInt)node_marker.size(); ++v_id){
            if (node_marker[v_id]){
                m_vertex_ids.push_back(v_id);
            }
        }

        this->finalize();
    }

    void MeshEntity::updateUniverse()
    {
        m_vertex_ids.resize(m_mesh->getNumVertices());
        for (CfemInt i = 0; i < (CfemInt)m_mesh->getNumVertices(); ++i)
            m_vertex_ids[i] = i;

        m_cell_ids.resize(m_mesh->getNumCells());
        for (CfemInt i = 0; i < (CfemInt)m_mesh->getNumCells(); ++i)
            m_cell_ids[i] = i;

        m_facet_ids.resize(m_mesh->getNumFacets());
        for (CfemInt i = 0; i < (CfemInt)m_mesh->getNumFacets(); ++i)
            m_facet_ids[i] = i;

        // No need to finalize (sort or unique) since we inserted indices in order
    }

    void MeshEntity::updateBoundary()
    {
        m_facet_ids.clear();
        m_vertex_ids.clear();

        std::vector<bool> node_marker(m_mesh->getNumVertices(), false);

        CfemInt numFacets = static_cast<CfemInt>(m_mesh->getNumFacets());
        for (CfemInt i = 0; i < numFacets; ++i)
        {
            const auto &facet = m_mesh->getFacet(i);

            // Logical check for boundary (RightCell == -1)
            if (facet.isBoundary()){
                m_facet_ids.push_back(i);

                auto f_vids = m_mesh->getFacetConnectivity(i);

                for (CfemInt v_id : f_vids){
                    node_marker[v_id] = true;
                }
            }
        }

        for (CfemInt v_id = 0; v_id < (CfemInt)node_marker.size(); ++v_id)
        {
            if (node_marker[v_id]){
                m_vertex_ids.push_back(v_id);
            }
        }

        this->finalize();
    }

    /**
     * @brief Ensures all ID vectors are sorted and contain no duplicates.
     * Should be called after initial construction or significant changes.
     */
    void MeshEntity::finalize()
    {
        auto clean = [](std::vector<CfemInt> &vec){
            std::sort(vec.begin(), vec.end());
            vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
        };
        clean(m_vertex_ids);
        clean(m_facet_ids);
        clean(m_cell_ids);
    }

    /**
     * @brief Replaces the vertex list by copying an existing vector.
     *        Use this version when you need to keep the original source vector intact
     *        after the call. It performs a deep copy of the input data.
     * @param newVertexIds A reference to the vector containing the new Node IDs.
     * @note Sorting and duplicate removal are performed in-place after the copy.
     * @complexity O(N) memory allocation + O(N log N) worst-case sorting.
     */
    void MeshEntity::replaceVertexList(const std::vector<CfemInt> &newVertexIds)
    {
        if (newVertexIds.empty()){
            m_vertex_ids.clear();
            return;
        }

        m_vertex_ids = newVertexIds;
        processAndStore(m_vertex_ids);
    }

    /**
     * @brief Replaces the vertex list by transferring ownership (Move semantics).
     * * This is the high-performance (HPC) version. Use it with temporary vectors
     * (r-values) or by using std::move() on an existing vector that is no longer needed.
     * It avoids expensive memory allocations by "stealing" the source vector's data.
     * @param newVertexIds An r-value reference to the source vector.
     * @complexity O(1) memory transfer + O(N log N) worst-case sorting.
     */
    void MeshEntity::replaceVertexList(std::vector<CfemInt> &&newVertexIds)
    {
        if (newVertexIds.empty()){
            m_vertex_ids.clear();
            return;
        }

        m_vertex_ids = std::move(newVertexIds);
        processAndStore(m_vertex_ids);
    }

} // namespace