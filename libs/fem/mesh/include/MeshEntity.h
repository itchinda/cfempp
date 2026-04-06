#ifndef CFEM_MESHENTITY_
#define CFEM_MESHENTITY_

#include <vector>
#include <string>
#include <algorithm>
#include "cfemutils.h"


namespace cfem {

class Mesh;

class MeshEntity {
public:
    enum class Type { Vertex, Facet, Cell, Mixed };

    explicit MeshEntity(const Mesh* mesh, 
                        CfemInt id, 
                        const std::string& name, 
                        CfemInt dim);

    const Mesh& getMesh() const { return *m_mesh; }
    CfemInt getDimension() const { return m_dim; }
    CfemInt getId() const { return m_id; }
    const std::string& getName() const { return m_name; }

    void addVertex(CfemInt id) { m_vertex_ids.push_back(id); }
    void addFacet(CfemInt id)  { m_facet_ids.push_back(id); }
    void addCell(CfemInt id)   { m_cell_ids.push_back(id); }

    const std::vector<CfemInt>& getVertexIds() const { return m_vertex_ids; }
    const std::vector<CfemInt>& getFacetIds()  const { return m_facet_ids; }
    const std::vector<CfemInt>& getCellIds()   const { return m_cell_ids; }
    
    size_t getNumVertices() const { return m_vertex_ids.size(); }
    size_t getNumCells()    const { return m_cell_ids.size(); }
    size_t getNumFacets()   const { return m_facet_ids.size(); }

    /**
     * @brief Replaces the vertex list (used during mesh upgrade).
     * We use std::move to avoid expensive copies of large ID vectors.
     */
    void setVertexIds(std::vector<CfemInt>&& ids) { 
        m_vertex_ids = std::move(ids); 
        // We sort after setting to ensure binary search works later
        std::sort(m_vertex_ids.begin(), m_vertex_ids.end());
    }
    
    void setDimension(CfemInt p_dim){
        m_dim = p_dim;
    }

    bool containsCell(CfemInt cell_id) const {
        return std::binary_search(m_cell_ids.begin(), m_cell_ids.end(), cell_id);
    }

    void update();
    void updateUniverse();
    void updateBoundary();
    void finalize() ;
    void clear(){
        m_cell_ids.clear();
        m_facet_ids.clear();
        m_vertex_ids.clear();
    }

    std::string getEntityFormat(bool include_Ids = false) const;

    bool isEmpty() const { return m_cell_ids.empty() && m_facet_ids.empty() && m_vertex_ids.empty(); }
    bool containsCells() const { return !m_cell_ids.empty(); }
    bool containsFacets() const { return !m_facet_ids.empty(); }

    void reserveVertices(CfemInt size) { m_vertex_ids.reserve(size); }
    void reserveCells(CfemInt size) { m_facet_ids.reserve(size); }

    void replaceVertexList(const std::vector<CfemInt>& newVertexIds);
    void replaceVertexList(std::vector<CfemInt>&& newVertexIds);



private:
    const Mesh* m_mesh;
    CfemInt m_id;
    std::string m_name;
    CfemInt m_dim; 
    std::vector<CfemInt> m_cell_ids;
    std::vector<CfemInt> m_facet_ids;
    std::vector<CfemInt> m_vertex_ids;

    void processAndStore(std::vector<CfemInt>& vec) {
        if (!std::is_sorted(vec.begin(), vec.end())) {
            std::sort(vec.begin(), vec.end());
        }
        vec.erase(std::unique(vec.begin(), vec.end()), vec.end());
    }

};

} //namespace
#endif