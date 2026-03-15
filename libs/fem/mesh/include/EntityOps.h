#ifndef CFEM_ENTITY_OPS_H
#define CFEM_ENTITY_OPS_H

#include "Mesh.h"
#include "MeshEntity.h"
#include <algorithm>
#include <vector>
#include <string>

namespace cfem {
namespace EntityOps {

    /**
     * @brief Performs a boolean union of two mesh entities, saves the result in the mesh, 
     * and returns a pointer to the new entity.
     */
    inline MeshEntity* Union(Mesh& mesh, 
                             const MeshEntity& a, 
                             const MeshEntity& b, 
                             const std::string& new_name) 
    {
        CFEM_ASSERT(&a.getMesh()== &mesh && &b.getMesh()== &mesh, "EntityOps::Union Error: Entities must belong to the same parent mesh.");
        CfemInt new_dim = std::max(a.getDimension(), b.getDimension());
        CfemInt new_id = mesh.createEntity(new_name, new_dim);

        MeshEntity* result = mesh.getEntity(new_id); // Assuming this getter exists

        // 4. Helper lambda for linear-time O(N) set union on sorted vectors
        auto do_union = [](const std::vector<CfemInt>& v1, const std::vector<CfemInt>& v2, std::vector<CfemInt>& out) {
            out.reserve(v1.size() + v2.size());
            std::set_union(v1.begin(), v1.end(),
                           v2.begin(), v2.end(),
                           std::back_inserter(out));
        };

        // 5. Combine the sorted IDs
        std::vector<CfemInt> new_vertices, new_facets, new_cells;
        do_union(a.getVertexIds(), b.getVertexIds(), new_vertices);
        do_union(a.getFacetIds(), b.getFacetIds(), new_facets);
        do_union(a.getCellIds(), b.getCellIds(), new_cells);

        // 6. Inject the IDs using your fast move-semantics
        result->replaceVertexList(new_vertices);
        
        // (Assuming you add replaceFacetList and replaceCellList similarly to MeshEntity)
        for (auto id : new_facets) result->addFacet(id);
        for (auto id : new_cells)  result->addCell(id);

        // Ensure internal state is perfectly clean
        result->finalize();

        return result;
    }

} // namespace EntityOps
} // namespace cfem

#endif