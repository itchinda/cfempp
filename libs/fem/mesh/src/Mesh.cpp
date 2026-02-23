#include "Mesh.h"
#include "FESpace.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <iomanip>
#include <sstream>

#include "QuadMeshBuilder.h"
#include "LineMeshBuilder.h"
#include "DiskMeshBuilder.h"
#include "BlockMeshBuilder.h"
#include "SphereMeshBuilder.h"
#include "CylinderMeshBuilder.h"
#include "AnnulusMeshBuilder.h"
#include "HollowCylinderMeshBuilder.h"
#include "MeshConformityChecker.h"

namespace cfem
{

    Mesh::Mesh(ConstructorToken ltok)
    {
        CFEM_CHECK_INVARIANTS();
    }

    Mesh::Mesh(ConstructorToken lTok, CfemInt p_meshDim, MeshOrder p_meshOrder) : m_meshDim(p_meshDim), m_meshOrder(p_meshOrder)
    {
        CFEM_PRECONDITION(p_meshDim >= 1 && p_meshDim <= 3);
        CFEM_POSTCONDITION(m_meshDim == p_meshDim);
        CFEM_POSTCONDITION(m_meshOrder == p_meshOrder);
        CFEM_CHECK_INVARIANTS();
    }

    CfemInt Mesh::addVertex(CfemReal x, CfemReal y, CfemReal z)
    {
        CfemInt id = m_vertices.size();
        m_vertices.push_back({{x, y, z}, {}});

        CFEM_POSTCONDITION(static_cast<CfemInt>(m_vertices.size()) > id);
        return id;
    }

    CfemInt Mesh::addCell(CellType type, const std::vector<CfemInt> &v_ids)
    {
        // Logic Validation
        checkConsistency(type, v_ids); // Could also check, avoiding adding 2D element in 3D mesh is essential

#ifdef NDEBUG
        const auto& final_vids = v_ids;
#else
        // In Debug, we may need to modify a copy
        std::vector<CfemInt> final_vids = v_ids; 
        CfemReal measure = computeCellMeasure(type, final_vids);

        if (measure < -cfem::ZERO_TOLERANCE) 
        {
            fixElementOrientation(type, final_vids);
            measure = computeCellMeasure(type, final_vids);
            if (measure < -cfem::ZERO_TOLERANCE){
                CFEM_ERROR("Orientation fix failed. Element might be mangled.");
            }
            m_invertedCellCounts[type]++;
        }
        else if (std::abs(measure) <= cfem::ZERO_TOLERANCE)
        {
            CFEM_WARNING("Degenerate (zero area/volume) cell added to mesh.");
        }
#endif

        // --- Critical part: the usage must be valid for both paths ---
        CfemInt id = static_cast<CfemInt>(m_cells.size());
        size_t startPos = m_allCellVertexIds.size();
        
        // Copy the potentially fixed nodes into the big block
        m_allCellVertexIds.insert(m_allCellVertexIds.end(), final_vids.begin(), final_vids.end());
        
        // 4. Update the CSR-style offsets vector (for external CSR algorithms)
        if(m_cellOffsets.empty()) m_cellOffsets.push_back(0);
        m_cellOffsets.push_back(m_allCellVertexIds.size());

        // 5. Initialize the Cell object with indices, NOT pointers
        m_cells.emplace_back();
        auto& cell = m_cells.back();
        cell.id = id;
        cell.type = type;
        cell.offset = startPos;        // Store the integer index
        cell.node_count = final_vids.size(); // Store the integer count

        // 6. Manage unique element types
        int typeIdx = static_cast<int>(type);
        if (!m_typeFlags.test(typeIdx))
        {
            m_typeFlags.set(typeIdx);
            m_uniqueTypes.push_back(type);
        }

        return id;
    }

    void Mesh::checkConsistency(CellType type, const std::vector<CfemInt> &v_ids) const
    {
        MeshOrder cellOrder = getCellOrderFromType(type);
        if (cellOrder != m_meshOrder)
        {
            CFEM_ERROR("Cannot add " + cellTypeToString(type) +
                       " to a " + (m_meshOrder == MeshOrder::Linear ? "linear mesh." : m_meshOrder == MeshOrder::Quadratic ? "quadratic mesh."
                                                                                                                           : "mesh with unkmown order."));
        }

        CfemInt cellDim = getDimFromType(type);
        if (cellDim != m_meshDim)
        {
            CFEM_ERROR("Cannot add a " + std::to_string(cellDim) + "D element (" +
                       cellTypeToString(type) + ") to a " + std::to_string(m_meshDim) + "D mesh.");
        }

        if (v_ids.size() != getNodeCountPerCell(type))
        {
            CFEM_ERROR("Topology error: " + cellTypeToString(type) +
                       " expects " + std::to_string(getNodeCountPerCell(type)) +
                       " nodes, but " + std::to_string(v_ids.size()) + " were provided.");
        }

        for (CfemInt v_id : v_ids)
        {
            if (v_id < 0 || v_id >= static_cast<CfemInt>(m_vertices.size()))
            {
                CFEM_ERROR("Index out of bounds: vertex " + std::to_string(v_id) +
                           " does not exist in a mesh with " +
                           std::to_string(m_vertices.size()) + " vertices.");
            }
        }
    }

    void Mesh::fixElementOrientation(CellType type, std::vector<CfemInt> &vids)
    {
        switch (type)
        {
            case CellType::TRI3:
                std::swap(vids[0], vids[1]);
                break;
            case CellType::TRI6:
                std::swap(vids[0], vids[1]); // Corners
                std::swap(vids[4], vids[5]); // Mid-nodes on edges 0-1 and 0-2
                break;
            case CellType::QUAD4:
                std::swap(vids[1], vids[3]);
                break;
            case CellType::QUAD9:
                std::swap(vids[1], vids[3]); // Corners
                std::swap(vids[4], vids[7]); // Mid-nodes on edges
                std::swap(vids[5], vids[6]); 
                break;
            case CellType::TET4:
                std::swap(vids[0], vids[1]); // Corners
                break;
            case CellType::TET10:
                std::swap(vids[0], vids[1]); // Swap corners 0 and 1
                std::swap(vids[5], vids[6]); // Swap mid-nodes on edges 1-2 and 2-0
                std::swap(vids[8], vids[7]); // Swap mid-nodes on edges 1-3 and 0-3
                break;
            case CellType::PRISM6:
            case CellType::HEX8:{
                int half = getNodeCountPerCell(type) / 2;
                for(int i=0; i<half; ++i) std::swap(vids[i], vids[i+half]);
                break;
            }
            case CellType::HEX27:{
                // Swap Corners: Bottom (0-3) with Top (4-7)
                for (int i = 0; i < 4; ++i) {
                    std::swap(vids[i], vids[i + 4]);
                }

                // Swap Horizontal Edges: Bottom (8-11) with Top (12-15)
                for (int i = 8; i < 12; ++i) {
                    std::swap(vids[i], vids[i + 4]);
                }

                // Swap Face Centers: Bottom (24) with Top (25)
                std::swap(vids[24], vids[25]);
                // - Vertical Edges (16-19): They just flip upside down in place.
                // - Side Face Centers (20, 21, 22, 23): Dead center of the side walls.
                // - Body Center (26): Dead center of the volume.
                break;
            }
            case CellType::PRISM18: {
                // 1. Swap Corners: Bottom (0,1,2) with Top (3,4,5)
                for (int i = 0; i < 3; ++i) std::swap(vids[i], vids[i + 3]);

                // 2. Swap Horizontal Edges: Bottom (6,7,8) with Top (9,10,11)
                for (int i = 6; i < 9; ++i) std::swap(vids[i], vids[i + 3]);
                // Vertical edges 12,13,14 and face centers 15,16,17 stay put.
                break;
            }
            default: CFEM_ERROR("Unsupported cell type.");
        }
    }

    /**
     * @brief Establishes Vertex -> Cell connectivity using a Compressed Sparse Row (CSR) pattern.
     * This method builds an "Inverse Map" that allows for O(1) access to all cells sharing a
     * specific vertex. By using a single contiguous memory block (@ref m_globalConnectivity),
     * it avoids the memory fragmentation and allocation overhead associated with
     * storing a std::vector inside every Vertex object.
     * @note This is critical for large meshes (e.g., 8 million+ cells) to maintain cache locality.
     * @par Usage Example:
     * @code
     * for (CfemInt v_id : cell.vertex_ids) {
     *      const auto& vertex = m_vertices[v_id];
     *      // Looping over neighbors is a simple pointer increment (high cache hits)
     *      for (CfemInt i = 0; i < vertex.num_connected_cells; ++i) {
     *          CfemInt neighbor_id = vertex.connected_cells_ptr[i];
     *          // ... process neighbor_id ...
     *      }
     *  }
     * @endcode
     */
    void Mesh::buildInverseConnectivities()
    {
        PROFILE_SCOPE("Mesh::buildInverseConnectivities");

        CfemInt numVertices = static_cast<CfemInt>(m_vertices.size());
        CfemInt numCells = static_cast<CfemInt>(m_cells.size());

        // 1. Count connections per vertex
        // Using a flat vector for counts is efficient for 10^6+ vertices
        std::vector<CfemInt> counts(numVertices, 0);
        CfemInt totalEntries = 0;

        for (CfemInt c_id = 0; c_id < numCells; ++c_id)
        {
            // Use the new accessor to get the safe span of node IDs
            auto v_ids = getCellConnectivity(c_id); 
            for (CfemInt v_id : v_ids)
            {
                counts[v_id]++;
                totalEntries++;
            }
        }

        // 2. Compute Offsets (Prefix Sum)
        // offsets[v_id] tells us exactly where in m_globalConnectivity the cell list for v_id starts
        std::vector<CfemInt> offsets(numVertices + 1, 0);
        for (CfemInt i = 0; i < numVertices; ++i)
        {
            offsets[i + 1] = offsets[i] + counts[i];
        }

        // 3. Allocate the "Big Block" 
        // This stores all cell IDs connected to all vertices in one contiguous array
        m_globalConnectivity.assign(totalEntries, -1);

        // 4. Fill the Big Block
        // Reset tracker to zero so we can use it as an offset within each vertex's chunk
        std::fill(counts.begin(), counts.end(), 0);
        for (CfemInt c_id = 0; c_id < numCells; ++c_id)
        {
            auto v_ids = getCellConnectivity(c_id);
            for (CfemInt v_id : v_ids)
            {
                CfemInt writePos = offsets[v_id] + counts[v_id];
                m_globalConnectivity[writePos] = c_id;
                counts[v_id]++;
            }
        }

        // 5. Update Vertex pointers
        // This part is actually safe even if m_globalConnectivity reallocates, 
        // BUT only as long as you call this function AFTER all cells are loaded.
        for (CfemInt v_id = 0; v_id < numVertices; ++v_id)
        {
            m_vertices[v_id].connected_cells_ptr = &m_globalConnectivity[offsets[v_id]];
            m_vertices[v_id].num_connected_cells = static_cast<CfemInt>(offsets[v_id + 1] - offsets[v_id]);
        }
    }

   void Mesh::buildAdjacency() 
{
    PROFILE_SCOPE("Mesh::buildAdjacency");

    m_allCellNeighbors.clear();
    m_allCellFacetIds.clear();
    m_neighborOffsets.clear();

    const size_t numCells = m_cells.size();
    m_neighborOffsets.assign(numCells + 1, 0);
    
    size_t totalSlots = 0;
    for (size_t i = 0; i < numCells; ++i) {
        m_neighborOffsets[i] = totalSlots;
        totalSlots += getLocalFacesCount(m_cells[i].type);
    }
    m_neighborOffsets[numCells] = totalSlots; 

    m_allCellNeighbors.assign(totalSlots, -1);
    m_allCellFacetIds.assign(totalSlots, -1);

    for (const auto& f : m_facets) 
    {
        auto f_nodes = getFacetConnectivity(f.id);

        // --- UNIVERSAL 2D/3D MATCHER ---
        auto getLocalFaceIdx = [&](CfemInt cell_id) -> int {
            int num_faces = getLocalFacesCount(m_cells[cell_id].type);
            
            for (int i = 0; i < num_faces; ++i) {
                // Use your existing helper that knows the HEX8/QUAD4 face patterns
                LocalFace currentLocal = getLocalFaceData(cell_id, i);
                
                // Compare the set of nodes (order independent)
                if (static_cast<size_t>(currentLocal.size) == f_nodes.size()) {
                    bool match = true;
                    for(size_t fn : f_nodes) {
                        bool found = false;
                        for(int k=0; k < currentLocal.size; ++k) {
                            if(currentLocal.v[k] == (CfemInt)fn) { found = true; break; }
                        }
                        if(!found) { match = false; break; }
                    }
                    if(match) return i;
                }
            }
            return -1;
        };

        // Link Left Cell
        int localIdxL = getLocalFaceIdx(f.left_cell);
        if (localIdxL != -1) {
            m_allCellNeighbors[m_neighborOffsets[f.left_cell] + localIdxL] = f.right_cell;
            m_allCellFacetIds[m_neighborOffsets[f.left_cell] + localIdxL] = f.id;
        }

        // Link Right Cell
        if (!f.isBoundary()) {
            int localIdxR = getLocalFaceIdx(f.right_cell);
            if (localIdxR != -1) {
                m_allCellNeighbors[m_neighborOffsets[f.right_cell] + localIdxR] = f.left_cell;
                m_allCellFacetIds[m_neighborOffsets[f.right_cell] + localIdxR] = f.id;
            }
        }
    }
}

    void Mesh::buildFacets() 
    {
        PROFILE_SCOPE("Mesh::buildFacets");

        // ---------------------------------------------------------
        // STEP 1: RESET & RESERVE
        // ---------------------------------------------------------
        m_facets.clear();
        m_allFacetVertexIds.clear();
        // m_facetOffsets is actually redundant if we store offset/count in the Facet struct,
        // but it's useful for CSR-style algorithms, so we'll keep it.
        m_facetOffsets.clear();

        size_t est_facets = m_cells.size() * 4; 
        m_allFacetVertexIds.reserve(est_facets * 4); 
        m_facetOffsets.reserve(est_facets + 1);
        m_facetOffsets.push_back(0); 

        std::vector<CfemInt> cell_face_to_facet(m_cells.size() * 6, -1);

        // Metadata structure to hold temp info
        struct FacetMeta {
            CfemInt left;
            CfemInt right;
            CellType type;
            size_t offset;
            size_t count;
        };
        std::vector<FacetMeta> temp_meta;
        temp_meta.reserve(est_facets);

        // ---------------------------------------------------------
        // STEP 2: BUILD TOPOLOGY
        // ---------------------------------------------------------
        for (CfemInt c_id = 0; c_id < static_cast<CfemInt>(m_cells.size()); ++c_id) 
        {

            const int num_faces = getLocalFacesCount(m_cells[c_id].type);
            
            for (int f_idx = 0; f_idx < num_faces; ++f_idx) 
            {
                if (cell_face_to_facet[c_id * 6 + f_idx] != -1) continue;

                LocalFace current_face = getLocalFaceData(c_id, f_idx);
                
                // --- Neighbor Search ---
                CfemInt neighbor_id = -1;
                int neighbor_f_idx = -1;

                CfemInt v_anchor = current_face.v[0];
                const auto& anchor_vert = m_vertices[v_anchor];

                for (CfemInt i = 0; i < anchor_vert.num_connected_cells; ++i) 
                {
                    CfemInt other_c_id = anchor_vert.connected_cells_ptr[i];
                    if (other_c_id == c_id) continue;

                    const int other_num_faces = getLocalFacesCount(m_cells[other_c_id].type);
                    for (int j = 0; j < other_num_faces; ++j) {
                        LocalFace other_face = getLocalFaceData(other_c_id, j);
                        if (current_face.matches(other_face)) {
                            neighbor_id = other_c_id;
                            neighbor_f_idx = j;
                            goto found_neighbor; 
                        }
                    }
                }
                found_neighbor:;

                // --- Store Connectivity in Big Block ---
                size_t start_offset = m_allFacetVertexIds.size();
                for (int k = 0; k < current_face.size; ++k) {
                    m_allFacetVertexIds.push_back(current_face.v[k]);
                }

                // --- Record Metadata ---
                temp_meta.push_back({
                    c_id, 
                    neighbor_id, 
                    current_face.type, 
                    start_offset, 
                    (size_t)current_face.size
                });

                // Update mapping
                CfemInt new_facet_id = static_cast<CfemInt>(temp_meta.size() - 1);
                cell_face_to_facet[c_id * 6 + f_idx] = new_facet_id;
                if (neighbor_id != -1) {
                    cell_face_to_facet[neighbor_id * 6 + neighbor_f_idx] = new_facet_id;
                }

                // Sync the CSR offsets (optional but good practice)
                m_facetOffsets.push_back(m_allFacetVertexIds.size());
            }
        }

        // ---------------------------------------------------------
        // STEP 3: FINALIZE FACETS (INDEX-BASED)
        // ---------------------------------------------------------
        m_facets.resize(temp_meta.size());


        for (size_t i = 0; i < temp_meta.size(); ++i) 
        {
            Facet& f = m_facets[i];
            const auto& meta = temp_meta[i];

            f.id = static_cast<CfemInt>(i);
            f.type = meta.type;
            f.left_cell = meta.left;
            f.right_cell = meta.right;
            
            // No pointers! Just integers. 
            // This is now 100% safe from reallocations.
            f.offset = meta.offset;
            f.node_count = meta.count;
            
        }
    }

    CfemInt Mesh::createEntity(const std::string &name, CfemInt dim)
    {
        validateEntityName(name);

        CfemInt newId = static_cast<CfemInt>(m_entities.size());
        if (getEntity(cfem::SYS_BOUNDARY_ID) && getEntity(cfem::SYS_BOUNDARY_ID))
        {
            newId -= 2;
        }
        else if (getEntity(cfem::SYS_BOUNDARY_ID) || getEntity(cfem::SYS_BOUNDARY_ID))
        {
            newId -= 1;
        }

        // 2. Double check just in case something weird happened to m_entity_ids
        while (m_entity_ids.find(newId) != m_entity_ids.end())
        {
            newId++;
        }

        m_entity_ids.insert(newId);
        m_entities.push_back(std::make_unique<MeshEntity>(this, newId, name, dim));

        return newId;
    }

    MeshEntity *Mesh::getEntity(const std::string &name)
    {
        for (auto &entity_ptr : m_entities)
        {
            if (entity_ptr->getName() == name)
                return entity_ptr.get();
        }
        return nullptr;
    }

    MeshEntity *Mesh::getEntity(CfemInt id)
    {
        for (auto &entity_ptr : m_entities)
        {
            if (entity_ptr->getId() == id)
                return entity_ptr.get();
        }
        return nullptr;
    }

    void Mesh::setOrder(MeshOrder newOrder)
    {
        CFEM_ASSERT(newOrder == MeshOrder::Linear || newOrder == MeshOrder::Quadratic);
        if (newOrder == m_meshOrder)
            return;

        if (!m_cells.empty())
        {
            if (newOrder == MeshOrder::Quadratic && m_meshOrder == MeshOrder::Linear)
            {
                this->upgradeToQuadratic();
            }
            else if (newOrder == MeshOrder::Linear && m_meshOrder == MeshOrder::Quadratic)
            {
                this->downgradeToLinear();
            }
        }
        m_meshOrder = newOrder;
    }

    void Mesh::setDimension(int newDim)
    {
        CFEM_ASSERT(newDim >= 1 || newDim <= 3);
        if (!(this->m_cells.empty()))
        {
            CFEM_ERROR("ARE YOU TRYING TO SET THE DIMENSION OF A MESH THAT ALREADY CONTAINS ELEMENTS? MAKE SURE THE MESH IS EMPTY.");
        }
        m_meshDim = newDim;
    }

    void Mesh::setReferenceNormal(const cfem::Vector3D& norm) {
        m_global_reference_normal = norm;
    }

    void Mesh::clearReferenceNormal() {
        m_global_reference_normal = std::nullopt;
    }

    MeshEntity &Mesh::createInternalEntity(CfemInt id, const std::string &name, CfemInt dim)
    {
        CFEM_PRECONDITION(m_entity_ids.find(id) == m_entity_ids.end());

        m_entity_ids.insert(id);
        m_entities.push_back(std::make_unique<MeshEntity>(this, id, name, dim));

        return *m_entities.back();
    }


    void Mesh::buildMesh()
    {
        if (m_cells.empty() || m_vertices.size() == 0 || m_meshDim == 0) {
            CFEM_ERROR("buildMesh failed: Mesh is empty.");
            return;
        }

        // 1. Build topological structures (The Golden Sequence)
        compactVertices();
        buildInverseConnectivities();
        buildFacets(); 
        buildAdjacency();

        // 2. Initialize Core Entities
        MeshEntity* universe = getEntity(cfem::UNIVERSE_ID);
        if (!universe) universe = &createInternalEntity(cfem::UNIVERSE_ID, "Universe", m_meshDim);

        MeshEntity* boundary = getEntity(cfem::SYS_BOUNDARY_ID);
        if (!boundary) boundary = &createInternalEntity(cfem::SYS_BOUNDARY_ID, "SYS_BOUNDARY", m_meshDim - 1);

        // 3. Populate Universe
        if (universe->isEmpty()) {
            universe->reserveVertices(m_vertices.size());
            universe->reserveCells(m_cells.size());
            for (CfemInt i = 0; i < (CfemInt)m_vertices.size(); ++i) universe->addVertex(i);
            for (CfemInt i = 0; i < (CfemInt)m_cells.size(); ++i)    universe->addCell(i);
        }

        // 4. Populate Facets and Identify Boundary
        // Use a bitset/bool vector to avoid massive duplicate insertions into 'boundary'
        std::vector<bool> vertex_on_boundary(m_vertices.size(), false);

        for (const auto &facet : m_facets) {
            universe->addFacet(facet.id);

            if (facet.isBoundary()) {
                boundary->addFacet(facet.id);
                auto f_vids = getFacetConnectivity(facet.id); 
                for (CfemInt v_id : f_vids) {
                    vertex_on_boundary[v_id] = true;
                }
            }
        }

        // Single-pass addition of boundary nodes
        for (CfemInt v_id = 0; v_id < (CfemInt)m_vertices.size(); ++v_id) {
            if (vertex_on_boundary[v_id]) boundary->addVertex(v_id);
        }

        // 5. Finalize
        universe->finalize();
        boundary->finalize();

        // 6. Validation and Reporting
    #ifndef NDEBUG
        ConformityReport report = MeshConformityChecker::check(*this, true);
        if (!report.is_valid) {
            CFEM_INFO(report.toString().c_str());
            CFEM_ERROR("Non conforming or degenerated mesh detected in debug build.");
        }
    #endif
        CFEM_INFO(("Mesh built: " + std::to_string(m_vertices.size()) + " Nodes, " 
                + std::to_string(m_cells.size()) + " Cells, " 
                + std::to_string(m_facets.size()) + " Facets.\n").c_str());

        if (!m_invertedCellCounts.empty())
        {
            size_t totalFixed = 0;
            std::string details = "";

            for (auto const& [type, count] : m_invertedCellCounts)
            {
                totalFixed += count;
                details += "\n  - " + cellTypeToString(type) + ": " + std::to_string(count);
            }

            CFEM_INFO("Orientation Check: " + std::to_string(totalFixed) + 
                    " inverted cells found and successfully fixed." + details);
            
            m_invertedCellCounts.clear();
        }
    }



    std::string Mesh::info() const
    {
        std::stringstream ss;
        ss << "\n================ MESH SUMMARY ================" << std::endl;
        ss << " Dimension: " << getDimension() << "D" << std::endl;
        ss << " Order:     " << (getOrder() == MeshOrder::Linear ? "Linear" : (getOrder() == MeshOrder::Quadratic) ? "Quadratic"
                                                                                                                   : "None")
           << std::endl;
        ss << " Vertices:  " << getNumVertices() << std::endl;
        ss << " Cells:     " << getNumCells() << std::endl;
        ss << " Facets:    " << getNumFacets() << std::endl;
        ss << "---------------------------------------------------" << std::endl;
        ss << std::left << std::setw(20) << "Entity Name"
           << std::setw(10) << "ID"
           << std::setw(15) << "Elements"
           << "Vertices" << std::endl;

        for (const auto &entity : getEntities())
        {
            std::string name = entity->getName();
            size_t count = (entity->getDimension() == getDimension())
                               ? entity->getNumCells()
                               : entity->getNumFacets();

            ss << std::left << std::setw(20) << name
               << std::setw(10) << entity->getId()
               << std::setw(15) << count
               << entity->getNumVertices() << std::endl;
        }
        ss << "===================================================\n"
           << std::endl;

        return ss.str();
    }


    std::string Mesh::getMeshFormat() const
{
    std::stringstream ss;
    ss << info();

    // --- Section 2: Nodes ---
    ss << "$Nodes\n";
    ss << m_vertices.size() << " # total number of nodes\n";
    for (size_t i = 0; i < m_vertices.size(); ++i)
    {
        const auto &v = m_vertices[i];
        ss << i << "  "
           << std::fixed << std::setprecision(6)
           << std::setw(10) << v.coords[0] << " "
           << std::setw(10) << v.coords[1] << " "
           << std::setw(10) << v.coords[2] << "\n";
    }
    ss << "$EndNodes\n";

    // --- Section 3: Cells (Updated with Neighbors) ---
    ss << "$Cells\n";
    ss << m_cells.size() << " # total number of cells\n";
    for (size_t i = 0; i < m_cells.size(); ++i)
    {
        const auto &cell = m_cells[i];
        auto cell_vids = getCellConnectivity(cell.id); 

        ss << i << " "
           << std::setw(8) << cellTypeToString(cell.type) << " "
           << cell_vids.size() << "  ";

        // Print Node IDs
        for (auto v_id : cell_vids) {
            ss << v_id << " ";
        }

        // --- NEW: Print Neighbors ---
        ss << "| Neighbors: ";
        auto neighbors = getCellNeighbors(cell.id);
        for (auto n_id : neighbors) {
            if (n_id == -1) ss << "B "; // B for Boundary
            else ss << n_id << " ";
        }
        ss << "\n";
    }
    ss << "$EndCells\n";

    // --- Section 4: Facets ---
    ss << "$Facets\n";
    ss << m_facets.size() << " # total number of facets\n";
    for (const auto &f : m_facets)
    {
        auto facet_vids = getFacetConnectivity(f.id);
        ss << f.id << " " << facet_vids.size() << "  ";

        for (auto v_id : facet_vids) {
            ss << v_id << " ";
        }

        ss << "| L: " << std::setw(8) << f.left_cell
           << " R: " << std::setw(8) << f.right_cell << "\n";
    }
    ss << "$EndFacets\n";

    return ss.str();
}

    void Mesh::upgradeToQuadratic() {
#ifndef NDEBUG
        CFEM_INFO( "Upgrading the mesh to quadratic... ");
#endif
        if (this->isQuadratic()) return;
        PROFILE_SCOPE("Mesh::upgradeToQuadratic");

        // Setup high-speed edge tracking
        struct Edge {
            CfemInt v1, v2;
            bool operator==(const Edge& o) const { 
                return (v1 == o.v1 && v2 == o.v2) || (v1 == o.v2 && v2 == o.v1); 
            }
        };
        struct EdgeHash {
            size_t operator()(const Edge& e) const {
                CfemInt a = std::min(e.v1, e.v2), b = std::max(e.v1, e.v2);
                return std::hash<long long>{}((static_cast<long long>(a) << 32) | b);
            }
        };

        struct QuadFace {
            CfemInt v1, v2, v3, v4; // Fixed typo: v3, v3 -> v3, v4

            // Equality operator: true if both faces have the same 4 vertices
            bool operator==(const QuadFace& o) const {
                std::array<CfemInt, 4> a = {v1, v2, v3, v4};
                std::array<CfemInt, 4> b = {o.v1, o.v2, o.v3, o.v4};
                
                std::sort(a.begin(), a.end());
                std::sort(b.begin(), b.end());
                
                return a == b;
            }
        };

        struct QuadFaceHash {
            size_t operator()(const QuadFace& f) const {
                std::array<CfemInt, 4> nodes = {f.v1, f.v2, f.v3, f.v4};
                std::sort(nodes.begin(), nodes.end());

                size_t seed = 0;
                // Standard hash combination (based on boost::hash_combine)
                // This efficiently mixes the sorted node IDs into a single size_t
                for (CfemInt v : nodes) {
                    seed ^= std::hash<CfemInt>{}(v) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                }
                return seed;
            }
        };
        std::unordered_map<Edge, CfemInt, EdgeHash> edge_to_midpoint;
        std::unordered_map<QuadFace, CfemInt, QuadFaceHash> quadface_to_midpoint;

        // Prepare the NEW Big Block for connectivity
        // Estimate: Quadratic connectivity is ~2.5x larger than linear
        std::vector<CfemInt> newAllCellVertexIds;
        newAllCellVertexIds.reserve(m_allCellVertexIds.size() * 3);

        std::vector<CfemInt> quad_buffer; 
        quad_buffer.reserve(27);

        // Loop over cells and upgrade
        for (size_t i = 0; i < m_cells.size(); ++i) {
            Cell &cell = m_cells[i];
            auto linear_vids = getCellConnectivity(cell.id); // Current safe span

            switch (cell.type) {
                case CellType::LINE2:  fillLine3Nodes(linear_vids, edge_to_midpoint, quad_buffer); break;
                case CellType::TRI3:   fillTri6Nodes(linear_vids, edge_to_midpoint, quad_buffer); break;
                case CellType::QUAD4:  fillQuad9Nodes(linear_vids, edge_to_midpoint, quad_buffer); break;
                case CellType::TET4:   fillTet10Nodes(linear_vids, edge_to_midpoint, quad_buffer); break;
                case CellType::PRISM6: fillPrism18Nodes(linear_vids, edge_to_midpoint, quadface_to_midpoint, quad_buffer); break;
                case CellType::HEX8:   fillHex27Nodes(linear_vids, edge_to_midpoint, quadface_to_midpoint, quad_buffer); break;
                default: CFEM_ERROR("UNSUPPORTED CELL TYPE."); break; 
            }

            // Update the Cell view to point to the NEW Big Block
            cell.offset = newAllCellVertexIds.size();
            cell.node_count = quad_buffer.size();
            cell.type = Cell::getQuadraticEquivalent(cell.type);

            // Append to the block
            newAllCellVertexIds.insert(newAllCellVertexIds.end(), quad_buffer.begin(), quad_buffer.end());
        }

        m_typeFlags.reset();
        for (auto &type : m_uniqueTypes)
        {
            type = Cell::getQuadraticEquivalent(type);
            int typeIdx = static_cast<int>(type);
            if (!m_typeFlags.test(typeIdx))
            {
                m_typeFlags.set(typeIdx);
            }
        }

        // Atomic Switch: Old linear block is deleted here
        m_allCellVertexIds = std::move(newAllCellVertexIds);
        m_meshOrder = MeshOrder::Quadratic;

        // Rebuild everything dependent on connectivity
        buildInverseConnectivities();
        buildFacets(); 
        buildAdjacency();
        updateEntities();

        for (auto it = m_observers.begin(); it != m_observers.end();)
        {
            if (auto space = it->lock())
            { // On vérifie si l'espace existe encore
                space->updateAfterMeshChange();
                ++it;
            }
            else
            {
                it = m_observers.erase(it); // Nettoyage automatique
            }
        }

#ifndef NDEBUG
        CFEM_INFO( "Mesh successfully upgraded to quadratic." );
        CFEM_INFO( this->info() );
#endif
    }


    void Mesh::addObserver(std::shared_ptr<FESpace> space)
    {
        m_observers.push_back(space);

        // Si le maillage est déjà quadratique, on force l'espace
        // qui vient d'arriver à se mettre à jour immédiatement.
        if (this->isQuadratic())
        {
            space->updateAfterMeshChange();
        }
    }

    void Mesh::downgradeToLinear()
    {
        if (this->isLinear())
            return;

        PROFILE_SCOPE("Mesh::downgradeToLinear");

        // Prepare a new Big Block for the linear connectivity
        // Linear connectivity is much smaller, so this is memory-safe
        std::vector<CfemInt> newLinearBlock;
        newLinearBlock.reserve(m_cells.size() * 8); // Heuristic reserve

        // Loop through cells and "Clip" the connectivity
        for (auto &cell : m_cells)
        {
            // Get the current quadratic nodes
            auto quad_nodes = getCellConnectivity(cell.id);
            
            // Determine how many corner nodes to keep
            int linearCount = 0;
            CellType linearType = cell.type;

            switch (cell.type)
            {
                case CellType::LINE3:   linearCount = 2; linearType = CellType::LINE2; break;
                case CellType::TRI6:    linearCount = 3; linearType = CellType::TRI3;  break;
                case CellType::QUAD9:   linearCount = 4; linearType = CellType::QUAD4; break;
                case CellType::TET10:   linearCount = 4; linearType = CellType::TET4;  break;
                case CellType::PRISM18: linearCount = 6; linearType = CellType::PRISM6; break;
                case CellType::HEX27:   linearCount = 8; linearType = CellType::HEX8;  break;
                default: 
                    linearCount = cell.node_count; // Already linear or unknown
                    break;
            }

            // Update the Cell Metadata
            cell.offset = newLinearBlock.size();
            cell.node_count = linearCount;
            cell.type = linearType;

            // Copy only the corners into the new block
            for (int i = 0; i < linearCount; ++i) {
                newLinearBlock.push_back(quad_nodes[i]);
            }
        }

        // Swap the blocks
        m_allCellVertexIds = std::move(newLinearBlock);
        m_meshOrder = MeshOrder::Linear;

        // Rebuild Topology
        // This will automatically "ignore" the mid-nodes still sitting in m_vertices
        buildInverseConnectivities();
        buildFacets(); 
        updateEntities();

        // Cleanup Observers
        for (auto it = m_observers.begin(); it != m_observers.end();) {
            if (auto space = it->lock()) {
                space->updateAfterMeshChange();
                ++it;
            } else {
                it = m_observers.erase(it);
            }
        }
    }

    /**
     * @brief
     *
     */
    void Mesh::updateEntities()
    {
        for (auto &entity_ptr : m_entities)
        {
            if (entity_ptr->getName() == "Universe")
            {
                entity_ptr->updateUniverse();
            }
            else if (entity_ptr->getName() == "Boundary")
            {
                entity_ptr->updateBoundary();
            }
            else
            {
                entity_ptr->update(); // Standard cell/facet-based update
            }
        }
    }

    //------------------------------------------------------------------------------------
    //------------------------- HELPER FUMCTIOMS -----------------------------------------
    //------------------------------------------------------------------------------------

    std::string Mesh::cellTypeToString(CellType type) const
    {
        switch (type)
        {
        case CellType::LINE2:
            return "LINE2";
        case CellType::TRI3:
            return "TRI3";
        case CellType::QUAD4:
            return "QUAD4";
        case CellType::TET4:
            return "TET4";
        case CellType::PRISM6:
            return "PRISM6";
        case CellType::HEX8:
            return "HEX8";
        case CellType::LINE3:
            return "LINE3";
        case CellType::TRI6:
            return "TRI6";
        case CellType::QUAD9:
            return "QUAD9";
        case CellType::TET10:
            return "TET10";
        case CellType::PRISM18:
            return "PRISM18";
        case CellType::HEX27:
            return "HEX27";
        // ... others
        default:
            return "UNKNOWN_CELL_TYPE";
        }
    }

    CfemReal Mesh::computeCellQuality(CfemInt cellId) const
    {
        const auto &cell = getCell(cellId);

        switch (cell.type)
        {
        case CellType::TRI3:
        case CellType::TRI6:
            return computeTriQuality(getCellConnectivity(cellId));

        case CellType::QUAD4:
        case CellType::QUAD9:
            return computeQuadQuality(getCellConnectivity(cellId));

        case CellType::TET4:
        case CellType::TET10:
            return computeTetQuality(getCellConnectivity(cellId));

        case CellType::PRISM6:
        case CellType::PRISM18:
            return computePrismQuality(getCellConnectivity(cellId));

        default:
            return 1.0; // Fallback for types not yet implemented
        }
    }

    CfemReal Mesh::computeQuadQuality(const std::span<const CfemInt>& v) const
    {
        if (v.size() < 4) return 0.0;

        CfemInt tri1[3] = { v[0], v[1], v[2] };
        CfemInt tri2[3] = { v[0], v[2], v[3] };

        CfemReal q1 = computeTriQuality(tri1);
        CfemReal q2 = computeTriQuality(tri2);
        
        return (q1 + q2) / 2.0;
    }

    CfemReal Mesh::computePrismQuality(const std::span<const CfemInt>& v) const
    {
        if (v.size() < 4) return 0.0;

        CfemInt tet1[4] = { v[0], v[1], v[2], v[3] };
        CfemInt tet2[4] = { v[1], v[4], v[2], v[3] };
        CfemInt tet3[4] = { v[4], v[5], v[2], v[3] };

        CfemReal q1 = computeTetQuality(tet1);
        CfemReal q2 = computeTetQuality(tet2);
        CfemReal q3 = computeTetQuality(tet3);
        
        return std::min({q1, q2, q3});
    }

    CfemReal Mesh::computeTriQuality(
        const std::span<const CfemInt>& v) const
    {
        const auto& p0 = getVertex(v[0]).coords;
        const auto& p1 = getVertex(v[1]).coords;
        const auto& p2 = getVertex(v[2]).coords;

        CfemReal a = (p1 - p2).norm();
        CfemReal b = (p0 - p2).norm();
        CfemReal c = (p0 - p1).norm();

        if (a <= ZERO_TOLERANCE ||
            b <= ZERO_TOLERANCE ||
            c <= ZERO_TOLERANCE)
            return 0.0;

        CfemReal quality = ((b + c - a) * (c + a - b) * (a + b - c)) / (a * b * c);

        return std::clamp(quality, CfemReal(0.0), CfemReal(1.0));
    }

    CfemReal Mesh::computeTetQuality(const std::span<const CfemInt> &v) const
    {
        const auto &p0 = getVertex(v[0]).coords;
        const auto &p1 = getVertex(v[1]).coords;
        const auto &p2 = getVertex(v[2]).coords;
        const auto &p3 = getVertex(v[3]).coords;

        // Edge vectors
        Vector3D e01 = p1 - p0;
        Vector3D e02 = p2 - p0;
        Vector3D e03 = p3 - p0;
        Vector3D e12 = p2 - p1;
        Vector3D e13 = p3 - p1;
        Vector3D e23 = p3 - p2;

        // Signed volume
        CfemReal vol = (e01).dot(e02.cross(e03)) / 6.0;

        if (vol <= 0.0)
            return 0.0;

        // Sum of squared edge lengths
        CfemReal sum_l2 = e01.normSquare() + e02.normSquare() + e03.normSquare() + e12.normSquare() + e13.normSquare() + e23.normSquare();

        CfemReal quality =
            (12.0 * std::pow(3.0 * vol, 2.0/3.0)) / sum_l2;

        return std::min(1.0, quality);
    }

    CfemInt Mesh::getDimFromType(CellType type) const
    {
        switch (type)
        {
        case CellType::LINE2:
        case CellType::LINE3:
            return 1;
        case CellType::TRI3:
        case CellType::TRI6:
        case CellType::QUAD4:
        case CellType::QUAD9:
            return 2;
        case CellType::TET4:
        case CellType::TET10:
        case CellType::PRISM6:
        case CellType::PRISM18:
        case CellType::HEX8:
        case CellType::HEX27:
            return 3;
        default:
            return 0;
        }
    }

    MeshOrder Mesh::getCellOrderFromType(CellType type) const
    {
        switch (type)
        {
        case CellType::LINE2:
        case CellType::TRI3:
        case CellType::QUAD4:
        case CellType::TET4:
        case CellType::PRISM6:
        case CellType::HEX8:
            return MeshOrder::Linear;
        case CellType::LINE3:
        case CellType::TRI6:
        case CellType::QUAD9:
        case CellType::TET10:
        case CellType::HEX27:
            return MeshOrder::Quadratic;
        default:
            return MeshOrder::None;
        }
    }

    void Mesh::validateEntityName(const std::string &name) const
    {
        if (name == "Universe" || name == "SYS_BOUNDARY")
        {
            CFEM_ERROR("Reserved Name Violation: The name '" + name +
                       "' is reserved for internal mesh topology. " +
                       "Please use a different name for the physical group.");
        }
    }


    /**
     * @brief
     *
     * @param order
     * @param subDiv
     * @param bounds
     * @param elemType
     * @return std::shared_ptr<Mesh>
     */
    std::shared_ptr<Mesh> Mesh::quad(
                                        MeshOrder order,
                                        const QuadSubdivisions& subDiv,
                                        const QuadBounds& bounds,
                                        QuadElementType elemType,
                                        const QuadBoundaryNames& names)
    {
        const double dx = std::abs(bounds.max.x - bounds.min.x);
        const double dy = std::abs(bounds.max.y - bounds.min.y);
        const double dz = std::abs(bounds.max.z - bounds.min.z);

        int collapsed = 0;
        if (dx < ZERO_TOLERANCE) collapsed++;
        if (dy < ZERO_TOLERANCE) collapsed++;
        if (dz < ZERO_TOLERANCE) collapsed++;

        if (collapsed != 1)
            CFEM_ERROR("QuadBounds must collapse exactly one dimension.");

        QuadCorners corners =
            (dz < ZERO_TOLERANCE)
                ? QuadCorners({
                    Vector3D{bounds.min.x, bounds.min.y, bounds.min.z},
                    Vector3D{bounds.max.x, bounds.min.y, bounds.min.z},
                    Vector3D{bounds.max.x, bounds.max.y, bounds.min.z},
                    Vector3D{bounds.min.x, bounds.max.y, bounds.min.z}})
            : (dy < ZERO_TOLERANCE)
                ? QuadCorners({
                    Vector3D{bounds.min.x, bounds.min.y, bounds.min.z},
                    Vector3D{bounds.min.x, bounds.min.y, bounds.max.z},
                    Vector3D{bounds.max.x, bounds.min.y, bounds.max.z},
                    Vector3D{bounds.max.x, bounds.min.y, bounds.min.z}})
            : QuadCorners({
                Vector3D{bounds.min.x, bounds.min.y, bounds.min.z},
                Vector3D{bounds.min.x, bounds.max.y, bounds.min.z},
                Vector3D{bounds.min.x, bounds.max.y, bounds.max.z},
                Vector3D{bounds.min.x, bounds.min.y, bounds.max.z}});

        return QuadMeshBuilder::create(order, subDiv, corners, elemType, names);
    }

    std::shared_ptr<Mesh> Mesh::quad(MeshOrder order,
                                        const QuadSubdivisions& subDiv,
                                        const QuadBounds& bounds,
                                        QuadElementType elemType)
    {
        return quad(order, subDiv, bounds, elemType, QuadBoundaryNames{});
    }


    /**
     * @brief 
     * 
     * @param order 
     * @param subDiv 
     * @param corners 
     * @param elemType 
     * @return std::shared_ptr<Mesh> 
     */
    std::shared_ptr<Mesh> Mesh::quad(MeshOrder order,
                                        const QuadSubdivisions& subDiv,
                                        const QuadCorners& corners,
                                        QuadElementType elemType,
                                        const QuadBoundaryNames& names)
    {
        return QuadMeshBuilder::create(order, subDiv, corners, elemType, names);
    }

    std::shared_ptr<Mesh> Mesh::quad(MeshOrder order,
                                        const QuadSubdivisions& subDiv,
                                        const QuadCorners& corners,
                                        QuadElementType elemType)
    {
        return quad(order, subDiv, corners, elemType, QuadBoundaryNames{});
    }


    /**
     * @brief 
     * 
     * @param order 
     * @param mesh_size 
     * @param center 
     * @param radius 
     * @param elemtype 
     * @return std::shared_ptr<Mesh> 
     */
    std::shared_ptr<Mesh> Mesh::disk(
                                MeshOrder order,
                                CfemReal mesh_size,
                                const Vector3D &center,
                                CfemReal radius,
                                const cfem::Vector3D &normal,
                                DiskElementType elemtype)
    {
        return DiskMeshBuilder::create(order, mesh_size, center, radius, normal, elemtype );
    }

    std::shared_ptr<Mesh> Mesh::disk(
                                MeshOrder order,
                                CfemReal mesh_size,
                                const Vector3D &center,
                                CfemReal radius,
                                DiskElementType elemtype)
    {
        return DiskMeshBuilder::create(order, mesh_size, center, radius, {0,0,1}, elemtype );
    }

    std::shared_ptr<Mesh> Mesh::annulus(
                                    MeshOrder order,
                                    CfemReal h,
                                    const cfem::Vector3D &center,
                                    CfemReal radius_in,
                                    CfemReal radius_out,
                                    DiskElementType type,
                                    const AnnulusBoundaryNames &names
                                )
    {
        return AnnulusMeshBuilder::create( order, h, center, radius_in, radius_out, DiskElementType::Triangle, names);
    }

    std::shared_ptr<Mesh> Mesh::annulus(
                                    MeshOrder order,
                                    CfemReal h,
                                    const cfem::Vector3D &center,
                                    CfemReal radius_in,
                                    CfemReal radius_out,
                                    DiskElementType type
                                )
    {
        return AnnulusMeshBuilder::create( order, h, center, radius_in, radius_out, DiskElementType::Triangle, AnnulusBoundaryNames{});
    }

    /**
     * @brief 
     * 
     * @param order 
     * @param mesh_size 
     * @param center 
     * @param radius 
     * @param elemtype 
     * @return std::shared_ptr<Mesh> 
     */
    std::shared_ptr<Mesh> Mesh::sphere(
                                    MeshOrder order,
                                    CfemReal mesh_size,
                                    const Vector3D &center,
                                    CfemReal radius,
                                    SphereElementType elemtype)
    {
        return SphereMeshBuilder::create(order, mesh_size, center, radius, elemtype );
    }


    /**
     * @brief 
     * 
     * @param order 
     * @param subDiv 
     * @param bounds 
     * @param elemType 
     * @param names 
     * @return std::shared_ptr<Mesh> 
     */
    std::shared_ptr<Mesh> Mesh::block(
                                MeshOrder order,
                                const BlockSubdivisions& subDiv,
                                const BlockBounds& bounds,
                                BlockElementType elemType,
                                const BlockBoundaryNames& names)
    {
        const double dx = std::abs(bounds.max.x - bounds.min.x);
        const double dy = std::abs(bounds.max.y - bounds.min.y);
        const double dz = std::abs(bounds.max.z - bounds.min.z);

        // Ensure this is a proper 3D volume
        if (dx < ZERO_TOLERANCE ||
            dy < ZERO_TOLERANCE ||
            dz < ZERO_TOLERANCE)
        {
            CFEM_ERROR("BlockBounds must define a non-degenerate 3D volume.");
        }

        // Standard hexahedron node ordering (right-handed, FEM-safe)
        //
        //   7 -------- 6
        //  /|         /|
        // 4 -------- 5 |
        // | |        | |
        // | 3 -------- 2
        // |/         |/
        // 0 -------- 1
        //
        // Bottom face: 0,1,2,3  (z = min)
        // Top face:    4,5,6,7  (z = max)

        BlockCorners corners({
            Vector3D{bounds.min.x, bounds.min.y, bounds.min.z}, // 0
            Vector3D{bounds.max.x, bounds.min.y, bounds.min.z}, // 1
            Vector3D{bounds.max.x, bounds.max.y, bounds.min.z}, // 2
            Vector3D{bounds.min.x, bounds.max.y, bounds.min.z}, // 3
            Vector3D{bounds.min.x, bounds.min.y, bounds.max.z}, // 4
            Vector3D{bounds.max.x, bounds.min.y, bounds.max.z}, // 5
            Vector3D{bounds.max.x, bounds.max.y, bounds.max.z}, // 6
            Vector3D{bounds.min.x, bounds.max.y, bounds.max.z}  // 7
        });

        return BlockMeshBuilder::create(order, subDiv, corners, elemType, names);
    }

    std::shared_ptr<Mesh> Mesh::block(MeshOrder order,
                                       const BlockSubdivisions& subDiv, 
                                       const BlockBounds& bounds,
                                       BlockElementType elemType)
    {
        return block(order, subDiv, bounds, elemType, BlockBoundaryNames{});
    }

    /**
     * @brief 
     * 
     * @param order 
     * @param subDiv 
     * @param corners 
     * @param elemType 
     * @param names 
     * @return std::shared_ptr<Mesh> 
     */
    std::shared_ptr<Mesh> Mesh::block(MeshOrder order,
                                       const BlockSubdivisions& subDiv, 
                                       const BlockCorners& corners,
                                       BlockElementType elemType,
                                       const BlockBoundaryNames& names)
    {
        return BlockMeshBuilder::create(order, subDiv, corners, elemType, names);

    }

    std::shared_ptr<Mesh> Mesh::block(MeshOrder order,
                                       const BlockSubdivisions& subDiv, 
                                       const BlockCorners& corners,
                                       BlockElementType elemType)
    {
        return BlockMeshBuilder::create(order, subDiv, corners, elemType, BlockBoundaryNames{});
    }

    std::shared_ptr<Mesh> Mesh::cylinder(
                                MeshOrder order,
                                CfemReal target_mesh_size,
                                const cfem::Vector3D &center,
                                CfemReal radius,
                                CfemReal height,
                                const cfem::Vector3D &axis,
                                CylinderElementType type,
                                const CylinderBoundaryNames &names)
    {
        return CylinderMeshBuilder::create(order, target_mesh_size, center, radius, height, axis, type, names);

    }

    std::shared_ptr<Mesh> Mesh::cylinder(MeshOrder order,
                                        CfemReal target_mesh_size,
                                        const cfem::Vector3D &center,
                                        CfemReal radius,
                                        CfemReal height,
                                        const cfem::Vector3D &axis,
                                        CylinderElementType type)
    {
        return CylinderMeshBuilder::create(order, target_mesh_size, center, radius, height, axis, type, CylinderBoundaryNames{});
    }


    std::shared_ptr<Mesh> Mesh::hollow_cylinder(
                                                    MeshOrder order,
                                                    CfemReal target_mesh_size,
                                                    const cfem::Vector3D &center,
                                                    CfemReal r_in,
                                                    CfemReal r_out,
                                                    CfemReal height,
                                                    const cfem::Vector3D &axis,
                                                    CylinderElementType type,
                                                    const HollowCylinderBoundaryNames &names
                                                )
    {
        return HollowCylinderMeshBuilder::create(order, target_mesh_size, center, r_in, r_out, height, axis, type, names);
    }
                                            
    std::shared_ptr<Mesh> Mesh::hollow_cylinder(
                                                    MeshOrder order,
                                                    CfemReal target_mesh_size,
                                                    const cfem::Vector3D &center,
                                                    CfemReal r_in,
                                                    CfemReal r_out,
                                                    CfemReal height,
                                                    const cfem::Vector3D &axis,
                                                    CylinderElementType type
                                                )
    {
        return HollowCylinderMeshBuilder::create(order, target_mesh_size, center, r_in, r_out, height, axis, type, HollowCylinderBoundaryNames{});
    }



    void Mesh::compactVertices()
    {
        const CfemInt oldVertexCount = static_cast<CfemInt>(m_vertices.size());
        if (oldVertexCount == 0) return;

        // ------------------------------------------------------------
        // 1) Mark used vertices
        // ------------------------------------------------------------
        std::vector<bool> used(oldVertexCount, false);

        for (CfemInt cellId = 0; cellId < static_cast<CfemInt>(m_cells.size()); ++cellId)
        {
            const auto conn = getCellConnectivity(cellId);
            for (CfemInt vid : conn)
            {
                if (vid >= 0 && vid < oldVertexCount)
                    used[vid] = true;
            }
        }

        // Count used
        CfemInt usedCount = 0;
        for (bool u : used)
            if (u) ++usedCount;

        // Nothing to do
        if (usedCount == oldVertexCount)
            return;

        // ------------------------------------------------------------
        // 2) Build old -> new index map
        // ------------------------------------------------------------
        std::vector<CfemInt> oldToNew(oldVertexCount, -1);

        std::vector<Vertex> newVertices;
        newVertices.reserve(usedCount);

        CfemInt newIndex = 0;
        for (CfemInt oldId = 0; oldId < oldVertexCount; ++oldId)
        {
            if (used[oldId])
            {
                oldToNew[oldId] = newIndex;
                newVertices.push_back(m_vertices[oldId]);
                ++newIndex;
            }
        }

        // ------------------------------------------------------------
        // 3) Remap cell connectivity (flat storage)
        // ------------------------------------------------------------
        for (CfemInt& vid : m_allCellVertexIds)
        {
            if (vid >= 0)
            {
                vid = oldToNew[vid];
            }
        }

        // ------------------------------------------------------------
        // 4) Remap facet connectivity (if stored explicitly)
        // ------------------------------------------------------------
        for (CfemInt& vid : m_allFacetVertexIds)
        {
            if (vid >= 0)
            {
                vid = oldToNew[vid];
            }
        }

        // ------------------------------------------------------------
        // 5) Remap entity vertex sets
        // ------------------------------------------------------------
        for (const auto& entity : m_entities)
        {
            std::vector<CfemInt> newEntityVertices;
            newEntityVertices.reserve(entity->getVertexIds().size());

            for (CfemInt oldVid : entity->getVertexIds())
            {
                if (oldVid >= 0 && oldVid < oldVertexCount)
                {
                    CfemInt mapped = oldToNew[oldVid];
                    if (mapped >= 0)
                        newEntityVertices.push_back(mapped);
                }
            }

            entity->replaceVertexList(newEntityVertices);
        }

        // ------------------------------------------------------------
        // 6) Replace vertex array
        // ------------------------------------------------------------
        m_vertices = std::move(newVertices);

        // ------------------------------------------------------------
        // 7) Clear adjacency (will be rebuilt)
        // ------------------------------------------------------------
        // buildInverseConnectivities();
        // buildFacets();
        // buildAdjacency();
    }

} // namespace