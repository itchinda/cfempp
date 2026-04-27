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

#include "Mesh.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <iomanip>
#include <sstream>

#include "QuadMeshBuilder.h"
#include "LineMeshBuilder.h"
#include "DiskMeshBuilder.h"
#include "BlockMeshBuilder.h"
#include "SolidSphereBuilder.h"
#include "SurfaceSphereBuilder.h"
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

    /** ************************************************************************************************************************************
     * @brief Add a new vertex to the mesh.
     *
     * Creates a vertex with the given coordinates and appends it to the
     * internal vertex container. The returned index uniquely identifies
     * the vertex in the mesh.
     *
     * Vertices should typically be added before any cells are defined.
     *
     * @param x X-coordinate of the vertex.
     * @param y Y-coordinate of the vertex.
     * @param z Z-coordinate of the vertex.
     *
     * @return The index of the newly created vertex.
     *
     * @note In 2D meshes the z-coordinate is ignored but still stored.
     ***************************************************************************************************************************************/
    CfemInt Mesh::addVertex(CfemReal x, CfemReal y, CfemReal z)
    {
        CfemInt id = m_vertices.size();
        m_vertices.push_back({{x, y, z}, {}});

        CFEM_POSTCONDITION(static_cast<CfemInt>(m_vertices.size()) > id);
        return id;
    }

    /** ************************************************************************************************************************************
     * @brief Add a new cell (finite element) to the mesh.
     *
     * Creates a cell of the specified type using the given vertex indices.
     * The connectivity is stored internally in a compact contiguous buffer.
     *
     * @param type The element type (triangle, quadrilateral, tetrahedron, etc.).
     * @param v_ids Vertex indices defining the element connectivity.
     *
     * @return The index of the newly created cell.
     *
     * @throws CFEM_ERROR If the connectivity is inconsistent
     *         with the specified element type.
     *
     * @note The orientation of the element may be corrected internally
     *       to ensure consistent geometry.
     ***************************************************************************************************************************************/
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
                CFEM_ERROR << "Orientation fix failed. Element might be mangled.";
            }
            m_invertedCellCounts[type]++;
        }
        else if (std::abs(measure) <= cfem::ZERO_TOLERANCE)
        {
            CFEM_WARN << "Degenerate (zero area/volume) cell added to mesh.";
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
        if (!m_cellTypeFlags.test(typeIdx))
        {
            m_cellTypeFlags.set(typeIdx);
            m_uniqueCellTypes.push_back(type);
        }


        auto facetTypes = Cell::getFacetTypes(type); 
        
        for (CellType fType : facetTypes)
        {
            int fTypeIdx = static_cast<int>(fType);
            
            // Si ce type de facette n'a pas encore été rencontré dans tout le maillage
            if (!m_facetTypeFlags.test(fTypeIdx))
            {
                m_facetTypeFlags.set(fTypeIdx);
                m_uniqueFacetTypes.push_back(fType);
            }
        }

        return id;
    }

    /** ************************************************************************************************************************************
     * @brief Validate the consistency of a cell before insertion into the mesh.
     *
     * This routine verifies that the provided element definition is compatible
     * with the current mesh configuration. It is typically called internally
     * by `addCell()` before the cell is inserted into the mesh data structures.
     *
     * The following checks are performed:
     *
     * - **Polynomial order compatibility**  
     *   Ensures that the element type matches the mesh order
     *   (e.g., linear elements cannot be inserted into a quadratic mesh).
     *
     * - **Topological dimension compatibility**  
     *   Ensures that the element dimension matches the mesh dimension.
     *   For example, a 3D element cannot be added to a 2D mesh.
     *
     * - **Node count validation**  
     *   Ensures that the number of provided vertex indices matches the
     *   expected number of nodes for the given cell type.
     *
     * - **Vertex index bounds checking**  
     *   Ensures that all vertex indices refer to existing mesh vertices.
     *
     * If any inconsistency is detected, the function aborts with a
     * `CFEM_ERROR`.
     *
     * @param type Type of the cell being added.
     * @param v_ids List of vertex indices defining the cell connectivity.
     *
     * @throws Runtime error if:
     * - the cell order differs from the mesh order,
     * - the element dimension differs from the mesh dimension,
     * - the number of provided nodes is incorrect,
     * - a vertex index is invalid.
     *
     * @note
     * This function is part of the mesh integrity safeguards and should
     * normally only be used internally by the mesh construction routines.
     ***************************************************************************************************************************************/
    void Mesh::checkConsistency(CellType type, const std::vector<CfemInt> &v_ids) const
    {
        MeshOrder cellOrder = getCellOrderFromType(type);
        if (cellOrder != m_meshOrder)
        {
            CFEM_ERROR << "Cannot add " << cellTypeToString(type) <<
                       " to a " << (m_meshOrder == MeshOrder::Linear ? "linear mesh." 
                                 : m_meshOrder == MeshOrder::Quadratic ? "quadratic mesh."
                                                                       : "mesh with unkmown order.");
        }

        CfemInt cellDim = getDimFromType(type);
        if (cellDim != m_meshDim)
        {
            CFEM_ERROR << "Cannot add a " << std::to_string(cellDim) << "D element (" <<
                       cellTypeToString(type) << ") to a " << std::to_string(m_meshDim) << "D mesh.";
        }

        if (v_ids.size() != getNodeCountPerCell(type))
        {
            CFEM_ERROR << "Topology error: " << cellTypeToString(type) <<
                       " expects " << std::to_string(getNodeCountPerCell(type)) <<
                       " nodes, but " << std::to_string(v_ids.size()) << " were provided.";
        }

        for (CfemInt v_id : v_ids)
        {
            if (v_id < 0 || v_id >= static_cast<CfemInt>(m_vertices.size()))
            {
                CFEM_ERROR << "Index out of bounds: vertex " << std::to_string(v_id) <<
                           " does not exist in a mesh with " <<
                           std::to_string(m_vertices.size()) << " vertices.";
            }
        }
    }

    /** ************************************************************************************************************************************
     * @brief Correct the orientation of a cell by reordering its vertices.
     *
     * Some mesh generators or external mesh files may produce elements
     * with an orientation opposite to the convention used by CFEM.
     * Incorrect orientation can lead to:
     *
     * - negative Jacobians
     * - incorrect normal vectors
     * - invalid finite element mappings
     *
     * This function fixes the orientation by applying a permutation of the
     * node indices that flips the element orientation while preserving the
     * element topology.
     *
     * The reordering rules depend on the specific cell type and ensure that
     * both corner nodes and higher-order nodes (mid-edge or face nodes) are
     * swapped consistently.
     *
     * Supported element types include:
     *
     * - Triangles (`TRI3`, `TRI6`)
     * - Quadrilaterals (`QUAD4`, `QUAD9`)
     * - Tetrahedra (`TET4`, `TET10`)
     * - Prisms (`PRISM6`, `PRISM18`)
     * - Hexahedra (`HEX8`, `HEX27`)
     *
     * The vertex array is modified **in place**.
     *
     * @param type Type of the element whose orientation must be corrected.
     * @param vids Vertex indices defining the element connectivity.
     *
     * @throws Runtime error if the cell type is not supported.
     *
     * @note
     * This function does **not detect** whether an element is incorrectly
     * oriented. It simply applies the permutation that flips the orientation.
     * Detection of invalid orientation is handled elsewhere in the mesh
     * construction pipeline.
     *
     * @warning
     * The ordering rules must remain consistent with the element topology
     * defined in `cfem_topology.h`. Changing these permutations without
     * updating the topology definitions may corrupt the mesh.
     ***************************************************************************************************************************************/
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
                // Swap Corners: Bottom (0,1,2) with Top (3,4,5)
                for (int i = 0; i < 3; ++i) std::swap(vids[i], vids[i + 3]);

                // Swap Horizontal Edges: Bottom (6,7,8) with Top (9,10,11)
                for (int i = 6; i < 9; ++i) std::swap(vids[i], vids[i + 3]);
                // Vertical edges 12,13,14 and face centers 15,16,17 stay put.
                break;
            }
            default: CFEM_ERROR << "Unsupported cell type.";
        }
    }

    /** ************************************************************************************************************************************
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
     ***************************************************************************************************************************************/
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

    /** ************************************************************************************************************************************
     * @brief Construct cell-to-cell adjacency information.
     *
     * This routine builds the neighbor connectivity for all cells in the mesh.
     * For each cell face, the algorithm determines which neighboring cell shares the corresponding facet.
     *
     * The resulting adjacency information is stored in compressed contiguous arrays owned by the Mesh:
     *
     * - `m_allCellNeighbors` : neighbor cell indices
     * - `m_allCellFacetIds`  : facet identifiers associated with each face
     * - `m_neighborOffsets`  : offsets into the above arrays
     *
     * The layout follows a CSR-style structure:
     *
     * @code
     * cell i neighbors are stored in 
     * m_allCellNeighbors[m_neighborOffsets[i] ... m_neighborOffsets[i+1]-1]
     * @endcode
     *
     * Each entry corresponds to one local face of the cell. If the face lies
     * on the boundary, the neighbor index is set to `-1`.
     *
     * In addition, the method initializes pointers inside each `Cell` structure so that the cell can directly access its neighbors and
     * associated facets without additional lookups.
     *
     * Algorithm overview:
     *
     * 1. Determine the number of faces for each cell.
     * 2. Allocate contiguous storage for neighbor relations.
     * 3. Loop over all facets in the mesh.
     * 4. For each facet:
     *    - Identify the local face index in the left cell
     *    - Record the right cell as its neighbor
     *    - Repeat symmetrically for the right cell (if not boundary)
     * 5. Link the resulting arrays back to the `Cell` structures.
     *
     * @note
     * This routine assumes that `buildFacets()` has already been executed, since facet ownership (`left_cell`, `right_cell`) is required to
     * establish adjacency.
     *
     * @warning
     * This method overwrites any existing adjacency information.
     *
     * @complexity
     * Approximately linear in the number of facets.
     ***************************************************************************************************************************************/
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

            auto getLocalFaceIdx = [&](CfemInt cell_id) -> int {
                int num_faces = getLocalFacesCount(m_cells[cell_id].type);
                for (int i = 0; i < num_faces; ++i) {
                    LocalFace currentLocal = getLocalFaceData(cell_id, i);
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
                size_t pos = m_neighborOffsets[f.left_cell] + localIdxL;
                m_allCellNeighbors[pos] = f.right_cell;
                m_allCellFacetIds[pos] = f.id;
            }

            // Link Right Cell
            if (!f.isBoundary()) {
                int localIdxR = getLocalFaceIdx(f.right_cell);
                if (localIdxR != -1) {
                    size_t pos = m_neighborOffsets[f.right_cell] + localIdxR;
                    m_allCellNeighbors[pos] = f.left_cell;
                    m_allCellFacetIds[pos] = f.id;
                }
            }
        }

        // Link the pointers in the Cell structs (Consistency Pass)
        for (size_t i = 0; i < numCells; ++i) {
            Cell& cell = m_cells[i];
            size_t startIdx = m_neighborOffsets[i];
            cell.num_faces = static_cast<CfemInt>(getLocalFacesCount(cell.type));
            
            // Point to the data inside the Mesh's big blocks
            cell.neighbor_cells_ptr = &m_allCellNeighbors[startIdx];
            cell.facet_ids_ptr     = &m_allCellFacetIds[startIdx];
        }
    }

    /** ************************************************************************************************************************************
     * @brief Construct mesh facets (edges in 2D, faces in 3D).
     *
     * This routine detects and constructs all facets shared by mesh cells. A facet represents the interface between two neighboring cells, or
     * a boundary surface if only one cell is attached.
     *
     * The resulting facets are stored in the mesh as:
     *
     * - `m_facets`              : list of facet descriptors
     * - `m_allFacetVertexIds`   : contiguous storage of facet vertex indices
     * - `m_facetOffsets`        : CSR-style offsets into the connectivity block
     *
     * Each facet records:
     *
     * - the owning cells (`left_cell`, `right_cell`)
     * - the facet topology type
     * - the offset and node count into the global connectivity array
     *
     * Boundary facets are identified when no neighboring cell is found.
     *
     * ---
     *
     * ### Algorithm Overview
     *
     * The construction proceeds in three stages:
     *
     * **1. Initialization**
     *
     * Existing facet data structures are cleared and storage is reserved based on an estimate of the number of facets.
     *
     * Temporary metadata structures are also prepared to store facet information during construction.
     *
     * **2. Topology Construction**
     *
     * Each cell face is processed:
     *
     * - A local face representation is generated (`LocalFace`)
     * - Neighbor cells are searched by examining cells connected to one anchor vertex of the face
     * - If another cell contains a matching face, the two cells share the same facet
     *
     * The face connectivity is then stored in a global contiguous
     * connectivity block.
     *
     * **3. Finalization**
     *
     * Temporary metadata is converted into permanent `Facet` objects
     * stored in `m_facets`.
     *
     * Each facet records:
     *
     * - its ID
     * - its topology type
     * - its left and right cells
     * - its connectivity offset
     * - its node count
     *
     * ---
     *
     * ### Data Layout
     *
     * Facet connectivity is stored using a contiguous block format:
     *
     * @code
     * m_allFacetVertexIds
     * [v0 v1 v2 v3 | v4 v5 v6 | v7 v8 ...]
     * @endcode
     *
     * Each facet stores:
     *
     * - `offset`     : start position in the block
     * - `node_count` : number of vertices
     *
     * This avoids per-facet allocations and improves cache locality.
     *
     * ---
     *
     * @note
     * This routine requires that inverse connectivities (`Vertex -> connected cells`) have already been built.
     *
     * @warning
     * Calling this function invalidates any previously computed facet
     * information.
     *
     * @complexity
     * Approximately linear in the number of cells, assuming bounded vertex connectivity.
     ***************************************************************************************************************************************/
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
            CfemInt left_facet_id;
            CfemInt right_facet_id;
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
                    (size_t)current_face.size,
                    f_idx,
                    neighbor_f_idx
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

            f.left_cell_facet_id = meta.left_facet_id;
            f.right_cell_facet_id = meta.right_facet_id;
            
            // No pointers! Just integers. 
            // This is now 100% safe from reallocations.
            f.offset = meta.offset;
            f.node_count = meta.count;
            
        }
    }

    /** ************************************************************************************************************************************
     * @name Mesh Entities
     *
     * Mesh entities represent labeled geometric subsets of the mesh such as
     * boundaries, regions, or internal interfaces. Entities are commonly used
     * to assign boundary conditions, material regions, or source terms.
     *
     * Entities are identified by both:
     * - a unique integer identifier
     * - a unique string name
     *
     * Internally, entities are stored as `MeshEntity` objects owned by the mesh.
     *
     * @{
     ***************************************************************************************************************************************/

     /** ************************************************************************************************************************************
     * @brief Create a new user-defined mesh entity.
     *
     * This method registers a new mesh entity with the given name andtopological dimension. 
     * Entities are typically used to represent physical regions or boundaries on which conditions will be applied.
     *
     * The mesh automatically assigns a unique internal identifier to the entity. The identifier is guaranteed not
     * to collide with existing entities.
     *
     * Some identifiers may be reserved for internal system entities (e.g. default boundary markers), and the numbering is adjusted
     * accordingly to avoid conflicts.
     *
     * @param name Unique name of the entity.
     * @param dim Topological dimension of the entity:
     *            - 0 : vertices
     *            - 1 : edges
     *            - 2 : faces
     *            - 3 : volumetric regions
     *
     * @return The unique identifier assigned to the new entity.
     *
     * @throws Runtime error if the entity name is invalid or already used.
     *
     * @note
     * Entity names must satisfy the mesh naming policy enforced by `validateEntityName()`.
     ***************************************************************************************************************************************/
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

    /** ************************************************************************************************************************************
     * @brief Retrieve an entity by its name.
     *
     * Searches the mesh entity registry and returns the entity whose
     * name matches the provided string.
     *
     * @param name Name of the entity.
     *
     * @return Pointer to the corresponding entity, or `nullptr`
     *         if no entity with this name exists.
     *
     * @note
     * The returned pointer is owned by the mesh and must not be deleted.
     ***************************************************************************************************************************************/
    MeshEntity *Mesh::getEntity(const std::string &name)
    {
        for (auto &entity_ptr : m_entities)
        {
            if (entity_ptr->getName() == name)
                return entity_ptr.get();
        }
        return nullptr;
    }

    const MeshEntity& Mesh::getEntityConst(const std::string &name) const
    {
        for (auto &entity_ptr : m_entities)
        {
            if (entity_ptr->getName() == name)
                return *entity_ptr;
        }
        CFEM_ERROR << "The mesh does not hold any entity with name " << name;
        return *m_entities[0]; // Makes the compiler happy
    }

    /** ************************************************************************************************************************************
     * @brief Retrieve an entity by its identifier.
     *
     * Searches the mesh entity registry and returns the entity with the specified internal identifier.
     *
     * @param id Identifier of the entity.
     *
     * @return Pointer to the entity if found, otherwise `nullptr`.
     *
     * @note
     * The returned pointer remains owned by the mesh.
     ***************************************************************************************************************************************/
    MeshEntity *Mesh::getEntity(CfemInt id)
    {
        for (auto &entity_ptr : m_entities)
        {
            if (entity_ptr->getId() == id)
                return entity_ptr.get();
        }
        return nullptr;
    }

    const MeshEntity& Mesh::getEntityConst(CfemInt id) const
    {
        for (auto &entity_ptr : m_entities)
        {
            if (entity_ptr->getId() == id)
                return *entity_ptr;
        }
        CFEM_ERROR << "The mesh does not hold any entity with id " << id;
        return *m_entities[0]; // Makes the compiler happy
    }

    bool Mesh::isPartOfTheMeshEntities(const MeshEntity& entity) const{
        return isPartOfTheMeshEntities(entity.getId());
    }

    bool Mesh::isPartOfTheMeshEntities(const std::string& entity_name) const{
        for (auto &entity_ptr : m_entities)
        {
            if (entity_ptr->getName() == entity_name)
                return true;
        }
        return false;
    }

    bool Mesh::isPartOfTheMeshEntities(CfemInt entity_id) const{
        for (auto &entity_ptr : m_entities)
        {
            if (entity_ptr->getId() == entity_id)
                return true;
        }

        return false;
    }

    /** ************************************************************************************************************************************
     * @brief Create a system-level mesh entity with a predefined identifier.
     *
     * This function is used internally by the mesh infrastructure to create
     * special entities such as system-defined boundaries or regions.
     *
     * Unlike `createEntity()`, the identifier is explicitly specified and
     * must not conflict with existing entities.
     *
     * @param id Unique identifier assigned to the entity.
     * @param name Name of the entity.
     * @param dim Topological dimension of the entity.
     *
     * @return Reference to the newly created entity.
     *
     * @pre The provided identifier must not already exist.
     *
     * @warning
     * This function bypasses the automatic ID generation used by
     * `createEntity()` and should normally only be used by internal
     * mesh-building routines.
     ***************************************************************************************************************************************/
    MeshEntity &Mesh::createInternalEntity(CfemInt id, const std::string &name, CfemInt dim)
    {
        CFEM_PRECONDITION(m_entity_ids.find(id) == m_entity_ids.end());

        m_entity_ids.insert(id);
        m_entities.push_back(std::make_unique<MeshEntity>(this, id, name, dim));

        return *m_entities.back();
    }
    /// @}

    /** ************************************************************************************************************************************
     * @name Mesh Configuration
     *
     * Functions used to configure fundamental mesh properties such as
     * polynomial order, spatial dimension, and geometric conventions.
     *
     * @{
     ***************************************************************************************************************************************/

    /** ************************************************************************************************************************************
     * @brief Change the polynomial order of the mesh.
     *
     * This method switches the mesh between linear and quadratic representations.
     *
     * If the mesh already contains elements, the connectivity will be transformed accordingly:
     *
     * - **Linear → Quadratic** : mid-edge (and possibly face) nodes are inserted using `upgradeToQuadratic()`.
     *
     * - **Quadratic → Linear** : higher-order nodes are removed using  `downgradeToLinear()`.
     *
     * If the mesh is empty, the order is simply updated.
     *
     * @param newOrder Desired mesh order (`Linear` or `Quadratic`).
     *
     * @note
     * Changing the mesh order modifies element connectivity and may invalidate existing finite element spaces or fields defined on
     * the mesh.
     ***************************************************************************************************************************************/
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

    /** ************************************************************************************************************************************
     * @brief Set the spatial dimension of the mesh.
     *
     * This function defines whether the mesh represents a 1D, 2D, or 3D domain.
     *
     * The dimension must be specified **before elements are inserted** into the mesh.
     *
     * @param newDim Mesh dimension (1, 2, or 3).
     *
     * @throws Runtime error if the mesh already contains elements.
     *
     * @note
     * Changing the dimension after elements have been added would invalidate the topology and is therefore prohibited.
     ***************************************************************************************************************************************/
    void Mesh::setDimension(int newDim)
    {
        CFEM_ASSERT(newDim >= 1 || newDim <= 3);
        if (!(this->m_cells.empty()))
        {
            CFEM_ERROR << "ARE YOU TRYING TO SET THE DIMENSION OF A MESH THAT ALREADY CONTAINS ELEMENTS? MAKE SURE THE MESH IS EMPTY.";
        }
        m_meshDim = newDim;
    }

    /** ************************************************************************************************************************************
     * @brief Define a global reference normal vector.
     *
     * Some mesh operations require a consistent orientation convention,
     * particularly when importing external meshes or constructing surface
     * meshes embedded in 3D.
     *
     * This method defines a global reference normal that can be used to
     * orient elements consistently.
     *
     * @param norm Reference normal vector.
     ***************************************************************************************************************************************/
    void Mesh::setReferenceNormal(const cfem::Vector3D& norm) {
        m_global_reference_normal = norm;
    }

    /** ************************************************************************************************************************************
     * @brief Remove the global reference normal.
     *
     * After calling this function, the mesh will no longer enforce a specific global orientation constraint.
     ***************************************************************************************************************************************/
    void Mesh::clearReferenceNormal() {
        m_global_reference_normal = std::nullopt;
    }

    /** ************************************************************************************************************************************
     * @brief Register a finite element space as an observer of the mesh.
     *
     * This method adds a finite element space to the list of observers that must be notified when the mesh topology changes (for example after
     * upgrading the mesh order or modifying connectivity).
     *
     * Observers are stored internally as weak references and are notified whenever the mesh structure is modified. This mechanism ensures that
     * finite element spaces remain consistent with the underlying mesh.
     *
     * If the mesh is already quadratic when the observer is added, the space is immediately notified so it can update its internal data structures.
     * This ensures that newly registered spaces are synchronized with the current mesh state.
     *
     * Typical updates triggered by observers include:
     * - rebuilding DOF numbering
     * - updating element connectivity
     * - resizing internal data structures
     *
     * @param space Shared pointer to the finite element space observing the mesh.
     *
     * @note
     * The mesh only stores weak references to observers internally, allowing finite element spaces to be destroyed safely without requiring explicit
     * deregistration.
     *
     * @see MeshObserver::updateAfterMeshChange()
     ***************************************************************************************************************************************/
    void Mesh::addObserver(std::shared_ptr<MeshObserver> observer)
    {
        m_observers.push_back(observer);

        // Si le maillage est déjà quadratique, on force l'espace
        // qui vient d'arriver à se mettre à jour immédiatement.
        if (this->isQuadratic())
        {
            observer->updateAfterMeshChange();
        }
    }

    void Mesh::notifyObservers() {
        for (auto it = m_observers.begin(); it != m_observers.end();)
        {
            if (auto obs = it->lock())
            { // On vérifie si l'observateur existe encore
                obs->updateAfterMeshChange();
                ++it;
            }
            else {
                it = m_observers.erase(it); // Nettoyage automatique
            }
        }
    }

    /** ************************************************************************************************************************************
     * @brief Finalize mesh construction and build all topological structures.
     *
     * This routine completes the mesh construction process by generating all derived topological relationships required 
     * for finite element operations. It must be called after vertices and cells have been inserted into the mesh.
     *
     * The method performs the following steps:
     *
     * ---
     *
     * ### 1. Topology Construction
     *
     * Core connectivity structures are generated in the following order (often referred to internally as the *golden sequence*):
     *
     * 1. `compactVertices()`  
     *    Removes unused vertices and compacts vertex storage.
     *
     * 2. `buildInverseConnectivities()`  
     *    Constructs vertex → cell connectivity.
     *
     * 3. `buildFacets()`  
     *    Detects all mesh facets (edges in 2D, faces in 3D) and identifies
     *    cell-to-facet relationships.
     *
     * 4. `buildAdjacency()`  
     *    Builds cell-to-cell adjacency through shared facets.
     *
     * ---
     *
     * ### 2. Core Entity Initialization
     *
     * Two system entities are automatically created if they do not exist:
     *
     * - **Universe**  
     *   Represents the full mesh domain and contains all vertices, cells, and facets.
     *
     * - **SYS_BOUNDARY**  
     *   Represents the geometric boundary of the mesh (facets with only one adjacent cell).
     *
     * ---
     *
     * ### 3. Populate Universe Entity
     *
     * All mesh vertices and cells are added to the Universe entity.
     *
     * ---
     *
     * ### 4. Boundary Detection
     *
     * Boundary facets are detected by checking facets that have only one adjacent cell. These facets are added to the `SYS_BOUNDARY` entity.
     *
     * All vertices belonging to boundary facets are also registered as boundary vertices.
     *
     * ---
     *
     * ### 5. Entity Finalization
     *
     * After population, system entities are finalized to optimize their internal storage and connectivity structures.
     *
     * ---
     *
     * ### 6. Debug Validation
     *
     * In debug builds (`!NDEBUG`), a mesh conformity check is performed using `MeshConformityChecker`. The verification detects:
     *
     * - non-conforming mesh topology
     * - degenerated elements
     * - invalid connectivity
     *
     * If issues are detected, the build is aborted.
     *
     * ---
     *
     * ### 7. Orientation Fix Reporting
     *
     * If any inverted elements were detected and corrected during mesh construction, a summary report is printed indicating how many
     * elements were reoriented for each cell type.
     *
     * ---
     *
     * @pre
     * - The mesh must contain vertices and cells.
     * - The mesh dimension must be defined.
     *
     * @throws Runtime error if the mesh is empty or invalid.
     *
     * @note
     * This function should typically be called once after mesh construction. Subsequent modifications to vertices or cells may require rebuilding
     * the mesh topology.
     *
     * @warning
     * Calling this method invalidates any previously computed connectivity structures and mesh entities.
     *
     * @complexity
     * Approximately linear in the number of mesh cells and facets.
     ***************************************************************************************************************************************/
    void Mesh::buildMesh()
    {
        if (m_cells.empty() || m_vertices.size() == 0 || m_meshDim == 0) {
            CFEM_ERROR << "buildMesh failed: Mesh is empty.";
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
        updateAffineFlags();

        // 6. Validation and Reporting
    #ifndef NDEBUG
        ConformityReport report = MeshConformityChecker::check(*this, true);
        if (!report.is_valid) {
            CFEM_INFO << report.toString().c_str();
            CFEM_ERROR << "Non conforming or degenerated mesh detected in debug build.";
        }
    #endif
        CFEM_INFO << "Mesh built: " 
                         << std::to_string(m_vertices.size()) 
                         << " Nodes, " 
                         << std::to_string(m_cells.size()) + " Cells, " 
                         << std::to_string(m_facets.size()) + " Facets.\n";

        if (!m_invertedCellCounts.empty())
        {
            size_t totalFixed = 0;
            std::string details = "";

            for (auto const& [type, count] : m_invertedCellCounts)
            {
                totalFixed += count;
                details += "\n  - " + cellTypeToString(type) + ": " + std::to_string(count);
            }

            CFEM_INFO << "Orientation Check: " << std::to_string(totalFixed) << 
                    " inverted cells found and successfully fixed." << details;
            
            m_invertedCellCounts.clear();
        }
    }
    
    /** ************************************************************************************************************************************
     * @brief Remove unused vertices from the mesh and remap connectivity.
     *
     * This function "compacts" the vertex array by eliminating vertices that are
     * not referenced by any cell or facet. It updates all relevant connectivity
     * arrays, including:
     *   - Cell-to-vertex connectivity (m_allCellVertexIds)
     *   - Facet-to-vertex connectivity (m_allFacetVertexIds)
     *   - Vertex sets stored in MeshEntities
     *
     * After this function:
     *   - The m_vertices array contains only vertices actually used by the mesh.
     *   - All vertex IDs in cells, facets, and entities are remapped to the
     *     new compacted indices.
     *
     * This improves memory locality, reduces storage requirements, and ensures
     * that vertex IDs are contiguous from 0 to n-1.
     *
     * Note:
     *   - Adjacency structures (inverse connectivity, facets, cell neighbors) 
     *     are not rebuilt automatically in this function. Call
     *     `buildInverseConnectivities()`, `buildFacets()`, and `buildAdjacency()`
     *     afterward if needed.
     *
     * Steps performed:
     *   1) Mark all vertices used by any cell.
     *   2) Build an old-to-new vertex index map and a new vertex array.
     *   3) Remap cell connectivity according to the new indices.
     *   4) Remap facet connectivity if stored explicitly.
     *   5) Remap vertex sets in all mesh entities.
     *   6) Replace m_vertices with the compacted vertex array.
     *   7) Optionally, clear/rebuild adjacency structures.
     ***************************************************************************************************************************************/
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


    /** ************************************************************************************************************************************
     * @brief Generate a human-readable summary of the mesh.
     *
     * This function returns a formatted string containing key information about the mesh structure, including:
     *
     * - spatial dimension
     * - polynomial order
     * - number of vertices
     * - number of cells
     * - number of facets
     *
     * It also lists all registered mesh entities and reports, for each:
     *
     * - entity name
     * - entity identifier
     * - number of associated elements (cells or facets depending on dimension)
     * - number of associated vertices
     *
     * The output is intended for debugging, logging, and quick inspection of the mesh state.
     *
     * Example output:
     *
     * @code
     * ================ MESH SUMMARY ================
     *  Dimension: 3D
     *  Order:     Linear
     *  Vertices:  125
     *  Cells:     320
     *  Facets:    480
     * ---------------------------------------------------
     *  Entity Name         ID        Elements       Vertices
     *  Universe            0         320            125
     *  SYS_BOUNDARY        1         96             64
     * ===================================================
     * @endcode
     *
     * @return A formatted string describing the mesh.
     *
     * @note
     * This function does not modify the mesh and is safe to call at any time after mesh construction.
     ***************************************************************************************************************************************/
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

    /** ************************************************************************************************************************************
     * @brief Export the mesh in a detailed textual debug format.
     *
     * This function generates a human-readable representation of the mesh including its topology and connectivity. The output is primarily
     * intended for debugging, inspection, and development purposes.
     *
     * The output begins with the mesh summary produced by `info()`, followed by structured sections describing nodes, cells, and facets.
     *
     * ---
     *
     * ### Output Structure
     *
     * 1. **Mesh Summary**
     *    The first section is the formatted summary returned by `info()`, containing global statistics and entity information.
     *
     * 2. **Nodes Section**
     *
     * ```
     * $Nodes
     * <number_of_nodes>
     * node_id  x  y  z
     * ...
     * $EndNodes
     * ```
     *
     * Lists all mesh vertices with their coordinates.
     *
     * 3. **Cells Section**
     *
     * ```
     * $Cells
     * <number_of_cells>
     * cell_id  type  node_count  node_ids | Neighbors: ...
     * ...
     * $EndCells
     * ```
     *
     * Each cell entry contains:
     *
     * - element type
     * - vertex connectivity
     * - neighboring cells
     *
     * Neighbor values are reported as:
     *
     * - a cell index for interior neighbors
     * - `B` if the face lies on the boundary.
     *
     * 4. **Facets Section**
     *
     * ```
     * $Facets
     * <number_of_facets>
     * facet_id node_count node_ids | L: <left_cell> R: <right_cell>
     * ...
     * $EndFacets
     * ```
     *
     * Each facet entry records:
     *
     * - its vertex connectivity
     * - the indices of adjacent cells
     *
     * Boundary facets have `R = -1`.
     *
     * ---
     *
     * @return A formatted string containing the mesh summary and full connectivity description.
     *
     * @note
     * This format is intended for debugging and topology inspection and should not be considered a stable mesh file format.
     *
     * @warning
     * The layout of this output may change between versions of the library.
     ***************************************************************************************************************************************/
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

    /** ************************************************************************************************************************************
     * @brief Upgrade the mesh from linear to quadratic order.
     *
     * This routine converts a linear mesh into a quadratic mesh by inserting additional nodes on element edges (and faces where required). The
     * connectivity of all cells is updated accordingly, and element types are replaced by their quadratic equivalents.
     *
     * The upgrade process creates:
     *
     * - **mid-edge nodes** for all elements
     * - **mid-face nodes** for elements that require them (e.g. HEX27, PRISM18)
     *
     * Newly created nodes are shared between neighboring elements to ensure a conforming quadratic mesh.
     *
     * ---
     *
     * ### Supported Element Transformations
     *
     * | Linear Type | Quadratic Type |
     * |-------------|----------------|
     * | LINE2       | LINE3          |
     * | TRI3        | TRI6           |
     * | QUAD4       | QUAD9          |
     * | TET4        | TET10          |
     * | PRISM6      | PRISM18        |
     * | HEX8        | HEX27          |
     *
     * ---
     *
     * ### Algorithm Overview
     *
     * 1. **Edge midpoint generation**
     *
     *    A hash map is used to track edges and ensure that each edge midpoint is created only once, even if the edge is shared by multiple cells.
     *
     * 2. **Face midpoint generation**
     *
     *    For elements that require face nodes (e.g. HEX27), a second hash map  tracks quadrilateral faces to ensure consistent face center creation.
     *
     * 3. **Connectivity reconstruction**
     *
     *    A new contiguous connectivity block is created containing the quadratic node ordering for all cells.
     *
     * 4. **Cell update**
     *
     *    Each cell updates:
     *
     *    - its node offset
     *    - its node count
     *    - its element type
     *
     * 5. **Topology rebuild**
     *
     *    Since connectivity changes, the following structures are rebuilt:
     *
     *    - vertex → cell connectivity
     *    - facets
     *    - cell adjacency
     *
     * 6. **Entity update**
     *
     *    Mesh entities are updated to reflect the new vertex set.
     *
     * 7. **Observer notification**
     *
     *    Registered finite element spaces observing the mesh are notified via `updateAfterMeshChange()` so that they can update their internal
     *    data structures.
     *
     * ---
     *
     * @note
     * If the mesh is already quadratic, this function returns immediately.
     *
     * @warning
     * This operation modifies the global mesh topology and invalidates previously stored connectivity information.
     *
     * @post
     * - All cells use quadratic element types.
     * - Additional vertices are inserted at edge and face midpoints.
     * - Mesh topology structures are rebuilt.
     *
     * @complexity
     * Approximately linear in the number of cells.
     ***************************************************************************************************************************************/
    void Mesh::upgradeToQuadratic() {
#ifndef NDEBUG
        CFEM_INFO << "Upgrading the mesh to quadratic... ";
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
        // struct EdgeHash {
        //     size_t operator()(const Edge& e) const {
        //         CfemInt a = std::min(e.v1, e.v2), b = std::max(e.v1, e.v2);
        //         return std::hash<long long>{}((static_cast<long long>(a) << 32) | b);
        //     }
        // };

        struct EdgeHash {
            size_t operator()(const Edge& e) const {
                size_t seed = 0;
                CfemInt a = std::min(e.v1, e.v2);
                CfemInt b = std::max(e.v1, e.v2);
                seed ^= std::hash<CfemInt>{}(a) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                seed ^= std::hash<CfemInt>{}(b) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                return seed;
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
        std::vector<CfemInt> newAllCellVertexIds;
        newAllCellVertexIds.reserve(m_allCellVertexIds.size() * 3);

        m_cellOffsets.clear();
        m_cellOffsets.reserve(m_cells.size() + 1);
        m_cellOffsets.push_back(0);

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
                default: CFEM_ERROR << "UNSUPPORTED CELL TYPE."; break; 
            }

            // Update the Cell view to point to the NEW Big Block
            cell.offset = newAllCellVertexIds.size();
            cell.node_count = quad_buffer.size();
            cell.type = Cell::getQuadraticEquivalent(cell.type);

            // Append to the block
            newAllCellVertexIds.insert(newAllCellVertexIds.end(), quad_buffer.begin(), quad_buffer.end());
            m_cellOffsets.push_back(newAllCellVertexIds.size());
        }

        m_cellTypeFlags.reset();
        for (auto &type : m_uniqueCellTypes)
        {
            type = Cell::getQuadraticEquivalent(type);
            int typeIdx = static_cast<int>(type);
            if (!m_cellTypeFlags.test(typeIdx))
            {
                m_cellTypeFlags.set(typeIdx);
            }
        }

        m_facetTypeFlags.reset();
        for (auto &type : m_uniqueFacetTypes)
        {
            type = Cell::getQuadraticEquivalent(type);
            int typeIdx = static_cast<int>(type);
            if (!m_facetTypeFlags.test(typeIdx))
            {
                m_facetTypeFlags.set(typeIdx);
            }
        }

        auto itCell = std::unique(m_uniqueCellTypes.begin(), m_uniqueCellTypes.end());
        m_uniqueCellTypes.erase(itCell, m_uniqueCellTypes.end());

        auto itFacet = std::unique(m_uniqueFacetTypes.begin(), m_uniqueFacetTypes.end());
        m_uniqueFacetTypes.erase(itFacet, m_uniqueFacetTypes.end());

        // Atomic Switch: Old linear block is deleted here
        m_allCellVertexIds = std::move(newAllCellVertexIds);
        m_meshOrder = MeshOrder::Quadratic;

        // Rebuild everything dependent on connectivity
        buildInverseConnectivities();
        buildFacets(); 
        buildAdjacency();
        updateEntities();
        // updateAffineFlags(); Pas necessaire rendu la! on a simplement ajouter des points milieux. si l'element etait affine, il le demeure.

        notifyObservers();

        CFEM_INFO <<  "Mesh successfully upgraded to quadratic.";
        CFEM_INFO <<  this->info();

    }


    /** ************************************************************************************************************************************
     * @brief Downgrade the mesh from quadratic to linear order.
     *
     * This routine converts a quadratic mesh into a linear mesh by removing higher-order nodes (mid-edge, mid-face, and interior nodes) from the
     * cell connectivity. Only the corner vertices of each element are kept, and the element types are replaced with their linear equivalents.
     *
     * The downgrade process does not delete the higher-order vertices from the global vertex list. Instead, the cell connectivity is rebuilt so
     * that only the corner nodes are referenced. As a result, unused vertices may remain in the mesh until a later compaction step is performed.
     *
     * ---
     *
     * ### Supported Element Transformations
     *
     * | Quadratic Type | Linear Type |
     * |----------------|------------|
     * | LINE3          | LINE2      |
     * | TRI6           | TRI3       |
     * | QUAD9          | QUAD4      |
     * | TET10          | TET4       |
     * | PRISM18        | PRISM6     |
     * | HEX27          | HEX8       |
     *
     * ---
     *
     * ### Algorithm Overview
     *
     * 1. **Connectivity clipping**
     *
     *    For each cell, only the first nodes corresponding to the corner vertices are copied into a new connectivity block.
     *
     * 2. **Cell update**
     *
     *    Each cell updates:
     *
     *    - its connectivity offset
     *    - its node count
     *    - its element type
     *
     * 3. **Connectivity replacement**
     *
     *    The new linear connectivity block replaces the quadratic one.
     *
     * 4. **Topology rebuild**
     *
     *    Since connectivity has changed, the following structures are rebuilt:
     *
     *    - vertex → cell inverse connectivity
     *    - facets
     *    - mesh entities
     *
     * 5. **Observer notification**
     *
     *    All registered finite element spaces observing the mesh are notified via `updateAfterMeshChange()` so they can rebuild their internal
     *    data structures.
     *
     * ---
     *
     * @note
     * If the mesh is already linear, the function returns immediately.
     *
     * @warning
     * This operation invalidates existing connectivity and topology data.
     *
     * @post
     * - All cells use linear element types.
     * - Only corner nodes remain in cell connectivity.
     * - Higher-order vertices may still exist in `m_vertices` but are unused.
     *
     * @complexity
     * Approximately linear in the number of cells.
     *
     * @see upgradeToQuadratic()
     ***************************************************************************************************************************************/
    void Mesh::downgradeToLinear()
    {
        if (this->isLinear())
            return;

        PROFILE_SCOPE("Mesh::downgradeToLinear");

        // 1. Préparation des nouveaux blocs de connectivité
        std::vector<CfemInt> newAllCellVertexIds;
        newAllCellVertexIds.reserve(m_cells.size() * 4); // Heuristique (moyenne 4 nœuds/cell)

        m_cellOffsets.clear();
        m_cellOffsets.reserve(m_cells.size() + 1);
        m_cellOffsets.push_back(0);

        // 2. Boucle sur les cellules pour "tronquer" la connectivité
        for (auto &cell : m_cells)
        {
            // Récupère les nœuds quadratiques actuels
            auto quad_nodes = getCellConnectivity(cell.id);
            
            CellType linearType = Cell::getLinearEquivalent(cell.type);
            int linearCount = Cell::getNodeCount(linearType);

            cell.offset = newAllCellVertexIds.size();
            cell.node_count = linearCount;
            cell.type = linearType;

            for (int i = 0; i < linearCount; ++i) {
                newAllCellVertexIds.push_back(quad_nodes[i]);
            }

            // Mise à jour de l'offset CSR pour la cellule suivante
            m_cellOffsets.push_back(newAllCellVertexIds.size());
        }

        // 3. Mise à jour des Flags de Types (Cells & Facets)
        m_cellTypeFlags.reset();
        for (auto &type : m_uniqueCellTypes)
        {
            type = Cell::getLinearEquivalent(type);
            m_cellTypeFlags.set(static_cast<int>(type));
        }
        // Nettoyage des doublons après conversion
        std::sort(m_uniqueCellTypes.begin(), m_uniqueCellTypes.end());
        m_uniqueCellTypes.erase(std::unique(m_uniqueCellTypes.begin(), m_uniqueCellTypes.end()), m_uniqueCellTypes.end());

        m_facetTypeFlags.reset();
        for (auto &type : m_uniqueFacetTypes)
        {
            type = Cell::getLinearEquivalent(type);
            m_facetTypeFlags.set(static_cast<int>(type));
        }
        std::sort(m_uniqueFacetTypes.begin(), m_uniqueFacetTypes.end());
        m_uniqueFacetTypes.erase(std::unique(m_uniqueFacetTypes.begin(), m_uniqueFacetTypes.end()), m_uniqueFacetTypes.end());

        // 4. Finalisation de la structure de données
        m_allCellVertexIds = std::move(newAllCellVertexIds);
        m_meshOrder = MeshOrder::Linear;

        compactVertices(); 
        buildInverseConnectivities();
        buildFacets(); 
        
        updateEntities();
        
        m_isGloballyAffine = true;
        m_isCellAffine.assign(m_cells.size(), true);

        notifyObservers();
    }


    /** ************************************************************************************************************************************
     * @brief Update all mesh entities after a topology modification.
     *
     * This method refreshes the contents of all registered mesh entities to ensure they remain consistent with the current mesh topology.
     *
     * Mesh entities represent subsets of the mesh such as:
     * - the global domain ("Universe")
     * - the global boundary
     * - user-defined physical groups
     *
     * Since operations like mesh refinement, order upgrade/downgrade, or connectivity rebuild may alter vertices, facets, or cells, entities
     * must be updated to reflect the new mesh state.
     *
     * The update behavior depends on the type of entity:
     *
     * - **Universe entity**: updated via `updateUniverse()`, which rebuilds the full set of cells, facets, and vertices belonging to the mesh.
     *
     * - **Boundary entity**: updated via `updateBoundary()`, which detects boundary facets and their associated vertices.
     *
     * - **User-defined entities**: updated using the generic `update()` method, which recomputes their vertex and element sets based on
     *   the entity definition. This method is typically invoked after operations that modify the mesh connectivity, such as:
     *
     * - `upgradeToQuadratic()`
     * - `downgradeToLinear()`
     * - mesh generation or rebuilding
     *
     * @note
     * The update preserves entity identifiers and names while rebuilding their internal lists of vertices, facets, and cells.
     ***************************************************************************************************************************************/
    void Mesh::updateEntities()
    {
        for (auto &entity_ptr : m_entities)
        {
            if (entity_ptr->getName() == "Universe")
            {
                entity_ptr->updateUniverse();
            }
            else if (entity_ptr->getName() == "SYS_BOUNDARY")
            {
                entity_ptr->updateBoundary();
            }
            else
            {
                entity_ptr->update(); // Standard cell/facet-based update
            }
        }
    }

    /**
     * @brief Updates the affinity status for each cell in the mesh.
     *
     * This function analyzes the geometric mapping of every cell to determine if it is 
     * "affine" (linear). A cell is considered affine if its Jacobian matrix is constant 
     * throughout its domain.
     *
     * @details
     * For linear elements (P1), the mapping is affine by definition. For quadratic 
     * elements (P2), the function performs the following checks:
     * * 1. **Edge Check**: Every edge must be straight. For each edge, the position 
     * of the quadratic midpoint is compared to the theoretical linear midpoint 
     * between the two edge corners.
     * * 2. **Face Check**: For quadrilateral faces (Quad9, Hex27, Prism18), the 
     * face-center node must be located at the geometric barycenter of the 
     * four face corners.
     * * 3. **Global Status**: Updates the @c m_isGloballyAffine flag, which allows 
     * the integration engine to trigger fast-path optimizations if the entire 
     * mesh consists of straight-edged elements.
     *
     * @note 
     * This function should be called after any mesh transformation, such as 
     * @c upgradeToQuadratic() or mesh smoothing, to ensure the integration 
     * kernels use the correct mapping logic.
     *
     * @param tol Numerical tolerance for comparing node positions (default: 1e-10).
     */
    void Mesh::updateAffineFlags() {
        const CfemInt numCells = getNumCells();
        
        // Initialiser tous les flags à true
        m_isCellAffine.assign(numCells, true);
        m_isGloballyAffine = true;

        if (m_meshOrder == MeshOrder::Linear) return;

        const CfemReal tol = cfem::ZERO_TOLERANCE * 1e4;

        #pragma omp parallel for reduction(&& : m_isGloballyAffine)
        for (CfemInt c = 0; c < numCells; ++c) {
            const auto& cell = getCell(c);
            
            if (cell.type == Cell::getLinearEquivalent(cell.type)) {
                continue; 
            }

            DynamicVector<Vector3D> coords;
            getCellVerticesCoords(c, coords);

            bool cellIsAffine = true;

            // =========================================================
            // VÉRIFICATION DES ARÊTES (Edges)
            // =========================================================
            const auto& edgeDefs = Cell::getLocalEdgesIndices(cell.type);
            for (const auto& edge : edgeDefs) {
                // Pour un élément P2, l'arête a 3 nœuds : {sommet1, sommet2, milieu}
                if (edge.size() < 3) continue; 
                
                const Vector3D& p1 = coords[edge[0]];
                const Vector3D& p2 = coords[edge[1]];
                const Vector3D& pm = coords[edge[2]]; // Nœud milieu

                Vector3D theoreticalMid = (p1 + p2) * 0.5;
                
                if ((pm - theoreticalMid).norm() > tol) {
                    cellIsAffine = false;
                    break; // Inutile de vérifier le reste, la cellule est courbée
                }
            }

            // =========================================================
            // VÉRIFICATION DES FACES QUADRILATÉRALES (Quad Faces)
            // =========================================================
            if (cellIsAffine) {
                const auto& quadFaceDefs = Cell::getLocalQuadFacesIndices(cell.type);
                for (const auto& face : quadFaceDefs) {
                    // Pour une face de Quad9, Prism18 ou Hex27, on a 9 nœuds.
                    // L'indice 8 est le centre de la face.
                    if (face.size() < 9) continue; 
                    
                    const Vector3D& p0 = coords[face[0]];
                    const Vector3D& p1 = coords[face[1]];
                    const Vector3D& p2 = coords[face[2]];
                    const Vector3D& p3 = coords[face[3]];
                    const Vector3D& pCenter = coords[face[8]]; // Nœud central

                    Vector3D theoreticalCenter = (p0 + p1 + p2 + p3) * 0.25;
                    
                    if ((pCenter - theoreticalCenter).norm() > tol) {
                        cellIsAffine = false;
                        break;
                    }
                }
            }

            // =========================================================
            // VÉRIFICATION DU VOLUME CENTRAL (Hex27 uniquement)
            // =========================================================
            if (cellIsAffine && cell.type == CellType::HEX27 && coords.size() == 27) {
                // il faut aussi s'assurer que le nœud 26 (centre du volume) est bien au milieu.
                Vector3D theoreticalVolumeCenter = Vector3D(0, 0, 0);
                for (int i = 0; i < 8; ++i) {
                    theoreticalVolumeCenter += coords[i];
                }
                theoreticalVolumeCenter = theoreticalVolumeCenter * 0.125;

                if ((coords[26] - theoreticalVolumeCenter).norm() > tol) {
                    cellIsAffine = false;
                }
            }

            m_isCellAffine[c] = cellIsAffine;
            m_isGloballyAffine = m_isGloballyAffine && cellIsAffine;
        }
    }

    //------------------------------------------------------------------------------------
    //------------------------- HELPER FUMCTIOMS -----------------------------------------
    //------------------------------------------------------------------------------------

    /** ************************************************************************************************************************************
     * @brief Convert a cell type enum to a human-readable string.
     *
     * This utility is used for reporting, logging, and debugging. It covers all standard linear and quadratic types.
     *
     * @param type The cell type (CellType enum).
     * @return A string representing the cell type, e.g., "TET4" or "QUAD9".
     *         Returns "UNKNOWN_CELL_TYPE" if the type is unrecognized.
     ***************************************************************************************************************************************/
    std::string Mesh::cellTypeToString(CellType type) const
    {
        switch (type)
        {
        case CellType::LINE2:   return "LINE2";
        case CellType::TRI3:    return "TRI3";
        case CellType::QUAD4:   return "QUAD4";
        case CellType::TET4:    return "TET4";
        case CellType::PRISM6:  return "PRISM6";
        case CellType::HEX8:    return "HEX8";
        case CellType::LINE3:   return "LINE3";
        case CellType::TRI6:    return "TRI6";
        case CellType::QUAD9:   return "QUAD9";
        case CellType::TET10:   return "TET10";
        case CellType::PRISM18: return "PRISM18";
        case CellType::HEX27:   return "HEX27";
        // ... others?
        default:
            return "UNKNOWN_CELL_TYPE";
        }
    }


    /** ************************************************************************************************************************************
     * @brief Compute a quality metric for a specific cell.
     *
     * Quality is a dimensionless number between 0 and 1, where 1 represents
     * an ideal element (perfectly shaped) and 0 represents a degenerate element
     * (zero area/volume or inverted orientation).
     *
     * Supported types:
     * - 2D: TRI3, TRI6, QUAD4, QUAD9
     * - 3D: TET4, TET10, PRISM6, PRISM18
     *
     * @param cellId The global index of the cell.
     * @return A real number in [0,1] measuring the element quality.
     ***************************************************************************************************************************************/
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

    /** ************************************************************************************************************************************
     * @brief Compute quality of a quadrilateral element (QUAD4/QUAD9).
     *
     * The quadrilateral is split into two triangles, and the quality is averaged.
     *
     * @param v A span of vertex IDs for the quad (corners first, then mid-edge nodes for quadratic).
     * @return Quality metric in [0,1]. Returns 0 if vertex count < 4.
     ***************************************************************************************************************************************/
    CfemReal Mesh::computeQuadQuality(const std::span<const CfemInt>& v) const
    {
        if (v.size() < 4) return 0.0;

        CfemInt tri1[3] = { v[0], v[1], v[2] };
        CfemInt tri2[3] = { v[0], v[2], v[3] };

        CfemReal q1 = computeTriQuality(tri1);
        CfemReal q2 = computeTriQuality(tri2);
        
        return (q1 + q2) / 2.0;
    }

    /** ************************************************************************************************************************************
     * @brief Compute quality of a prism element (PRISM6/PRISM18).
     *
     * The prism is subdivided into three tetrahedra. The quality is the
     * minimum of the three tetrahedra.
     *
     * @param v A span of vertex IDs for the prism.
     * @return Quality metric in [0,1]. Returns 0 if vertex count < 4.
     ***************************************************************************************************************************************/
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

    /** ************************************************************************************************************************************
     * @brief Compute quality of a triangular element (TRI3/TRI6).
     *
     * Quality is based on edge lengths and triangle shape:
     * ideal equilateral triangle has quality 1, degenerate triangles 0.
     *
     * @param v A span of 3 vertex IDs (for linear) or 6 vertex IDs (for quadratic).
     * @return Quality metric in [0,1].
     ***************************************************************************************************************************************/
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

    /** ************************************************************************************************************************************
     * @brief Compute quality of a tetrahedral element (TET4/TET10).
     *
     * Quality metric is based on volume and edge lengths.
     *
     * @param v A span of 4 vertex IDs (linear) or 10 vertex IDs (quadratic).
     * @return Quality metric in [0,1]. Returns 0 for inverted or degenerate tetrahedra.
     ***************************************************************************************************************************************/
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

    /** ************************************************************************************************************************************
     * @brief Get the topological dimension of a given cell type.
     *
     * @param type The cell type (CellType enum).
     * @return Dimension as integer:
     *         1 for lines, 2 for quads/triangles, 3 for tetra/prism/hex,
     *         0 for unknown/unsupported type.
     ***************************************************************************************************************************************/
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

    /** ************************************************************************************************************************************
     * @brief Determine the polynomial order of a cell type.
     *
     * @param type The cell type (CellType enum).
     * @return MeshOrder::Linear, MeshOrder::Quadratic, or MeshOrder::None.
     ***************************************************************************************************************************************/
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

    /** ************************************************************************************************************************************
     * @brief Check if a user-defined entity name is valid.
     *
     * Reserved names "Universe" and "SYS_BOUNDARY" are not allowed for
     * user entities. This helps avoid conflicts with internal mesh entities.
     *
     * @param name The proposed entity name.
     * @throw CFEM_ERROR if the name is reserved.
     ***************************************************************************************************************************************/
    void Mesh::validateEntityName(const std::string &name) const
    {
        if (name == "Universe" || name == "SYS_BOUNDARY")
        {
            CFEM_ERROR << "Reserved Name Violation: The name '" << name <<
                       "' is reserved for internal mesh topology. " <<
                       "Please use a different name for the physical group.";
        }
    }


    /** ************************************************************************************************************************************
     * @name Built-in Mesh Generators
     *
     * Factory functions that generate common computational domains.
     *
     * These methods construct a mesh, populate vertices and cells,
     * build connectivity structures, and automatically create
     * boundary entities.
     *
     * Generated meshes are immediately ready for use in finite
     * element computations.
     *
     * Supported geometries include:
     *
     * - Rectangular domains (`quad`)
     * - 3D blocks (`block`)
     * - Circular domains (`disk`)
     * - Annular regions (`annulus`)
     * - Spheres (`sphere`)
     * - Cylinders (`cylinder`)
     * - Hollow cylinders (`hollow_cylinder`)
     *
     * All generators support linear and quadratic meshes through
     * the `MeshOrder` parameter.
     *
     * Example:
     *
     * @code
     * auto mesh = Mesh::block(
     *     MeshOrder::Linear,
     *     Mesh::BlockSubdivisions{10,10,10},
     *     Mesh::BlockBounds({0,0,0},{1,1,1}),
     *     Mesh::BlockElementType::Tetrahedron
     * );
     * @endcode
     *
     * @{
     *************************************************************************************************************************************************************************************************************************************/

     /** ************************************************************************************************************************************
     * @brief Generate a rectangular mesh.
     *
     * Creates a structured two-dimensional mesh over a rectangular
     * domain. The domain can be defined either by axis-aligned bounds
     * or by specifying the coordinates of its four corners.
     *
     * @param order Polynomial order of the mesh.
     * @param subDiv Number of subdivisions along the x and y directions.
     * @param bounds Rectangular domain bounds.
     * @param elemType Element type used to discretize the domain.
     * @param names Optional boundary entity names.
     *
     * @return A newly created mesh.
     *
     * @note Boundary entities are automatically created for the
     *       domain edges (`xmin`, `xmax`, `ymin`, `ymax`).
     ***************************************************************************************************************************************/
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
            CFEM_ERROR << "QuadBounds must collapse exactly one dimension.";

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

    /** ************************************************************************************************************************************
     * @brief Generate a rectangular mesh using default boundary names.
     *
     * Same as the full overload but automatically assigns the
     * default boundary entity names.
     ***************************************************************************************************************************************/
    std::shared_ptr<Mesh> Mesh::quad(MeshOrder order,
                                        const QuadSubdivisions& subDiv,
                                        const QuadBounds& bounds,
                                        QuadElementType elemType)
    {
        return quad(order, subDiv, bounds, elemType, QuadBoundaryNames{});
    }

    /** ************************************************************************************************************************************
     * @brief Generate a quadrilateral mesh using explicit corner coordinates.
     *
     * This version allows defining non-axis-aligned quadrilateral domains
     * by specifying the coordinates of the four corners.
     ***************************************************************************************************************************************/
    std::shared_ptr<Mesh> Mesh::quad(MeshOrder order,
                                        const QuadSubdivisions& subDiv,
                                        const QuadCorners& corners,
                                        QuadElementType elemType,
                                        const QuadBoundaryNames& names)
    {
        return QuadMeshBuilder::create(order, subDiv, corners, elemType, names);
    }

    /** ************************************************************************************************************************************
     * @brief Generate a quadrilateral mesh using default boundary names.
     *
     * Same as the full overload but automatically assigns the
     * default boundary entity names.
     ***************************************************************************************************************************************/
    std::shared_ptr<Mesh> Mesh::quad(MeshOrder order,
                                        const QuadSubdivisions& subDiv,
                                        const QuadCorners& corners,
                                        QuadElementType elemType)
    {
        return quad(order, subDiv, corners, elemType, QuadBoundaryNames{});
    }


    /** ************************************************************************************************************************************
     * @brief Generate a circular mesh.
     *
     * Creates a two-dimensional mesh of a disk centered at the
     * specified point.
     *
     * @param order Polynomial order of the mesh.
     * @param mesh_size Target element size.
     * @param center Disk center.
     * @param radius Disk radius.
     * @param normal Optional normal vector defining the disk plane.
     * @param elemtype Element type used in the mesh.
     *
     * @return A newly created mesh.
     *
     * @note The boundary of the disk is automatically registered
     *       as a mesh entity.
     ***************************************************************************************************************************************/
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

    /** ************************************************************************************************************************************
     * @brief Generate an annular mesh.
     *
     * Creates a two-dimensional mesh representing a ring-shaped
     * domain between two concentric circles.
     *
     * @param order Polynomial order of the mesh.
     * @param h Target mesh size.
     * @param center Center of the annulus.
     * @param radius_in Inner radius.
     * @param radius_out Outer radius.
     * @param type Element type used for discretization.
     * @param names Optional boundary entity names.
     *
     * @return A newly created mesh.
     *
     * @note Two boundary entities are generated:
     *       the inner and outer circles.
     ***************************************************************************************************************************************/
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

    /** ************************************************************************************************************************************
     * @brief Generate a spherical mesh.
     *
     * Creates a three-dimensional mesh of a sphere.
     *
     * @param order Polynomial order of the mesh.
     * @param mesh_size Target mesh size.
     * @param center Sphere center.
     * @param radius Sphere radius.
     * @param elemtype Element type used for discretization.
     *
     * @return A newly created mesh.
     ***************************************************************************************************************************************/
    std::shared_ptr<Mesh> Mesh::sphere(
                                    MeshOrder order,
                                    CfemReal mesh_size,
                                    const Vector3D &center,
                                    CfemReal radius,
                                    SolidSphereElementType elemtype)
    {
        return SolidSphereBuilder::create(order, mesh_size, center, radius, elemtype);
    }

    std::shared_ptr<Mesh> Mesh::surface_sphere(
                                            MeshOrder order, 
                                            CfemReal mesh_size, 
                                            const Vector3D& center, 
                                            CfemReal radius
                                        )
    {
        return SurfaceSphereBuilder::create(order, mesh_size, center, radius );

    }


    /** ************************************************************************************************************************************
     * @brief Generate a 3D block mesh.
     *
     * Creates a structured three-dimensional mesh over a rectangular
     * prism domain.
     *
     * @param order Polynomial order of the mesh.
     * @param subDiv Number of subdivisions along each axis.
     * @param bounds Axis-aligned bounding box of the domain.
     * @param elemType Element type used for discretization.
     * @param names Optional boundary entity names.
     *
     * @return A newly created mesh.
     *
     * @note Boundary entities are automatically created for the
     *       six faces of the block (`xmin`, `xmax`, `ymin`,
     *       `ymax`, `zmin`, `zmax`).
     ***************************************************************************************************************************************/
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
            CFEM_ERROR << "BlockBounds must define a non-degenerate 3D volume.";
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

    /** ************************************************************************************************************************************
     * @brief Generate a 3D block mesh.
     *
     * Creates a structured three-dimensional mesh over a rectangular
     * prism domain.
     *
     * @param order Polynomial order of the mesh.
     * @param subDiv Number of subdivisions along each axis.
     * @param bounds Axis-aligned bounding box of the domain.
     * @param elemType Element type used for discretization.
     *
     * @return A newly created mesh.
     *
     * @note Boundary entities are automatically created for the
     *       six faces of the block using the default boundary names (`xmin`, `xmax`, `ymin`,
     *       `ymax`, `zmin`, `zmax`).
     ***************************************************************************************************************************************/
    std::shared_ptr<Mesh> Mesh::block(MeshOrder order,
                                       const BlockSubdivisions& subDiv, 
                                       const BlockBounds& bounds,
                                       BlockElementType elemType)
    {
        return block(order, subDiv, bounds, elemType, BlockBoundaryNames{});
    }


    /** ************************************************************************************************************************************
     * @brief Generate a block mesh using explicit corner coordinates.
     *
     * This version allows the block to be defined using the coordinates
     * of its eight vertices rather than axis-aligned bounds.
     *
     * @param order Polynomial order of the mesh.
     * @param subDiv Number of subdivisions along each axis.
     * @param corners The eight vertices coordinates.
     * @param elemType Element type used for discretization.
     * @param names Optional boundary entity names.
     *
     * @return A newly created mesh.
     *
     * @note Boundary entities are automatically created for the
     *       six faces of the block using user-provided boundary names.
     ***************************************************************************************************************************************/
    std::shared_ptr<Mesh> Mesh::block(MeshOrder order,
                                       const BlockSubdivisions& subDiv, 
                                       const BlockCorners& corners,
                                       BlockElementType elemType,
                                       const BlockBoundaryNames& names)
    {
        return BlockMeshBuilder::create(order, subDiv, corners, elemType, names);

    }

    /** ************************************************************************************************************************************
     * @brief Generate a block mesh using explicit corner coordinates.
     *
     * This version allows the block to be defined using the coordinates
     * of its eight vertices rather than axis-aligned bounds.
     *
     * @param order Polynomial order of the mesh.
     * @param subDiv Number of subdivisions along each axis.
     * @param corners The eight vertices coordinates.
     * @param elemType Element type used for discretization.
     *
     * @return A newly created mesh.
     *
     * @note Boundary entities are automatically created for the
     *       six faces of the block using the default boundary names (`xmin`, `xmax`, `ymin`,
     *       `ymax`, `zmin`, `zmax`).
     ***************************************************************************************************************************************/
    std::shared_ptr<Mesh> Mesh::block(MeshOrder order,
                                       const BlockSubdivisions& subDiv, 
                                       const BlockCorners& corners,
                                       BlockElementType elemType)
    {
        return BlockMeshBuilder::create(order, subDiv, corners, elemType, BlockBoundaryNames{});
    }

    /** ************************************************************************************************************************************
     * @brief Generate a cylindrical mesh.
     *
     * Creates a mesh representing a solid cylinder.
     *
     * @param order Polynomial order of the mesh.
     * @param target_mesh_size Desired element size.
     * @param center Center of the cylinder base.
     * @param radius Cylinder radius.
     * @param height Cylinder height.
     * @param axis Direction of the cylinder axis.
     * @param type Element type used for discretization.
     * @param names Optional boundary entity names.
     *
     * @return A newly created mesh.
     *
     * @note Boundary entities are created for the lateral surface
     *       and the top and bottom caps.
     ***************************************************************************************************************************************/
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


    /** ************************************************************************************************************************************
     * @brief Generate a hollow cylindrical mesh.
     *
     * Creates a mesh representing a cylindrical shell defined by
     * inner and outer radii.
     *
     * @param order Polynomial order of the mesh.
     * @param target_mesh_size Desired element size.
     * @param center Cylinder center.
     * @param r_in Inner radius.
     * @param r_out Outer radius.
     * @param height Cylinder height.
     * @param axis Direction of the cylinder axis.
     * @param type Element type used for discretization.
     * @param names Optional boundary entity names.
     *
     * @return A newly created mesh.
     *
     * @note Boundary entities are created for the inner surface,
     *       outer surface, and cylinder caps.
     ***************************************************************************************************************************************/
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
    /// @}

} // namespace