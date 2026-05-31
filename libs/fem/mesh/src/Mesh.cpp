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

#include <algorithm>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>

#include "AnnulusMeshBuilder.h"
#include "BlockMeshBuilder.h"
#include "CylinderMeshBuilder.h"
#include "DiskMeshBuilder.h"
#include "HollowCylinderMeshBuilder.h"
#include "LineMeshBuilder.h"
#include "Mesh.h"
#include "MeshConformityChecker.h"
#include "QuadMeshBuilder.h"
#include "SolidSphereBuilder.h"
#include "SurfaceSphereBuilder.h"

namespace cfem
{

    Mesh::Mesh(ConstructorToken ltok)
    {
        m_cellCountPerType.fill(0);
        m_faceCountPerType.fill(0);
        CFEM_CHECK_INVARIANTS();
    }

    Mesh::Mesh(ConstructorToken lTok, CfemInt p_topologicalDimension, MeshOrder p_meshOrder) : m_topologicalDimension(p_topologicalDimension), m_meshOrder(p_meshOrder)
    {
        CFEM_PRECONDITION(p_topologicalDimension >= 1 && p_topologicalDimension <= 3);
        CFEM_POSTCONDITION(m_topologicalDimension == p_topologicalDimension);
        CFEM_POSTCONDITION(m_meshOrder == p_meshOrder);

        m_cellCountPerType.fill(0);
        m_faceCountPerType.fill(0);

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

        // update the cell count for this type
        ++m_cellCountPerType[static_cast<size_t>(type)];

        // Manage unique element types
        CfemInt typeIdx = static_cast<int>(type);
        if (!m_cellTypeFlags.test(typeIdx))
        {
            m_cellTypeFlags.set(typeIdx);
            m_uniqueCellTypes.push_back(type);
        }


        auto facetTypes = Cell::getFacetTypes(type); 
        
        for (CellType fType : facetTypes)
        {
            CfemInt fTypeIdx = static_cast<int>(fType);
            
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
            CFEM_ERROR << "Cannot add " << Cell::cellTypeToString(type) <<
                       " to a " << (m_meshOrder == MeshOrder::Linear ? "linear mesh." 
                                 : m_meshOrder == MeshOrder::Quadratic ? "quadratic mesh."
                                                                       : "mesh with unkmown order.");
        }

        CfemInt cellDim = getDimFromType(type);
        if (cellDim != m_topologicalDimension)
        {
            CFEM_ERROR << "Cannot add a " << std::to_string(cellDim) << "D element (" <<
                       Cell::cellTypeToString(type) << ") to a " << std::to_string(m_topologicalDimension) << "D mesh.";
        }

        if (v_ids.size() != getNodeCountPerCell(type))
        {
            CFEM_ERROR << "Topology error: " << Cell::cellTypeToString(type) <<
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
                CfemInt half = getNodeCountPerCell(type) / 2;
                for(CfemInt i=0; i<half; ++i) std::swap(vids[i], vids[i+half]);
                break;
            }
            case CellType::HEX27:{
                // Swap Corners: Bottom (0-3) with Top (4-7)
                for (CfemInt i = 0; i < 4; ++i) {
                    std::swap(vids[i], vids[i + 4]);
                }

                // Swap Horizontal Edges: Bottom (8-11) with Top (12-15)
                for (CfemInt i = 8; i < 12; ++i) {
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
                for (CfemInt i = 0; i < 3; ++i) std::swap(vids[i], vids[i + 3]);

                // Swap Horizontal Edges: Bottom (6,7,8) with Top (9,10,11)
                for (CfemInt i = 6; i < 9; ++i) std::swap(vids[i], vids[i + 3]);
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
     * @brief Builds the cell-to-cell adjacency graph based on shared facets.
     *
     * @details
     * Populates highly optimized, contiguous internal arrays (CSR format) to map each 
     * cell to its immediate neighbors and bounding facets. This establishes the spatial 
     * relationship between elements, which is strictly required for algorithms involving 
     * face-based flux computations, gradient reconstructions, or mesh smoothing.
     * 
     * Internal `Cell` structures are updated with direct pointers to this contiguous 
     * memory, guaranteeing $O(1)$ cache-friendly neighbor lookups during assembly.
     * Boundary faces are explicitly designated by a neighbor ID of `-1`.
     *
     * @pre `buildFacets()` MUST be executed prior to this function, as the adjacency 
     * logic strictly relies on the left/right cell ownership established by the facets.
     * 
     * @warning This is a destructive operation. It clears and rewrites all neighbor 
     * mapping arrays, invalidating any previously cached adjacency pointers.
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
            totalSlots += Cell::getFaceCount(m_cells[i].type);
        }
        m_neighborOffsets[numCells] = totalSlots; 

        m_allCellNeighbors.assign(totalSlots, -1);
        m_allCellFacetIds.assign(totalSlots, -1);

        for (const auto& f : m_facets) 
        {
            auto f_nodes = getFacetConnectivity(f.id);

            auto getLocalFaceIdx = [&](CfemInt cell_id) -> CfemInt {
                CfemInt num_faces = Cell::getFaceCount(m_cells[cell_id].type);
                for (CfemInt i = 0; i < num_faces; ++i) {
                    LocalFace currentLocal = getLocalFaceData(cell_id, i);
                    if (static_cast<size_t>(currentLocal.size) == f_nodes.size()) {
                        bool match = true;
                        for(size_t fn : f_nodes) {
                            bool found = false;
                            for(CfemInt k=0; k < currentLocal.size; ++k) {
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
            CfemInt localIdxL = getLocalFaceIdx(f.left_cell);
            if (localIdxL != -1) {
                size_t pos = m_neighborOffsets[f.left_cell] + localIdxL;
                m_allCellNeighbors[pos] = f.right_cell;
                m_allCellFacetIds[pos] = f.id;
            }

            // Link Right Cell
            if (!f.isBoundary()) {
                CfemInt localIdxR = getLocalFaceIdx(f.right_cell);
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
            cell.num_faces = static_cast<CfemInt>(Cell::getFaceCount(cell.type));
            
            // Point to the data inside the Mesh's big blocks
            cell.neighbor_cells_ptr = &m_allCellNeighbors[startIdx];
            cell.facet_ids_ptr     = &m_allCellFacetIds[startIdx];
        }
    }

    /** ************************************************************************************************************************************
     * @brief Detects and constructs all unique geometric facets (edges/faces) in the mesh.
     *
     * @details
     * Scans the connectivity of all elements to extract a globally unique set of internal 
     * interfaces and boundary surfaces. For each discovered facet, the algorithm determines 
     * its topological ownership, assigning the "left" and "right" sharing cells (or marking 
     * it as a boundary if unshared).
     * 
     * To maximize HPC performance during boundary and surface integrations, all facet 
     * connectivities are packed into a single, contiguous CSR memory block rather than 
     * individual heap allocations.
     *
     * @pre `buildInverseConnectivities()` MUST be called before this function. The facet-matching 
     * algorithm uses the `Vertex -> Cells` inverse map to perform ultra-fast, localized 
     * neighbor searches instead of global $O(N^2)$ checks.
     * 
     * @warning This completely destroys and recreates the global facet list. Any external 
     * integer IDs referencing facets, or `std::span` views of facet connectivities acquired 
     * prior to this call, will become dangling or invalid.
     *
     ***************************************************************************************************************************************/
    void Mesh::buildFacets() 
    {
        PROFILE_SCOPE("Mesh::buildFacets");

        // RESET & RESERVE
        m_facets.clear();
        m_allFacetVertexIds.clear();
        m_facetOffsets.clear();
        m_faceCountPerType = {};

        size_t est_facets = m_cells.size() * 4; 
        m_allFacetVertexIds.reserve(est_facets * 4); 
        m_facetOffsets.reserve(est_facets + 1);
        m_facetOffsets.push_back(0); 
        m_facets.reserve(est_facets);

        std::vector<CfemInt> cell_face_to_facet(m_cells.size() * 6, -1);

        // BUILD TOPOLOGY
        for (CfemInt c_id = 0; c_id < static_cast<CfemInt>(m_cells.size()); ++c_id) 
        {
            const CfemInt num_faces = Cell::getFaceCount(m_cells[c_id].type);
            
            for (CfemInt f_idx = 0; f_idx < num_faces; ++f_idx) 
            {
                if (cell_face_to_facet[c_id * 6 + f_idx] != -1) continue;

                LocalFace current_face = getLocalFaceData(c_id, f_idx);
                
                // Neighbor Search: Find Best Anchor
                CfemInt v_anchor = current_face.v[0];
                CfemInt min_degree = m_vertices[v_anchor].num_connected_cells;
                
                for (CfemInt k = 1; k < current_face.size; ++k) {
                    CfemInt candidate = current_face.v[k];
                    CfemInt degree = m_vertices[candidate].num_connected_cells;
                    if (degree < min_degree) {
                        min_degree = degree;
                        v_anchor = candidate;
                    }
                }

                const auto& anchor_vert = m_vertices[v_anchor];
                CfemInt neighbor_id = -1;
                CfemInt neighbor_f_idx = -1;

                for (CfemInt i = 0; i < anchor_vert.num_connected_cells; ++i) 
                {
                    CfemInt other_c_id = anchor_vert.connected_cells_ptr[i];
                    if (other_c_id == c_id) continue;

                    const CfemInt other_num_faces = Cell::getFaceCount(m_cells[other_c_id].type);
                    for (CfemInt j = 0; j < other_num_faces; ++j) {
                        LocalFace other_face = getLocalFaceData(other_c_id, j);
                        if (current_face.matches(other_face)) {
                            neighbor_id = other_c_id;
                            neighbor_f_idx = j;
                            goto found_neighbor; 
                        }
                    }
                }
                found_neighbor:;

                // Store Connectivity in Big Block (Bulk Insert)
                size_t start_offset = m_allFacetVertexIds.size();
                m_allFacetVertexIds.insert(m_allFacetVertexIds.end(), current_face.v, current_face.v + current_face.size);

                // Construction in m_facets
                CfemInt new_facet_id = static_cast<CfemInt>(m_facets.size());
                
                m_facets.emplace_back();
                Facet& f = m_facets.back();
                f.id = new_facet_id;
                f.type = current_face.type;
                f.left_cell = c_id;
                f.right_cell = neighbor_id;
                f.left_cell_facet_id = f_idx;
                f.right_cell_facet_id = neighbor_f_idx;
                f.offset = start_offset;
                f.node_count = static_cast<size_t>(current_face.size);

                ++m_faceCountPerType[static_cast<size_t>(f.type)];

                // Update mapping
                cell_face_to_facet[c_id * 6 + f_idx] = new_facet_id;
                if (neighbor_id != -1) {
                    cell_face_to_facet[neighbor_id * 6 + neighbor_f_idx] = new_facet_id;
                }

                m_facetOffsets.push_back(m_allFacetVertexIds.size());
            }
        }
    }

    /** ************************************************************************************************************************************
     * @name Mesh Entities
     *
     * Mesh entities map physical groups (e.g., "inlet", "wall", "steel_region") to 
     * specific topological subsets of the mesh. They act as the bridge between the 
     * raw geometry and the finite element formulation (Boundary Conditions, Materials).
     * @{
     ***************************************************************************************************************************************/

     /** ************************************************************************************************************************************
     * @brief Registers a new user-defined physical group in the mesh.
     *
     * @details
     * Creates a named entity of a specific topological dimension (0=Vertex, 1=Edge, 
     * 2=Face, 3=Volume). The mesh automatically safely assigns a unique internal ID 
     * that is guaranteed not to collide with critical system entities (like `Universe` 
     * or `SYS_BOUNDARY`).
     *
     * @param name Unique string identifier (must pass `validateEntityName()` policy).
     * @param dim Topological dimension of the geometry it represents.
     * @return The auto-generated safe integer ID for this entity.
     *
     ***************************************************************************************************************************************/
    CfemInt Mesh::createEntity(const std::string &name, CfemInt dim)
    {
        validateEntityName(name);

        CfemInt newId = static_cast<CfemInt>(m_entities.size());
        if (getEntity(cfem::UNIVERSE_ID) && getEntity(cfem::SYS_BOUNDARY_ID))
        {
            newId -= 2;
        }
        else if (getEntity(cfem::SYS_BOUNDARY_ID) || getEntity(cfem::UNIVERSE_ID))
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
     * @brief Retrieves a non-owning pointer to an entity by its name.
     * @param name Name of the entity.
     * @return Pointer to the corresponding entity, or `nullptr` if no entity with this name exists.
     *
     * @note Performs a linear search `O(N)` over registered entities. Since the number 
     * of physical groups is typically very small, this is fast enough for initialization 
     * phases. Do not call this inside hot assembly loops.
     * 
     * @warning
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

    /** ************************************************************************************************************************************
     * @brief Retrieves a const reference to an entity by its name.
     * @param name The name of the entity.
     * @return const reference to the corresponding entity.
     * @throw CFEM_ERROR  if no entity with this name exists
     *
     * @note Performs a linear search `O(N)` over registered entities. Since the number 
     * of physical groups is typically very small, this is fast enough for initialization 
     * phases. Do not call this inside hot assembly loops.
     * 
     ***************************************************************************************************************************************/
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
     * @brief Retrieves a non-owning pointer to an entity by its ID.
     * @param id Internal ID of the entity.
     * @return Pointer to the corresponding entity, or `nullptr` if no entity with this ID exists.
     *
     * @note Performs a linear search `O(N)` over registered entities. Since the number 
     * of physical groups is typically very small, this is fast enough for initialization 
     * phases. Do not call this inside hot assembly loops.
     * 
     * @warning
     * The returned pointer is owned by the mesh and must not be deleted.
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

    /** ************************************************************************************************************************************
     * @brief Retrieves a const reference to an entity by its ID.
     * @param id Internal ID of the entity.
     * @return const reference to the corresponding entity.
     * @throw CFEM_ERROR  if no entity with this ID exists
     *
     * @note Performs a linear search `O(N)` over registered entities. Since the number 
     * of physical groups is typically very small, this is fast enough for initialization 
     * phases. Do not call this inside hot assembly loops.
     * 
     ***************************************************************************************************************************************/
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
     * @brief Internal factory for system-critical entities.
     *
     * @details
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
     * @warning **DO NOT USE for physical groups.** This bypasses collision safety checks 
     * and forces a specific ID. It is strictly reserved for the core mesh builder to 
     * construct the `Universe` and `SYS_BOUNDARY` singletons.
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
     * @brief This method switches the mesh between linear and quadratic representations.
     *
     * @details
     * Acts as a high-level wrapper. If the mesh is populated, this immediately triggers 
     * a heavy topological rebuild (`upgradeToQuadratic()` or `downgradeToLinear()`).
     *
     * If the mesh is empty, the order is simply updated.
     *
     * @param newOrder Desired mesh order (`Linear` or `Quadratic`).
     * 
     * @warning Modifying the order of a populated mesh is a destructive operation. 
     * It invalidates all previously acquired connectivity spans and automatically 
     * triggers updates on all registered FE Space observers.
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
     * @brief Sets the spatial dimension (1D, 2D, or 3D) of the domain.
     *
     * @pre **The mesh MUST be completely empty.** Calling this after vertices or 
     * cells have been added will trigger a fatal assertion, as it would render the 
     * existing topology geometrically invalid.
     *
     * @param newDim Mesh tological dimension (1, 2, or 3).
     *
     * @throws CFEM_ERROR if the mesh already contains elements.
     *
     * @note
     * Changing the dimension after elements have been added would invalidate the topology and is therefore prohibited.
     ***************************************************************************************************************************************/
    void Mesh::setTopologicalDimension(CfemInt newDim)
    {
        CFEM_ASSERT(newDim >= 1 || newDim <= 3);
        if (!(this->m_cells.empty()))
        {
            CFEM_ERROR << "ARE YOU TRYING TO SET THE DIMENSION OF A MESH THAT ALREADY CONTAINS ELEMENTS? MAKE SURE THE MESH IS EMPTY.";
        }
        m_topologicalDimension = newDim;
    }

    /** ************************************************************************************************************************************
     * @brief Defines a global reference normal for 2D manifolds embedded in 3D.
     * 
     * @details
     * Crucial for surface meshes (e.g., a sphere surface or a 2D sheet in 3D space) 
     * to ensure consistent winding rules and positive signed areas across all elements.
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
     * @brief Registers a Finite Element space to listen for topological mutations.
     *
     * @details
     * Uses a weak-referencing Observer pattern. The mesh does not take ownership of 
     * the FE space, allowing spaces to be safely destroyed without dangling pointer 
     * issues. 
     * 
     * If the mesh is already in a non-default state (e.g., Quadratic) when the 
     * observer is added, the observer is immediately notified to synchronize its DOFs.
     *
     * @param space Shared pointer to the finite element space observing the mesh.
     *
     * @note
     * The mesh only stores weak references to observers internally, allowing finite element spaces
     * to be destroyed safely without requiring explicit deregistration.
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

    /** ************************************************************************************************************************************
     * @brief Broadcasts a topology change event to all living observers.
     * 
     * @details
     * Automatically performs garbage collection during the broadcast: any observer 
     * (FE Space) that has been destroyed by the user is silently removed from the 
     * notification list.
     ***************************************************************************************************************************************/
    void Mesh::notifyObservers() {
        for (auto it = m_observers.begin(); it != m_observers.end();)
        {
            if (auto obs = it->lock()) { // On vérifie si l'observateur existe encore
                obs->updateAfterMeshChange();
                ++it;
            }
            else {
                it = m_observers.erase(it); // Nettoyage automatique
            }
        }
    }

    /** ************************************************************************************************************************************
     * @brief Finalizes mesh construction by orchestrating the complete topological rebuild.
     *
     * @details
     * Acts as the main compiler for the mesh data structure. It takes the raw, unstructured 
     * list of vertices and cells, and weaves the complete topological web required for 
     * finite element operations.
     * 
     * This method executes the "Golden Sequence" of topology construction:
     * 1. Garbage collects unreferenced vertices.
     * 2. Builds `Vertex -> Cell` inverse mappings.
     * 3. Extracts unique boundary and internal Facets.
     * 4. Establishes `Cell -> Cell` CSR adjacency.
     * 
     * Additionally, it automatically generates and populates the core system entities: 
     * `Universe` (containing all elements) and `SYS_BOUNDARY` (containing all unshared 
     * external facets and their associated vertices).
     *
     * ### Operational Details
     * - **Safe Execution:** Triggers `updateAffineFlags()` to cache Jacobian mapping 
     *   properties for the integration engine.
     * - **Validation:** In debug builds (`!NDEBUG`), automatically runs the 
     *   `MeshConformityChecker` to catch degenerate elements or twisted topologies 
     *   before they crash the solver.
     *
     * @pre The mesh must possess a defined dimension (`m_topologicalDimension > 0`) and contain 
     * at least one cell and one vertex.
     * @warning This is a heavy, destructive operation. It will completely overwrite 
     * any existing topological maps, adjacency structures, and system entity lists. 
     * It should typically be called only once after the initial geometric generation 
     * or after a major topological mutation.
     *
     * @complexity
     * Approximately linear in the number of mesh cells and facets.
     ***************************************************************************************************************************************/
    void Mesh::buildMesh()
    {
        if (m_cells.empty() || m_vertices.size() == 0 || m_topologicalDimension == 0) {
            CFEM_ERROR << "buildMesh failed: Mesh is empty.";
            return;
        }

        // Build topological structures (The Golden Sequence)
        compactVertices();
        buildInverseConnectivities();
        buildFacets(); 
        buildAdjacency();

        // Initialize Core Entities
        MeshEntity* universe = getEntity(cfem::UNIVERSE_ID);
        if (!universe) universe = &createInternalEntity(cfem::UNIVERSE_ID, "Universe", m_topologicalDimension);

        MeshEntity* boundary = getEntity(cfem::SYS_BOUNDARY_ID);
        if (!boundary) boundary = &createInternalEntity(cfem::SYS_BOUNDARY_ID, "SYS_BOUNDARY", m_topologicalDimension - 1);

        // Populate Universe
        if (universe->isEmpty()) {
            universe->reserveVertices(m_vertices.size());
            universe->reserveCells(m_cells.size());
            for (CfemInt i = 0; i < (CfemInt)m_vertices.size(); ++i) universe->addVertex(i);
            for (CfemInt i = 0; i < (CfemInt)m_cells.size(); ++i)    universe->addCell(i);
        }

        // Populate Facets and Identify Boundary
        // Use a bitset/bool vector to avoid massive duplicate insertions into 'boundary'
        std::vector<uint8_t> vertex_on_boundary(m_vertices.size(), 0);

        for (const auto &facet : m_facets) {
            universe->addFacet(facet.id);

            if (facet.isBoundary()) {
                boundary->addFacet(facet.id);
                auto f_vids = getFacetConnectivity(facet.id); 
                for (CfemInt v_id : f_vids) {
                    vertex_on_boundary[v_id] = 1;
                }
            }
        }

        // Single-pass addition of boundary nodes
        for (CfemInt v_id = 0; v_id < (CfemInt)m_vertices.size(); ++v_id) {
            if (vertex_on_boundary[v_id]) boundary->addVertex(v_id);
        }

        // Finalize
        universe->finalize();
        boundary->finalize();        
        updateAffineFlags();

        // Validation and Reporting
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
                details += "\n  - " + Cell::cellTypeToString(type) + ": " + std::to_string(count);
            }

            CFEM_INFO << "Orientation Check: " << std::to_string(totalFixed) << 
                    " inverted cells found and successfully fixed." << details;
            
            m_invertedCellCounts.clear();
        }
    }
    
    /** ************************************************************************************************************************************
     * @brief Performs garbage collection on the mesh by removing unreferenced vertices.
     *
     * @details
     * Analyzes the active cell connectivity to identify and discard orphaned nodes. 
     * It then automatically remaps all global identifiers across the entire mesh 
     * architecture, ensuring a contiguous, memory-tight Vertex array from 0 to N-1.
     *
     * ### Operational Details
     * - **Global Remapping:** Safely updates `m_allCellVertexIds`, `m_allFacetVertexIds`, 
     *   and the internal vertex lists of all physical `MeshEntities`.
     * - **Memory Locality:** Eliminating gaps in the vertex array significantly improves 
     *   cache hit rates during finite element matrix assembly.
     * - **HPC Optimized:** Uses byte-arrays instead of bit-masks for the marking phase 
     *   to ensure loop vectorization and maximum CPU throughput.
     *
     * @note This function solely cleans and remaps data; it does NOT automatically 
     * rebuild adjacency maps (like cell neighbors). The caller must invoke 
     * `buildInverseConnectivities()`, `buildFacets()`, and `buildAdjacency()` 
     * after compaction if the topology needs to be fully restored.
     *
     * @warning **Destructive ID Shift:** This operation shifts internal IDs. Any integer 
     * referencing a vertex ID acquired before calling this function will likely point 
     * to the wrong node (or out of bounds) afterward.
     *
     ***************************************************************************************************************************************/
    void Mesh::compactVertices()
    {
        const CfemInt oldVertexCount = static_cast<CfemInt>(m_vertices.size());
        if (oldVertexCount == 0) return;

        // Mark used vertices (Using uint8_t instead of bool for HPC memory speed)
        std::vector<uint8_t> used(oldVertexCount, 0);

        for (CfemInt cellId = 0; cellId < static_cast<CfemInt>(m_cells.size()); ++cellId)
        {
            const auto conn = getCellConnectivity(cellId);
            for (CfemInt vid : conn)
            {
                if (vid >= 0 && vid < oldVertexCount)
                    used[vid] = 1;
            }
        }

        // Count used vertices
        CfemInt usedCount = 0;
        for (uint8_t u : used) {
            if (u) ++usedCount;
        }

        // Early exit if no compaction is needed
        if (usedCount == oldVertexCount) return;

        // Build old -> new index map
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

        // Remap cell connectivity (flat storage)
        for (CfemInt& vid : m_allCellVertexIds) {
            if (vid >= 0) vid = oldToNew[vid];
        }

        // Remap facet connectivity (if stored explicitly)
        for (CfemInt& vid : m_allFacetVertexIds) {
            if (vid >= 0) vid = oldToNew[vid];
        }

        // Remap entity vertex sets
        for (const auto& entity : m_entities)
        {
            std::vector<CfemInt> newEntityVertices;
            newEntityVertices.reserve(entity->getVertexIds().size());

            for (CfemInt oldVid : entity->getVertexIds())
            {
                if (oldVid >= 0 && oldVid < oldVertexCount)
                {
                    CfemInt mapped = oldToNew[oldVid];
                    if (mapped >= 0) newEntityVertices.push_back(mapped);
                }
            }
            entity->replaceVertexList(newEntityVertices);
        }

        // Finalize
        m_vertices = std::move(newVertices);
    }


    /** ************************************************************************************************************************************
     * @brief Generates a concise, human-readable summary of the mesh.
     *
     * @details
     * Returns a formatted string detailing the global state of the mesh. This is 
     * highly useful for quick topological inspection, logging, and verifying the 
     * success of mesh transformations (e.g., confirming order upgrades).
     *
     * The summary includes:
     * - **Global Statistics:** Dimension, polynomial order, and total counts for 
     *   vertices, cells, and facets.
     * - **Cell Type Distribution:** A precise breakdown of the elements present 
     *   in the mesh (e.g., distinguishing the number of HEX27 vs. PRISM18).
     * - **Physical Entities:** A table of all registered mesh entities (Universe, 
     *   Boundaries, user groups), including their internal IDs and topological sizes.
     *
     * @return A formatted multi-line string containing the mesh summary.
     *
     * @note This method is `const` and relies entirely on pre-computed $O(1)$ 
     * metadata counters. It is computationally free and completely safe to call 
     * at any time, even within high-frequency simulation loops or debug assertions.
     *
     ***************************************************************************************************************************************/
    std::string Mesh::info() const
    {
        std::stringstream ss;
        ss << "\n================ MESH SUMMARY ================" << std::endl;
        ss << " Dimension: " << getTopologicalDimension() << "D" << std::endl;
        ss << " Order:     " << (getOrder() == MeshOrder::Linear ? "Linear" : 
                                (getOrder() == MeshOrder::Quadratic) ? "Quadratic" : "None") << std::endl;
        ss << " Vertices:  " << getNumVertices() << std::endl;
        ss << " Cells:     " << getNumCells() << std::endl;
        ss << " Facets:    " << getNumFacets() << std::endl;

        if (!m_uniqueCellTypes.empty()) 
        {
            ss << "---------------------------------------------------\n";
            ss << std::left << std::setw(20) << "Cell Type" 
               << "Count\n";
            for (CellType type : m_uniqueCellTypes) 
            {
                size_t count = m_cellCountPerType[static_cast<size_t>(type)];
                if (count > 0) 
                {
                    ss << std::left << std::setw(20) << Cell::cellTypeToString(type) 
                       << count << "\n";
                }
            }
        }

        ss << "---------------------------------------------------" << std::endl;
        ss << std::left << std::setw(20) << "Entity Name"
           << std::setw(10) << "ID"
           << std::setw(15) << "Elements"
           << "Vertices" << std::endl;

        for (const auto &entity : getEntities())
        {
            std::string name = entity->getName();
            size_t count = (entity->getTopologicalDimension() == getTopologicalDimension())
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
     * @brief Generates a detailed, human-readable string representation of the mesh topology.
     *
     * @details
     * Dumps the complete internal state of the mesh into a formatted text string. This 
     * is primarily designed for debugging, deep inspection, and verifying connectivity 
     * algorithms during development.
     *
     * The generated string is divided into four main sections:
     * - **1. Mesh Summary:** Global statistics and physical entity listings (via `info()`).
     * - **2. Nodes:** A list of all vertices and their physical 3D coordinates.
     * - **3. Cells:** Element types, node connectivity arrays, and adjacent neighbor IDs 
     *   (with 'B' indicating a boundary face).
     * - **4. Facets:** Edge/face connectivity and their Left/Right bounding cell IDs.
     *
     * @return A formatted multi-line string describing the full mesh connectivity.
     *
     * @note This is a diagnostic tool, not a standardized file exporter. To export 
     * meshes for external solvers or visualization software (e.g., ParaView), use 
     * dedicated writers (like VTK or Gmsh exporters) instead.
     * @warning The specific string layout and formatting is unstable by design and 
     * subject to change across library versions. Do not parse this output programmatically.
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
            << std::setw(8) << Cell::cellTypeToString(cell.type) << " "
            << cell_vids.size() << "  ";

            // Print Node IDs
            for (auto v_id : cell_vids) {
                ss << v_id << " ";
            }

            // Print Neighbors
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
     * @brief Upgrades the mesh from linear (P1) to quadratic (P2) order.
     *
     * @details
     * Converts all linear elements into their quadratic counterparts by inserting 
     * new vertices at edge midpoints and face centers.
     * 
     * The upgrade process guarantees a conforming mesh: nodes added on edges or faces 
     * shared by multiple cells are uniquely created and correctly referenced across 
     * neighboring elements.
     *
     * ### Operational Details
     * - **Memory Efficient:** Uses exact pre-allocation based on the current cell type 
     *   histogram to avoid any dynamic memory reallocation overhead.
     * - **Topology Rebuild:** Automatically rebuilds all connectivity structures 
     *   (inverse connectivities, facets, and cell adjacency).
     * - **Entity Update:** Refreshes all physical groups (MeshEntities) to include 
     *   the newly created higher-order vertices.
     * - **Observer Notification:** Triggers `updateAfterMeshChange()` on all attached 
     *   FE spaces so they can resize their internal DOF mappings.
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
     * @note If the mesh is already quadratic, this function returns immediately.
     * @warning This operation modifies the global mesh topology. It invalidates any 
     * previously acquired `std::span` views or pointers to connectivity arrays.
     *
     * @complexity
     * Approximately linear in the number of cells.
     ***************************************************************************************************************************************/
    void Mesh::upgradeToQuadratic() {
        if (this->isQuadratic()) return;
        PROFILE_SCOPE("Mesh::upgradeToQuadratic");
        
        CFEM_INFO << "Upgrading the mesh to quadratic... ";

        // FAST EDGE TRACKING
        struct Edge {
            CfemInt v1, v2;
            Edge(CfemInt a, CfemInt b) : v1(std::min(a, b)), v2(std::max(a, b)) {}
            
            bool operator==(const Edge& o) const { 
                return v1 == o.v1 && v2 == o.v2; 
            }
        };

        struct EdgeHash {
            size_t operator()(const Edge& e) const {
                if constexpr (sizeof(CfemInt) <= 4) { // If using 32 bits int
                    return (static_cast<uint64_t>(e.v1) << 32) | static_cast<uint64_t>(e.v2);
                } 
                else { // if using 64 bits int
                    size_t seed = 0;
                    seed ^= std::hash<CfemInt>{}(e.v1) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                    seed ^= std::hash<CfemInt>{}(e.v2) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                    return seed;
                }
            }
        };

        // FAST QUAD FACE TRACKING (Manual Sorting Network)
        struct QuadFace {
            CfemInt v[4];

            QuadFace(CfemInt a, CfemInt b, CfemInt c, CfemInt d) {
                v[0] = a; v[1] = b; v[2] = c; v[3] = d;
                
                // Avoid the STL overhead of std::sort call.
                if (v[0] > v[1]) std::swap(v[0], v[1]);
                if (v[2] > v[3]) std::swap(v[2], v[3]);
                if (v[0] > v[2]) std::swap(v[0], v[2]);
                if (v[1] > v[3]) std::swap(v[1], v[3]);
                if (v[1] > v[2]) std::swap(v[1], v[2]);
            }

            bool operator==(const QuadFace& o) const {
                return v[0] == o.v[0] && v[1] == o.v[1] && v[2] == o.v[2] && v[3] == o.v[3];
            }
        };

        struct QuadFaceHash {
            size_t operator()(const QuadFace& f) const {
                size_t seed = 0;
                // Boucle déroulée manuellement pour mixer les 4 noeuds déjà triés.
                seed ^= std::hash<CfemInt>{}(f.v[0]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                seed ^= std::hash<CfemInt>{}(f.v[1]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                seed ^= std::hash<CfemInt>{}(f.v[2]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                seed ^= std::hash<CfemInt>{}(f.v[3]) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
                return seed;
            }
        };

        std::unordered_map<Edge, CfemInt, EdgeHash> edge_to_midpoint;
        std::unordered_map<QuadFace, CfemInt, QuadFaceHash> quadface_to_midpoint;
    
        // Prepare the new Big Block for connectivity
        std::vector<CfemInt> newAllCellVertexIds;

        size_t totalNodes = 0;
        size_t totalEdges = 0;
        size_t totalQuadFaces = 0;

        for (size_t t = 0; t < static_cast<size_t>(CellType::COUNT); ++t) {
            size_t count = m_cellCountPerType[t];
            if (count == 0) continue;

            CellType type = static_cast<CellType>(t);

            totalNodes += count * Cell::getNodeCount(Cell::getQuadraticEquivalent(type));
            totalEdges += count * Cell::getEdgeCount(type);
            totalQuadFaces += count * Cell::getLocalQuadFacesCount(type);
        }

        newAllCellVertexIds.reserve(totalNodes);
        edge_to_midpoint.reserve(totalEdges);
        quadface_to_midpoint.reserve(totalQuadFaces);

        m_cellOffsets.clear();
        m_cellOffsets.reserve(m_cells.size() + 1);
        m_cellOffsets.push_back(0);

        std::vector<CfemInt> quad_buffer; 
        quad_buffer.reserve(27);

        // MAIN UPGRADE LOOP
        for (size_t i = 0; i < m_cells.size(); ++i) {
            Cell &cell = m_cells[i];
            auto linear_vids = getCellConnectivity(cell.id);

            quad_buffer.clear(); // Make sure the buffer is empty at each iteration loop

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

        // UPDATE FLAGS & FIX STD::UNIQUE BUG
        m_cellTypeFlags.reset();
        for (auto &type : m_uniqueCellTypes) {
            type = Cell::getQuadraticEquivalent(type);
            CfemInt typeIdx = static_cast<int>(type);
            if (!m_cellTypeFlags.test(typeIdx)) {
                m_cellTypeFlags.set(typeIdx);
            }
        }

        m_facetTypeFlags.reset();
        for (auto &type : m_uniqueFacetTypes) {
            type = Cell::getQuadraticEquivalent(type);
            CfemInt typeIdx = static_cast<int>(type);
            if (!m_facetTypeFlags.test(typeIdx)) {
                m_facetTypeFlags.set(typeIdx);
            }
        }

        std::array<size_t, static_cast<size_t>(CellType::COUNT)> newCellCountPerType{};
        newCellCountPerType.fill(0);

        for (size_t i = 0; i < m_cellCountPerType.size(); ++i)
        {
            const auto count = m_cellCountPerType[i];
            if (count == 0) continue;

            const auto qt = static_cast<size_t>( Cell::getQuadraticEquivalent(static_cast<CellType>(i)) );
            newCellCountPerType[qt] += count;
        }

        m_cellCountPerType = newCellCountPerType;

        std::sort(m_uniqueCellTypes.begin(), m_uniqueCellTypes.end());
        auto itCell = std::unique(m_uniqueCellTypes.begin(), m_uniqueCellTypes.end());
        m_uniqueCellTypes.erase(itCell, m_uniqueCellTypes.end());

        std::sort(m_uniqueFacetTypes.begin(), m_uniqueFacetTypes.end());
        auto itFacet = std::unique(m_uniqueFacetTypes.begin(), m_uniqueFacetTypes.end());
        m_uniqueFacetTypes.erase(itFacet, m_uniqueFacetTypes.end());

        // FINALIZE (Atomic Switch)
        m_allCellVertexIds = std::move(newAllCellVertexIds);
        m_meshOrder = MeshOrder::Quadratic;

        buildInverseConnectivities();
        buildFacets(); 
        buildAdjacency();
        updateEntities();

        notifyObservers();

#ifndef CFEM_TESTING

        CFEM_INFO << "Mesh successfully upgraded to quadratic.";
        CFEM_INFO << this->info();
#endif
    }

    /** ************************************************************************************************************************************
     * @brief Downgrades the mesh from quadratic (P2) to linear (P1) order.
     *
     * @details
     * Converts all quadratic elements into their linear counterparts by discarding 
     * higher-order nodes (mid-edge, face-center, and volume-center nodes) from the 
     * cell connectivity.
     * 
     * The downgrade process safely preserves the core geometric shape while executing 
     * a complete garbage collection of orphaned higher-order nodes.
     *
     * ### Operational Details
     * - **Exact Memory Allocation:** Calculates the exact required capacity for the 
     *   new connectivity block to guarantee zero dynamic reallocations.
     * - **Topology Compaction & Rebuild:** Automatically compacts the vertex array 
     *   to remove unused nodes, then strictly rebuilds all adjacency structures 
     *   (inverse connectivities, facets, and neighbors).
     * - **Entity Update:** Refreshes all physical groups (MeshEntities) to reflect 
     *   the new, compacted vertex set.
     * - **Observer Notification:** Triggers `updateAfterMeshChange()` on all attached 
     *   FE spaces so they can resize and remap their internal DOFs.
     *
     * ### Supported Element Transformations
     *
     * | Quadratic Type | Linear Type |
     * |----------------|-------------|
     * | LINE3          | LINE2       |
     * | TRI6           | TRI3        |
     * | QUAD9          | QUAD4       |
     * | TET10          | TET4        |
     * | PRISM18        | PRISM6      |
     * | HEX27          | HEX8        |
     *
     * @note If the mesh is already linear, this function returns immediately.
     * @warning This operation is highly destructive to the global mesh topology. 
     * Because it triggers `compactVertices()`, **existing Vertex IDs may change**. 
     * It also invalidates any previously acquired `std::span` views or pointers to 
     * connectivity arrays.
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

        // Prepare the new Big Bloc for connectivities
        std::vector<CfemInt> newAllCellVertexIds;
        size_t totalNodes = 0;

        for (size_t t = 0; t < static_cast<size_t>(CellType::COUNT); ++t) {
            size_t count = m_cellCountPerType[t];
            if (count == 0) continue;

            CellType type = static_cast<CellType>(t);
            totalNodes += count * Cell::getNodeCount(Cell::getLinearEquivalent(type));
        }
        newAllCellVertexIds.reserve(totalNodes);

        m_cellOffsets.clear();
        m_cellOffsets.reserve(m_cells.size() + 1);
        m_cellOffsets.push_back(0);

        // Loop over cells to truncate quadratic nodes
        for (auto &cell : m_cells)
        {
            // Get the current quadratic cell connectivity
            auto quad_nodes = getCellConnectivity(cell.id);
            
            CellType linearType = Cell::getLinearEquivalent(cell.type);
            CfemInt linearCount = Cell::getNodeCount(linearType);

            cell.offset = newAllCellVertexIds.size();
            cell.node_count = linearCount;
            cell.type = linearType;
            
            newAllCellVertexIds.insert( 
                newAllCellVertexIds.end(),
                quad_nodes.begin(),
                quad_nodes.begin() + linearCount
            );

            // CSR offset update for the next cell
            m_cellOffsets.push_back(newAllCellVertexIds.size());
        }

        // Update cell count per type
        std::array<size_t, static_cast<size_t>(CellType::COUNT)> newCellCountPerType{};
        newCellCountPerType.fill(0);

        for (size_t i = 0; i < m_cellCountPerType.size(); ++i)
        {
            const auto count = m_cellCountPerType[i];
            if (count == 0) continue;

            const auto qt = static_cast<size_t>( Cell::getLinearEquivalent(static_cast<CellType>(i)) );
            newCellCountPerType[qt] += count;
        }
        m_cellCountPerType = newCellCountPerType;

        // Update types flags
        m_cellTypeFlags.reset();
        for (auto &type : m_uniqueCellTypes)
        {
            type = Cell::getLinearEquivalent(type);
            m_cellTypeFlags.set(static_cast<int>(type));
        }
        // Clean up after upgrading
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

        // Data structure finalization
        m_allCellVertexIds = std::move(newAllCellVertexIds);
        m_meshOrder = MeshOrder::Linear;

        compactVertices(); 
        buildInverseConnectivities();
        buildFacets(); 
        buildAdjacency();
        
        updateEntities();
        
        m_isGloballyAffine = true;
        m_isCellAffine.assign(m_cells.size(), true);

        notifyObservers();
    }


    /** ************************************************************************************************************************************
     * @brief Synchronizes all physical groups (entities) with the current mesh topology.
     *
     * @details
     * Refreshes the internal geometric sets (vertices, facets, cells) of all registered 
     * mesh entities, including the global domain ("Universe"), the geometric boundary 
     * ("SYS_BOUNDARY"), and user-defined physical groups.
     * 
     * This must be invoked after any operation that modifies the underlying connectivity 
     * (e.g., order upgrades, downgrades, or manual vertex compaction) to ensure that 
     * material assignments and boundary conditions remain mapped to valid elements.
     *
     * ### Operational Details
     * - **Entity Preservation:** Entity IDs and names remain strictly intact; only their 
     *   underlying element lists are recomputed.
     * - **System-Aware:** Automatically triggers specialized rebuild routines for system 
     *   entities (Universe and Boundary) before updating standard user entities.
     *
     * @warning Failing to call this after a topological modification will result in 
     * physical groups referencing orphaned, invalid, or mismatched element IDs, leading 
     * to catastrophic failures during finite element assembly.
     *
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

    /** ************************************************************************************************************************************
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
     * @c upgradeToQuadratic() or mesh smoothing ( @c relaxMesh() ), to ensure the integration 
     * kernels use the correct mapping logic.
     *
     * @param tol Numerical tolerance for comparing node positions (default: 1e-10).
     ************************************************************************************************************************************ */
    void Mesh::updateAffineFlags() {
        const CfemInt numCells = getNumCells();
        
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

            // Check if cell Edges are affine
            const auto& edgeDefs = Cell::getLocalEdgesIndices(cell.type);
            for (const auto& edge : edgeDefs) {
                // For a quadratic cell, an edge has 3 vertices : {p1, p2, pm}
                if (edge.size() < 3) continue; 
                
                const Vector3D& p1 = coords[edge[0]];
                const Vector3D& p2 = coords[edge[1]];
                const Vector3D& pm = coords[edge[2]]; // Mid-node

                Vector3D theoreticalMid = (p1 + p2) * 0.5;
                
                if ((pm - theoreticalMid).norm() > tol) {
                    cellIsAffine = false;
                    break; // Unecessary to check the rest... the cell is curved.
                }
            }

            // Check if Quad Faces are affine (At this stage, triangular faces are affine: no mid-node)
            if (cellIsAffine) {
                const auto& quadFaceDefs = Cell::getLocalQuadFacesIndices(cell.type);
                for (const auto& face : quadFaceDefs) {
                    // Quad9, Prism18 and Hex27 have 9 nodes quad faces.
                    // The mid-node index is the last in the facet connectivity (node at index 8)
                    if (face.size() < 9) continue; 
                    
                    const Vector3D& p0 = coords[face[0]];
                    const Vector3D& p1 = coords[face[1]];
                    const Vector3D& p2 = coords[face[2]];
                    const Vector3D& p3 = coords[face[3]];
                    const Vector3D& pCenter = coords[face[8]]; // mid-node

                    Vector3D theoreticalCenter = (p0 + p1 + p2 + p3) * 0.25;
                    
                    if ((pCenter - theoreticalCenter).norm() > tol) {
                        cellIsAffine = false;
                        break;
                    }
                }
            }

            // ISOPARAMETRIC CENTER VERIFICATION (Hex27 only)
            if (cellIsAffine && cell.type == CellType::HEX27 && coords.size() == 27) {
                // For an affine transformation, the central node (node 26) must coincide 
                // exactly with the geometric barycenter of the 8 corner vertices.
                Vector3D theoreticalBarycenter = Vector3D(0, 0, 0);
                for (CfemInt i = 0; i < 8; ++i) {
                    theoreticalBarycenter += coords[i];
                }
                theoreticalBarycenter = theoreticalBarycenter * 0.125;

                if ((coords[26] - theoreticalBarycenter).norm() > tol) {
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
     * @brief Compute a quality metric for a specific facet.
     *
     * Quality is a dimensionless number between 0 and 1, where 1 represents
     * an ideal element (perfectly shaped) and 0 represents a degenerate element
     * (zero area/volume or inverted orientation).
     *
     * Supported types:
     * - TRI3, TRI6, QUAD4, QUAD9
     *
     * @param facetId The global index of the facet.
     * @return A real number in [0,1] measuring the element quality.
     ***************************************************************************************************************************************/
    CfemReal Mesh::computeFacetQuality(CfemInt facetId) const
    {
        const auto &facet = getFacet(facetId);

        switch (facet.type)
        {
        case CellType::TRI3:
        case CellType::TRI6:
            return computeTriQuality(getFacetConnectivity(facetId));

        case CellType::QUAD4:
        case CellType::QUAD9:
            return computeQuadQuality(getFacetConnectivity(facetId));

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

        CfemInt collapsed = 0;
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