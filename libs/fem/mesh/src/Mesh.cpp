#include "Mesh.h"
#include "FESpace.h"
#include <iostream>
#include <algorithm>
#include <map>
#include <iomanip>
#include <sstream>

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
        validateConnectivity(type, v_ids);
        CfemInt id = m_cells.size();

        m_cells.emplace_back();
        auto& cell = m_cells.back();
        cell.type = type;
        cell.vertex_ids = v_ids;
       // m_cells.push_back({type, v_ids});

        int typeIdx = static_cast<int>(type);
        if (!m_typeFlags.test(typeIdx))
        {
            m_typeFlags.set(typeIdx);
            m_uniqueTypes.push_back(type);
        }

        CFEM_POSTCONDITION(!m_cells.empty());
        return id;
    }

    // CfemInt Mesh::addTri3CellDirect(CfemInt v0, CfemInt v1, CfemInt v2) {
    //     // Emplace back to avoid copying the whole Cell struct
    //     CFEM_ASSERT()
    //     m_cells.emplace_back(); 
    //     Cell& cell = m_cells.back();
    //     cell.type = type;
        
    //     // Crucial: Use the size from getNodeCountPerCell to reserve
    //     // This prevents the internal vertex_ids from reallocating later
    //     cell.vertex_ids = ids; 
    // }
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
        std::vector<CfemInt> counts(numVertices, 0);
        CfemInt totalEntries = 0;
        for (const auto &cell : m_cells)
        {
            for (CfemInt v_id : cell.vertex_ids)
            {
                counts[v_id]++;
                totalEntries++;
            }
        }

        // Compute Offsets (Prefix Sum)
        // We store where each vertex's data starts in the global array
        std::vector<CfemInt> offsets(numVertices + 1, 0);
        for (CfemInt i = 0; i < numVertices; ++i)
        {
            offsets[i + 1] = offsets[i] + counts[i];
        }

        // Allocate the "Big Block"
        m_globalConnectivity.assign(totalEntries, -1);

        // Fill the Big Block
        // We reuse the 'counts' vector as a local tracker for the current insertion point
        std::fill(counts.begin(), counts.end(), 0);
        for (CfemInt c_id = 0; c_id < numCells; ++c_id)
        {
            for (CfemInt v_id : m_cells[c_id].vertex_ids)
            {
                CfemInt writePos = offsets[v_id] + counts[v_id];
                m_globalConnectivity[writePos] = c_id;
                counts[v_id]++;
            }
        }

        // 5. Update Vertex pointers
        // This allows your existing cell loops to work without knowing about CSR math
        for (CfemInt v_id = 0; v_id < numVertices; ++v_id)
        {
            m_vertices[v_id].connected_cells_ptr = &m_globalConnectivity[offsets[v_id]];
            m_vertices[v_id].num_connected_cells = offsets[v_id + 1] - offsets[v_id];
        }
    }

    /**
     * @brief Identifies unique facets (edges/faces) and sets up neighbor relationships.
     */
    void Mesh::buildFacets() {
        PROFILE_SCOPE("Mesh::buildFacets");
        m_facets.clear();
        m_facets.reserve(m_cells.size() * 2); // Heuristic

        // Mapping: cell_id -> local_face_index -> global_facet_id
        // Max faces is 6 (for HEX). Using a flat vector for cache efficiency.
        std::vector<CfemInt> cell_face_to_facet(m_cells.size() * 6, -1);

        for (CfemInt c_id = 0; c_id < static_cast<CfemInt>(m_cells.size()); ++c_id) {
            const int num_faces = getLocalFacesCount(m_cells[c_id].type);

            for (int f_idx = 0; f_idx < num_faces; ++f_idx) {
                // 1. Skip if a neighbor already defined this facet
                if (cell_face_to_facet[c_id * 6 + f_idx] != -1) continue;

                LocalFace current_face = getLocalFaceData(c_id, f_idx);
                CfemInt v_anchor = current_face.v[0];
                const auto& anchor_vert = m_vertices[v_anchor];

                CfemInt neighbor_id = -1;
                int neighbor_f_idx = -1;

                // 2. Search neighbors via CSR Inverse Connectivity
                for (CfemInt i = 0; i < anchor_vert.num_connected_cells; ++i) {
                    CfemInt other_c_id = anchor_vert.connected_cells_ptr[i];
                    if (other_c_id == c_id) continue;

                    // Check faces of the neighbor
                    const int other_num_faces = getLocalFacesCount(m_cells[other_c_id].type);
                    for (int j = 0; j < other_num_faces; ++j) {
                        LocalFace other_face = getLocalFaceData(other_c_id, j);
                        if (current_face.matches(other_face)) {
                            neighbor_id = other_c_id;
                            neighbor_f_idx = j;
                            break;
                        }
                    }
                    if (neighbor_id != -1) break;
                }

                // 3. Register the Facet
                Facet f;
                f.id = static_cast<CfemInt>(m_facets.size());
                f.type = current_face.type;
                f.vertex_ids.assign(current_face.v, current_face.v + current_face.size);
                f.left_cell = c_id;
                f.right_cell = neighbor_id;
                f.number_of_vertices = static_cast<CfemInt>(current_face.size);

                cell_face_to_facet[c_id * 6 + f_idx] = f.id;
                if (neighbor_id != -1) {
                    cell_face_to_facet[neighbor_id * 6 + neighbor_f_idx] = f.id;
                }

                m_facets.push_back(std::move(f));
            }
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

    MeshEntity &Mesh::createInternalEntity(CfemInt id, const std::string &name, CfemInt dim)
    {
        CFEM_PRECONDITION(m_entity_ids.find(id) == m_entity_ids.end());

        m_entity_ids.insert(id);
        m_entities.push_back(std::make_unique<MeshEntity>(this, id, name, dim));

        return *m_entities.back();
    }

    void Mesh::buildMesh()
    {
        if (m_cells.empty() || m_vertices.empty() || m_meshDim == 0)
        {
            CFEM_ERROR("buildMesh failed: Mesh is empty.");
        }

        buildInverseConnectivities();
        buildFacets();

        MeshEntity *universe = this->getEntity(cfem::UNIVERSE_ID);
        MeshEntity *boundary = this->getEntity(cfem::SYS_BOUNDARY_ID);

        // Total mesh universe entity
        if (!universe)
        {
            universe = &(createInternalEntity(cfem::UNIVERSE_ID, "Universe", m_meshDim));
        }
        if (!boundary)
        {
            boundary = &(createInternalEntity(cfem::SYS_BOUNDARY_ID, "SYS_BOUNDARY", m_meshDim - 1));
        }

        if (universe->getCellIds().empty())
        {
            universe->clear();
            for (CfemInt i = 0; i < static_cast<CfemInt>(m_vertices.size()); ++i)
            {
                universe->addVertex(i);
            }
            for (CfemInt i = 0; i < static_cast<CfemInt>(m_cells.size()); ++i)
            {
                universe->addCell(i);
            }
        }
        else if (universe->getVertexIds().empty())
        {
            for (CfemInt cell_id : universe->getVertexIds())
            {
                for (CfemInt v_id : getCell(cell_id).vertex_ids)
                {
                    universe->addVertex(v_id);
                }
            }
        }

        for (const auto &facet : m_facets)
        {
            universe->addFacet(facet.id);
            if (facet.right_cell == -1)
            {
                boundary->addFacet(facet.id);
                for (CfemInt i = 0; i < facet.number_of_vertices; ++i)
                {
                    boundary->addVertex(facet.vertex_ids[i]);
                }
            }
        }

        universe->finalize();
        boundary->finalize();

#ifndef NDEBUG
        std::cout << "[CFEM INFO] Mesh built successfully: "
                  << m_vertices.size() << " Nodes, "
                  << m_cells.size() << " Cells, "
                  << m_facets.size() << " Facets." << std::endl;
#endif
    }

    std::string Mesh::info() const
    {
        std::stringstream ss;
        ss << "\n================ MESH SUMMARY ================" << std::endl;
        ss << " Dimension: " << getDimension() << "D" << std::endl;
        ss << " Order:     " << (getOrder() == MeshOrder::Linear ? "Lineae" : (getOrder() == MeshOrder::Quadratic) ? "Quadratic"
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
               << std::setw(10) << v.coords[0] << " "
               << std::setw(10) << v.coords[1] << " "
               << std::setw(10) << v.coords[2] << "\n";
        }
        ss << "$EndNodes\n";

        // --- Section 3: Cells ---
        ss << "$Cells\n";
        ss << m_cells.size() << " # total number of cells\n";
        for (size_t i = 0; i < m_cells.size(); ++i)
        {
            const auto &cell = m_cells[i];

            // Output: ID Type NodeCount [Nodes...]
            ss << i << " "
               << std::setw(8) << cellTypeToString(cell.type) << " "
               << cell.vertex_ids.size() << "  ";

            for (auto v_id : cell.vertex_ids)
            {
                ss << v_id << " ";
            }
            ss << "\n";
        }
        ss << "$EndCells\n";

        // --- Section 4: Facets ---
        ss << "$Facets\n";
        ss << m_facets.size() << " # total number of facets\n";
        for (const auto &f : m_facets)
        {
            // ID | NodeCount | Nodes... | LeftCell | RightCell
            ss << f.id << " " << f.number_of_vertices << "  ";

            for (CfemInt i = 0; i < f.number_of_vertices; ++i)
            {
                ss << f.vertex_ids[i] << " ";
            }

            // Use -1 (or another marker) to clearly show Boundary vs Internal
            ss << "| L: " << std::setw(5) << f.left_cell
               << " R: " << std::setw(5) << f.right_cell << "\n";
        }
        ss << "$EndFacets\n";

        return ss.str();
    }

    void Mesh::upgradeToQuadratic()
    {
        if (this->isQuadratic())
            return;

        PROFILE_SCOPE("Mesh::upgradeToQuadratic");

        // 1. Map to track unique edges: {min_id, max_id} -> new_vertex_id
        std::map<std::pair<CfemInt, CfemInt>, CfemInt> edge_to_midpoint;

        for (auto &cell : m_cells)
        {
            //std::vector<CfemInt> v = cell.vertex_ids; // Copy original linear nodes

            if (cell.type == CellType::TRI3)
            {
                upgradeTri3ToTri6(cell, edge_to_midpoint);
            }
            else if (cell.type == CellType::QUAD4)
            {
                upgradeQuad4ToQuad9(cell, edge_to_midpoint);
            }
            else if (cell.type == CellType::TET4)
            {
                // Edges: 01, 12, 20, 03, 13, 23
                upgradeTetra4ToTetra10(cell, edge_to_midpoint);
            }
            else if (cell.type == CellType::PRISM6)
            {
                upgradePrism6ToPrism18(cell, edge_to_midpoint);
            }
            else if (cell.type == CellType::HEX8)
            {
                upgradeHex8ToHex27(cell, edge_to_midpoint);
            }
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

        m_meshOrder = MeshOrder::Quadratic;

        buildInverseConnectivities();

        buildFacets();
        // Rebuild facets to include mid-nodes on boundaries

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

        for (auto &cell : m_cells)
        {
            if (cell.type == CellType::LINE3)
            {
                cell.type = CellType::LINE2;
                cell.vertex_ids.resize(2); // Keep only corners
            }
            else if (cell.type == CellType::TRI6)
            {
                cell.type = CellType::TRI3;
                cell.vertex_ids.resize(3);
            }
            else if (cell.type == CellType::QUAD9)
            {
                cell.type = CellType::QUAD4;
                cell.vertex_ids.resize(4);
            }
            else if (cell.type == CellType::TET10)
            {
                cell.type = CellType::TET4;
                cell.vertex_ids.resize(4);
            }
            else if (cell.type == CellType::PRISM18)
            {
                cell.type = CellType::PRISM18;
                cell.vertex_ids.resize(6);
            }
            else if (cell.type == CellType::HEX27)
            {
                cell.type = CellType::HEX8;
                cell.vertex_ids.resize(8);
            }
        }

        // Optionally: I could 'vacuum' here to remove
        // orphaned vertices, but that requires re-indexing everything.
        // Usually, it's enough to just leave them as isolated points.

        m_meshOrder = MeshOrder::Linear;
        buildInverseConnectivities();
        buildFacets(); // Rebuild facets to remove mid-nodes
        updateEntities();
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

    CfemReal Mesh::computeCellMeasure(CellType type, const std::vector<CfemInt> &v_ids) const
    {
        switch (type)
        {
        case CellType::LINE2:
            return compute2NodesLineLength(v_ids[0], v_ids[1]);

        case CellType::TRI3:
        {
            cfem::Vector3D lReferenceNormale{0.0, 0.0, 0.0};
            return computeSigned3NodesTriangleArea(v_ids[0], v_ids[1], v_ids[2], lReferenceNormale);
        }

        case CellType::QUAD4:
        {
            cfem::Vector3D lReferenceNormale{0.0, 0.0, 0.0};
            CfemReal A1 = computeSigned3NodesTriangleArea(v_ids[0], v_ids[1], v_ids[2], lReferenceNormale);
            CfemReal A2 = computeSigned3NodesTriangleArea(v_ids[0], v_ids[2], v_ids[3], lReferenceNormale);

            if ((std::abs(A1) <= cfem::ZERO_TOLERANCE) || (std::abs(A2) <= cfem::ZERO_TOLERANCE))
            {
                CFEM_WARNING("Sub-triangle in QUAD4 has zero area. Check if nodes are colinear.");
            }

            if ((A1 > cfem::ZERO_TOLERANCE && A2 < -cfem::ZERO_TOLERANCE) ||
                (A1 < -cfem::ZERO_TOLERANCE && A2 > cfem::ZERO_TOLERANCE))
            {
                CFEM_WARNING("Poor quality QUAD4 cell: Self-intersecting geometry detected.");
            }

            CfemReal totalArea = std::abs(A1 + A2);
            if (totalArea < cfem::ZERO_TOLERANCE)
            {
                CFEM_WARNING("Poor quality QUAD4 cell: Inverted or degenerate total area.");
            }
            return totalArea;
        }

        case CellType::TET4:
        case CellType::TET10:
            return computeSigned4NodesTetraVolume(v_ids[0], v_ids[1], v_ids[2], v_ids[3]);

        case CellType::PRISM6:
        {
            CfemReal lTotalVolume{0.0};
            CfemInt positiveTets = 0;
            CfemInt negativeTets = 0;

            std::vector<std::vector<int>> subsets = {
                {v_ids[0], v_ids[1], v_ids[2], v_ids[3]},
                {v_ids[1], v_ids[2], v_ids[3], v_ids[4]},
                {v_ids[2], v_ids[3], v_ids[4], v_ids[5]}};

            for (const auto &tetra_v_ids : subsets)
            {
                CfemReal tetra_vol = computeSigned4NodesTetraVolume(tetra_v_ids[0], tetra_v_ids[1], tetra_v_ids[2], tetra_v_ids[3]);

                if (std::abs(tetra_vol) <= cfem::ZERO_TOLERANCE)
                {
                    CFEM_WARNING("Poor quality PRISM6 cell: Sub-Tetra in PRISM6 has zero area. Check if nodes are coplanar.");
                }

                if (tetra_vol > 0)
                    positiveTets++;
                else if (tetra_vol < 0)
                    negativeTets++;

                lTotalVolume += tetra_vol;
            }

            if (positiveTets > 0 && negativeTets > 0)
            {
                CFEM_WARNING("Invalid PRISM6: Twisted/Self-intersecting geometry detected.");
            }

            return std::abs(lTotalVolume);
        }

        case CellType::HEX8:
        case CellType::HEX27:
        {
            CfemReal lTotalVolume{0.0};
            CfemInt positiveTets = 0;
            CfemInt negativeTets = 0;

            std::vector<std::vector<int>> subsets = {
                {v_ids[0], v_ids[1], v_ids[3], v_ids[4]},
                {v_ids[1], v_ids[3], v_ids[4], v_ids[5]},
                {v_ids[4], v_ids[5], v_ids[3], v_ids[7]},
                {v_ids[1], v_ids[3], v_ids[5], v_ids[2]},
                {v_ids[3], v_ids[7], v_ids[5], v_ids[2]},
                {v_ids[5], v_ids[7], v_ids[6], v_ids[2]}};

            for (const auto &tetra_v_ids : subsets)
            {
                CfemReal tetra_vol = computeSigned4NodesTetraVolume(tetra_v_ids[0], tetra_v_ids[1], tetra_v_ids[2], tetra_v_ids[3]);

                if (std::abs(tetra_vol) <= cfem::ZERO_TOLERANCE)
                {
                    CFEM_WARNING("Invalid HEX8: Sub-tetra in HEX8 has zero volume. Check if nodes are coplanar.");
                }

                if (tetra_vol > 0)
                    positiveTets++;
                else if (tetra_vol < 0)
                    negativeTets++;

                lTotalVolume += tetra_vol;
            }

            if (positiveTets > 0 && negativeTets > 0)
            {
                CFEM_WARNING("Invalid HEX8: Twisted/Self-intersecting geometry detected.");
            }

            return std::abs(lTotalVolume);
        }

        default:
            CFEM_ERROR("Measure calculation not implemented for " + cellTypeToString(type));
        }
    }

    CfemReal Mesh::compute2NodesLineLength(CfemInt p_idx0, CfemInt p_idx1) const
    {
        const auto &x0 = m_vertices[p_idx0].coords;
        const auto &x1 = m_vertices[p_idx1].coords;

        const CfemReal dx = x1[0] - x0[0];
        const CfemReal dy = x1[1] - x0[1];
        const CfemReal dz = x1[2] - x0[2];

        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    CfemReal Mesh::computeSigned3NodesTriangleArea(CfemInt p_idx0,
                                                   CfemInt p_idx1,
                                                   CfemInt p_idx2,
                                                   cfem::Vector3D &pReferenceNormale) const
    {
        const auto &A = m_vertices[p_idx0].coords;
        const auto &B = m_vertices[p_idx1].coords;
        const auto &C = m_vertices[p_idx2].coords;

        cfem::Vector3D AB = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
        cfem::Vector3D AC = {C[0] - A[0], C[1] - A[1], C[2] - A[2]};

        CfemReal cp_x = AB[1] * AC[2] - AB[2] * AC[1];
        CfemReal cp_y = AB[2] * AC[0] - AB[0] * AC[2];
        CfemReal cp_z = AB[0] * AC[1] - AB[1] * AC[0];

        CfemReal double_area_mag = std::sqrt(cp_x * cp_x + cp_y * cp_y + cp_z * cp_z);

        if (double_area_mag < cfem::ZERO_TOLERANCE)
        {
            return 0.0;
        }

        bool is_ref_zero = (std::abs(pReferenceNormale[0]) < cfem::ZERO_TOLERANCE &&
                            std::abs(pReferenceNormale[1]) < cfem::ZERO_TOLERANCE &&
                            std::abs(pReferenceNormale[2]) < cfem::ZERO_TOLERANCE);

        if (is_ref_zero)
        {
            pReferenceNormale[0] = cp_x / double_area_mag;
            pReferenceNormale[1] = cp_y / double_area_mag;
            pReferenceNormale[2] = cp_z / double_area_mag;

            return 0.5 * double_area_mag;
        }

        // Project the current normal onto the Reference Normal (Dot Product)
        // This provides the signed magnitude relative to the bootstrap orientation
        CfemReal signed_double_area = (pReferenceNormale[0] * cp_x +
                                       pReferenceNormale[1] * cp_y +
                                       pReferenceNormale[2] * cp_z);

        return 0.5 * signed_double_area;
    }

    // Internal helper for measure calculation
    CfemReal Mesh::computeSigned4NodesTetraVolume(CfemInt p_idx0, CfemInt p_idx1, CfemInt p_idx2, CfemInt p_idx3) const
    {
        const auto &p0 = m_vertices[p_idx0].coords;
        const auto &p1 = m_vertices[p_idx1].coords;
        const auto &p2 = m_vertices[p_idx2].coords;
        const auto &p3 = m_vertices[p_idx3].coords;

        // Standard Scalar Triple Product (1/6)
        // (p1-p0) . ((p2-p0) x (p3-p0))
        return (
                   (p1[0] - p0[0]) * ((p2[1] - p0[1]) * (p3[2] - p0[2]) - (p2[2] - p0[2]) * (p3[1] - p0[1])) +
                   (p1[1] - p0[1]) * ((p2[2] - p0[2]) * (p3[0] - p0[0]) - (p2[0] - p0[0]) * (p3[2] - p0[2])) +
                   (p1[2] - p0[2]) * ((p2[0] - p0[0]) * (p3[1] - p0[1]) - (p2[1] - p0[1]) * (p3[0] - p0[0]))) /
               6.0;
    }

    CfemReal Mesh::computeCellQuality(CfemInt cellId) const
    {
        const auto &cell = getCell(cellId);

        switch (cell.type)
        {
        case CellType::TRI3:
        case CellType::TRI6:
            return computeTriQuality(cell.vertex_ids);

        case CellType::QUAD4:
        case CellType::QUAD9:
            return computeQuadQuality(cell.vertex_ids);

        case CellType::TET4:
        case CellType::TET10:
            return computeTetQuality(cell.vertex_ids);

        default:
            return 1.0; // Fallback for types not yet implemented
        }
    }

    CfemReal Mesh::computeQuadQuality(const std::vector<CfemInt> &v) const
    {
        CfemReal q1 = computeTriQuality({v[0], v[1], v[2]});
        CfemReal q2 = computeTriQuality({v[0], v[2], v[3]});
        return (q1 + q2) / 2.0;
    }

    CfemReal Mesh::computeTriQuality(const std::vector<CfemInt> &v) const
    {
        const cfem::Vector3D &p0 = getVertex(v[0]).coords;
        const cfem::Vector3D &p1 = getVertex(v[1]).coords;
        const cfem::Vector3D &p2 = getVertex(v[2]).coords;

        // Calculate edge lengths
        CfemReal a = (p1 - p2).norm();
        CfemReal b = (p0 - p2).norm();
        CfemReal c = (p0 - p1).norm();

        // Semi-perimeter
        CfemReal s = (a + b + c) / 2.0;

        // Area (Heron's Formula)
        CfemReal area = std::sqrt(s * (s - a) * (s - b) * (s - c));

        // Radius Ratio = 2 * (Inradius / Circumradius)
        // Simplified formula: (b+c-a)*(c+a-b)*(a+b-c) / (a*b*c)
        if (a * b * c == 0)
            return 0.0;
        CfemReal quality = ((b + c - a) * (c + a - b) * (a + b - c)) / (a * b * c);

        return std::max(0.0, std::min(1.0, quality));
    }

    CfemReal Mesh::computeTetQuality(const std::vector<CfemInt> &v) const
    {
        const cfem::Vector3D &p0 = getVertex(v[0]).coords;
        const cfem::Vector3D &p1 = getVertex(v[1]).coords;
        const cfem::Vector3D &p2 = getVertex(v[2]).coords;
        const cfem::Vector3D &p3 = getVertex(v[3]).coords;

        // Volume of Tetrahedron
        CfemReal vol = std::abs((p1 - p0).dot((p2 - p0).cross(p3 - p0))) / 6.0;
        if (vol < 1e-15)
            return 0.0;

        // Sum of face areas (simplified approximation for quality)
        CfemReal a0 = (p1 - p0).cross(p2 - p0).norm() * 0.5;
        CfemReal a1 = (p1 - p0).cross(p3 - p0).norm() * 0.5;
        CfemReal a2 = (p2 - p0).cross(p3 - p0).norm() * 0.5;
        CfemReal a3 = (p2 - p1).cross(p3 - p1).norm() * 0.5;

        CfemReal rms_edge = ((p1 - p0).norm() + (p2 - p0).norm() + (p3 - p0).norm() +
                             (p2 - p1).norm() + (p3 - p1).norm() + (p3 - p2).norm()) /
                            6.0;

        // Scaled Ratio: 12 * (3 * Vol / AreaSum) / rms_edge
        // This is a standard normalized metric where 1.0 is ideal.
        return std::min(1.0, (124.7 * vol) / std::pow(rms_edge, 3));
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

    void Mesh::validateConnectivity(CellType type, const std::vector<CfemInt> &v_ids) const
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

#ifndef NDEBUG
        if (computeCellMeasure(type, v_ids) <= cfem::ZERO_TOLERANCE)
        {
            CFEM_WARNING("Zero or negative volume element detected! Mesh is likely invalid.");
        }
#endif
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

} // namespace