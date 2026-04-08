#ifndef CFEM_MESH_TPP
#define CFEM_MESH_TPP
#include "Mesh.h"
#include "DynamicVector.h"
#endif

namespace cfem
{

    /** ************************************************************************************************************************************
     * @brief Retrieves the 3D coordinates of all vertices of a given cell.
     *
     * This function populates a user-provided dynamic vector with the coordinates
     * of each vertex that forms the specified cell. The order of coordinates
     * corresponds to the cell's connectivity.
     *
     * This is useful for geometric computations such as:
     * - Calculating cell volumes or areas
     * - Computing face normals
     * - Visualization and mesh export
     *
     * @param cell_id The ID of the cell whose vertices are to be retrieved. Must satisfy 0 <= cell_id < getNumCells().
     * @param v_coords Reference to a DynamicVector<Vector3D> that will be resized to the number of vertices in the cell and filled with their
     *                  coordinates.
     * @pre cell_id must be a valid index in the mesh.
     * @post v_coords.size() == number of vertices in the cell.
     *
     * @note The function works for any cell type (linear or quadratic) since it relies
     *       on the cell's connectivity array.
     ***************************************************************************************************************************************/
    inline void Mesh::getCellVerticesCoords(CfemInt cell_id, DynamicVector<Vector3D> &v_coords) const
    {
        CFEM_ASSERT(cell_id >= 0 && cell_id < getNumCells());
        auto node_ids = getCellConnectivity(cell_id);
        const size_t numVertices = node_ids.size();
        v_coords.resize(numVertices);

        for (size_t i = 0; i < numVertices; ++i)
        {
            v_coords[i] = m_vertices[node_ids[i]].coords;
        }
    }

    /** ************************************************************************************************************************************
     * @brief Retrieves the 3D coordinates of all vertices of a given facet.
     *
     * This function populates a user-provided dynamic vector with the coordinates
     * of each vertex that forms the specified facet. The order of coordinates
     * corresponds to the facet's connectivity.
     *
     * This is useful for geometric computations such as:
     * - Calculating facet lenght areas
     * - Visualization and mesh export
     *
     * @param facet_id The ID of the facet whose vertices are to be retrieved. Must satisfy 0 <= facet_id < getNumFacets().
     * @param v_coords Reference to a DynamicVector<Vector3D> that will be resized to the number of vertices in the facet and filled with their
     *                  coordinates.
     * @pre facet_id must be a valid index in the mesh.
     * @post v_coords.size() == number of vertices in the facet.
     *
     * @note The function works for any facet type (linear or quadratic) since it relies
     *       on the facet's connectivity array.
     ***************************************************************************************************************************************/
    inline void Mesh::getFacetVerticesCoords(CfemInt cell_id, DynamicVector<Vector3D> &v_coords) const
    {
        CFEM_ASSERT(cell_id >= 0 && cell_id < getNumFacets());
        auto node_ids = getFacetConnectivity(cell_id);
        const size_t numVertices = node_ids.size();
        v_coords.resize(numVertices);

        for (size_t i = 0; i < numVertices; ++i)
        {
            v_coords[i] = m_vertices[node_ids[i]].coords;
        }
    }

    inline CfemInt getNodeCountPerCell(CellType type)
    {
        switch (type)
        {
        case CellType::LINE2:   return 2;
        case CellType::TRI3:    return 3;
        case CellType::QUAD4:   return 4;
        case CellType::TET4:    return 4;
        case CellType::PRISM6:  return 6;
        case CellType::HEX8:    return 8;
        case CellType::LINE3:   return 3;
        case CellType::TRI6:    return 6;
        case CellType::QUAD9:   return 9;
        case CellType::TET10:   return 10;
        case CellType::PRISM18: return 18;
        case CellType::HEX27:   return 27;
        default:
            return 0;
        }
    }


    /** ************************************************************************************************************************************
     * @brief Retrieve the local face (facet) information of a cell.
     *
     * This function returns the connectivity and type of a specific face of a given cell. The face is identified by its **local face index**
     * within the cell.
     *
     * The returned structure contains:
     * - the face element type (e.g. LINE2, TRI3, QUAD4, TRI6, QUAD9)
     * - the number of vertices defining the face
     * - the vertex indices corresponding to the mesh vertex array
     *
     * The vertex ordering follows the internal CFEM convention and is consistent across cells to ensure correct facet matching when
     * building mesh adjacency.
     *
     * This routine supports both **linear** and **quadratic** elements. For higher-order cells, the returned face also includes mid-edge
     * or mid-face nodes when required.
     *
     * ---
     *
     * ### Example
     *
     * For a tetrahedral cell (TET4):
     *
     * | face_idx | Face vertices |
     * |----------|---------------|
     * | 0        | (v0, v1, v2)  |
     * | 1        | (v0, v1, v3)  |
     * | 2        | (v1, v2, v3)  |
     * | 3        | (v2, v0, v3)  |
     *
     * For a hexahedral cell (HEX8):
     *
     * | face_idx | Face vertices |
     * |----------|---------------|
     * | 0        | bottom face   |
     * | 1        | top face      |
     * | 2        | front face    |
     * | 3        | right face    |
     * | 4        | back face     |
     * | 5        | left face     |
     *
     * ---
     *
     * @param cell_id Index of the cell in the mesh.
     * @param face_idx Local index of the face within the cell.
     *
     * @return A LocalFace structure containing:
     *         - the face type
     *         - number of vertices
     *         - vertex indices defining the face
     *
     * @pre
     * - `cell_id` must be a valid cell index.
     * - `face_idx` must be within the number of faces for the cell type.
     *
     * @note
     * This function is heavily used internally when building mesh topology structures such as:
     * - facet lists
     * - cell adjacency
     * - boundary detection
     *
     * @see getLocalFacesCount()
     * @see buildFacets()
     * @see buildAdjacency()
     ***************************************************************************************************************************************/
    inline LocalFace Mesh::getLocalFaceData(CfemInt cell_id, int face_idx) const
    {
        CFEM_ASSERT(cell_id >= 0 && cell_id < getNumCells());
        const auto &cell = m_cells[cell_id];
        const auto &v = getCellConnectivity(cell_id);
        LocalFace lf;

        switch (cell.type)
        {
        case CellType::LINE2:
        case CellType::LINE3:
            lf.type = CellType::NODE;
            lf.size = 1;
            lf.v[0] = v[face_idx]; // face_idx 0 or 1
            break;

        case CellType::TRI3:
            lf.type = CellType::LINE2;
            lf.size = 2;
            if (face_idx == 0)
            {
                lf.v[0] = v[0]; lf.v[1] = v[1];
            }
            else if (face_idx == 1)
            {
                lf.v[0] = v[1]; lf.v[1] = v[2];
            }
            else
            {
                lf.v[0] = v[2]; lf.v[1] = v[0];
            }
            break;

        case CellType::QUAD4:
            lf.type = CellType::LINE2;
            lf.size = 2;
            lf.v[0] = v[face_idx];
            lf.v[1] = v[(face_idx + 1) % 4];
            break;

        case CellType::TET4:
            lf.type = CellType::TRI3;
            lf.size = 3;
            if (face_idx == 0)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[1];
                lf.v[2] = v[2];
            }
            else if (face_idx == 1)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[1];
                lf.v[2] = v[3];
            }
            else if (face_idx == 2)
            {
                lf.v[0] = v[1];
                lf.v[1] = v[2];
                lf.v[2] = v[3];
            }
            else
            {
                lf.v[0] = v[2];
                lf.v[1] = v[0];
                lf.v[2] = v[3];
            }
            break;

        case CellType::PRISM6:
            if (face_idx < 2)
            {
                lf.type = CellType::TRI3;
                lf.size = 3;
            }
            else
            {
                lf.type = CellType::QUAD4;
                lf.size = 4;
            }
            if (face_idx == 0)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[2];
                lf.v[2] = v[1];
            }
            else if (face_idx == 1)
            {
                lf.v[0] = v[3];
                lf.v[1] = v[4];
                lf.v[2] = v[5];
            }
            else if (face_idx == 2)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[1];
                lf.v[2] = v[4];
                lf.v[3] = v[3];
            }
            else if (face_idx == 3)
            {
                lf.v[0] = v[1];
                lf.v[1] = v[2];
                lf.v[2] = v[5];
                lf.v[3] = v[4];
            }
            else
            {
                lf.v[0] = v[2];
                lf.v[1] = v[0];
                lf.v[2] = v[3];
                lf.v[3] = v[5];
            }
            break;

        case CellType::HEX8:
            lf.type = CellType::QUAD4;
            lf.size = 4;
            if (face_idx == 0)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[3];
                lf.v[2] = v[2];
                lf.v[3] = v[1];
            }
            else if (face_idx == 1)
            {
                lf.v[0] = v[4];
                lf.v[1] = v[5];
                lf.v[2] = v[6];
                lf.v[3] = v[7];
            }
            else if (face_idx == 2)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[1];
                lf.v[2] = v[5];
                lf.v[3] = v[4];
            }
            else if (face_idx == 3)
            {
                lf.v[0] = v[1];
                lf.v[1] = v[2];
                lf.v[2] = v[6];
                lf.v[3] = v[5];
            }
            else if (face_idx == 4)
            {
                lf.v[0] = v[2];
                lf.v[1] = v[3];
                lf.v[2] = v[7];
                lf.v[3] = v[6];
            }
            else
            {
                lf.v[0] = v[3];
                lf.v[1] = v[0];
                lf.v[2] = v[4];
                lf.v[3] = v[7];
            }
            break;

        case CellType::TRI6:
            lf.type = CellType::LINE3;
            lf.size = 3;
            if (face_idx == 0)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[1];
                lf.v[2] = v[3];
            }
            else if (face_idx == 1)
            {
                lf.v[0] = v[1];
                lf.v[1] = v[2];
                lf.v[2] = v[4];
            }
            else
            {
                lf.v[0] = v[2];
                lf.v[1] = v[0];
                lf.v[2] = v[5];
            }
            break;

        case CellType::QUAD9:
            lf.type = CellType::LINE3;
            lf.size = 3;
            if (face_idx == 0)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[1];
                lf.v[2] = v[4];
            }
            else if (face_idx == 1)
            {
                lf.v[0] = v[1];
                lf.v[1] = v[2];
                lf.v[2] = v[5];
            }
            else if (face_idx == 2)
            {
                lf.v[0] = v[2];
                lf.v[1] = v[3];
                lf.v[2] = v[6];
            }
            else
            {
                lf.v[0] = v[3];
                lf.v[1] = v[0];
                lf.v[2] = v[7];
            }
            break;

        case CellType::TET10:
            lf.type = CellType::TRI6;
            lf.size = 6;
            if (face_idx == 0)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[1];
                lf.v[2] = v[2];
                lf.v[3] = v[4];
                lf.v[4] = v[5];
                lf.v[5] = v[6];
            }
            else if (face_idx == 1)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[1];
                lf.v[2] = v[3];
                lf.v[3] = v[4];
                lf.v[4] = v[8];
                lf.v[5] = v[7];
            }
            else if (face_idx == 2)
            {
                lf.v[0] = v[1];
                lf.v[1] = v[2];
                lf.v[2] = v[3];
                lf.v[3] = v[5];
                lf.v[4] = v[9];
                lf.v[5] = v[8];
            }
            else
            {
                lf.v[0] = v[2];
                lf.v[1] = v[0];
                lf.v[2] = v[3];
                lf.v[3] = v[6];
                lf.v[4] = v[7];
                lf.v[5] = v[9];
            }
            break;

        case CellType::PRISM18:
            if (face_idx < 2)
            {
                lf.type = CellType::TRI6;
                lf.size = 6;
            }
            else
            {
                lf.type = CellType::QUAD9;
                lf.size = 9;
            }
            if (face_idx == 0)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[2];
                lf.v[2] = v[1];
                lf.v[3] = v[8];
                lf.v[4] = v[7];
                lf.v[5] = v[6];
            }
            else if (face_idx == 1)
            {
                lf.v[0] = v[3];
                lf.v[1] = v[4];
                lf.v[2] = v[5];
                lf.v[3] = v[9];
                lf.v[4] = v[10];
                lf.v[5] = v[11];
            }
            else if (face_idx == 2)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[1];
                lf.v[2] = v[4];
                lf.v[3] = v[3];
                lf.v[4] = v[6];
                lf.v[5] = v[13];
                lf.v[6] = v[9];
                lf.v[7] = v[12];
                lf.v[8] = v[15];
            }
            else if (face_idx == 3)
            {
                lf.v[0] = v[1];
                lf.v[1] = v[2];
                lf.v[2] = v[5];
                lf.v[3] = v[4];
                lf.v[4] = v[7];
                lf.v[5] = v[14];
                lf.v[6] = v[10];
                lf.v[7] = v[13];
                lf.v[8] = v[16];
            }
            else if (face_idx == 4)
            {
                lf.v[0] = v[2];
                lf.v[1] = v[0];
                lf.v[2] = v[3];
                lf.v[3] = v[5];
                lf.v[4] = v[8];
                lf.v[5] = v[12];
                lf.v[6] = v[11];
                lf.v[7] = v[14];
                lf.v[8] = v[17];
            }
            else{
                CFEM_ERROR << "face_idx out of bound.";
            }
            break;

        case CellType::HEX27:
            lf.type = CellType::QUAD9;
            lf.size = 9;
            if (face_idx == 0)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[3];
                lf.v[2] = v[2];
                lf.v[3] = v[1];
                lf.v[4] = v[11];
                lf.v[5] = v[10];
                lf.v[6] = v[9];
                lf.v[7] = v[8];
                lf.v[8] = v[24];
            }
            else if (face_idx == 1)
            {
                lf.v[0] = v[4];
                lf.v[1] = v[5];
                lf.v[2] = v[6];
                lf.v[3] = v[7];
                lf.v[4] = v[12];
                lf.v[5] = v[13];
                lf.v[6] = v[14];
                lf.v[7] = v[15];
                lf.v[8] = v[25];
            }
            else if (face_idx == 2)
            {
                lf.v[0] = v[0];
                lf.v[1] = v[1];
                lf.v[2] = v[5];
                lf.v[3] = v[4];
                lf.v[4] = v[8];
                lf.v[5] = v[17];
                lf.v[6] = v[12];
                lf.v[7] = v[16];
                lf.v[8] = v[22];
            }
            else if (face_idx == 3)
            {
                lf.v[0] = v[1];
                lf.v[1] = v[2];
                lf.v[2] = v[6];
                lf.v[3] = v[5];
                lf.v[4] = v[9];
                lf.v[5] = v[18];
                lf.v[6] = v[13];
                lf.v[7] = v[17];
                lf.v[8] = v[21];
            }
            else if (face_idx == 4)
            {
                lf.v[0] = v[2];
                lf.v[1] = v[3];
                lf.v[2] = v[7];
                lf.v[3] = v[6];
                lf.v[4] = v[10];
                lf.v[5] = v[19];
                lf.v[6] = v[14];
                lf.v[7] = v[18];
                lf.v[8] = v[23];
            }
            else if (face_idx == 5)
            {
                lf.v[0] = v[3];
                lf.v[1] = v[0];
                lf.v[2] = v[4];
                lf.v[3] = v[7];
                lf.v[4] = v[11];
                lf.v[5] = v[16];
                lf.v[6] = v[15];
                lf.v[7] = v[19];
                lf.v[8] = v[20];
            }
            else
            {
                CFEM_ERROR << "Index out of bounds. HEX Cell has exactely 6 faces. valid indices: 0,...,5.";
            }
            break;

        default:
            break;
        }
        return lf;
    }

    /** ************************************************************************************************************************************
     * @brief Returns the number of local faces for a given cell type.
     *
     * This function provides the number of facets (faces in 3D, edges in 2D, nodes in 1D) associated with a specific cell type. It is primarily used
     * when iterating over the faces of a cell during mesh topology construction, such as building facet lists or computing cell adjacency.
     *
     * The returned value corresponds to the number of valid indices that can be passed to Mesh::getLocalFaceData().
     *
     * ---
     *
     * ### Face counts by cell type
     *
     * | Cell Type        | Faces                                    |
     * |------------------|------------------------------------------|
     * | LINE2 / LINE3    | 2 nodes                                  |
     * | TRI3 / TRI6      | 3 edges                                  |
     * | QUAD4 / QUAD9    | 4 edges                                  |
     * | TET4 / TET10     | 4 triangular faces                       |
     * | PRISM6 / PRISM18 | 5 faces (2 triangular + 3 quadrilateral) |
     * | HEX8 / HEX27     | 6 quadrilateral faces                    |
     *
     * ---
     *
     * @param type The cell type whose number of faces is requested.
     *
     * @return The number of local faces associated with the cell type.
     *
     * @pre
     * The cell type must be supported by the mesh topology system.
     *
     * @note
     * Linear and quadratic versions of the same element share the same number of faces since higher-order nodes only add interpolation
     * points but do not change the topological structure.
     *
     * @warning
     * An error is triggered if an unsupported cell type is provided.
     *
     * @see getLocalFaceData()
     * @see buildFacets()
     * @see buildAdjacency()
     ***************************************************************************************************************************************/
    inline int Mesh::getLocalFacesCount(CellType type) const
    {
        switch (type)
        {
            case CellType::LINE2:
            case CellType::LINE3:
                return 2;
            case CellType::TRI3:
            case CellType::TRI6:
                return 3;
            case CellType::QUAD4:
            case CellType::QUAD9:
                return 4;
            case CellType::TET4:
            case CellType::TET10:
                return 4;
            case CellType::PRISM6:
            case CellType::PRISM18:
                return 5;
            case CellType::HEX8:
            case CellType::HEX27:
                return 6;

            default:
                CFEM_ERROR << "Unsupported cell type."; 
                return -1;
        }
    }

    /** ************************************************************************************************************************************
     * @brief Finds the local face index of a given face within a cell.
     *
     * This function searches the faces of a specified cell and determines which local face corresponds to a given target face. The comparison
     * is performed using the vertex indices of the face, ignoring ordering.
     *
     * The function iterates through all local faces of the cell, retrieves their connectivity using getLocalFaceData(), and checks whether the
     * vertices of the target face are contained in the candidate face.
     *
     * If a match is found, the corresponding local face index is returned.
     *
     * ---
     *
     * ### Matching rule
     *
     * Two faces are considered identical if:
     *
     * - They have the same number of vertices.
     * - All vertices of the target face are contained in the candidate face.
     *
     * The vertex ordering is **not required to match**.
     *
     * ---
     *
     * @param cell_id The ID of the cell in which the face search is performed.
     * @param target  The face definition whose local index is being searched.
     *
     * @return
     * - The local face index within the cell if a matching face is found.
     * - `-1` if no matching face exists in the cell.
     *
     * @pre
     * - `cell_id` must be a valid cell index.
     * - `target` must represent a valid face topology.
     *
     * @note
     * This function is typically used when reconstructing mesh topology, such as mapping global facets back to the local face index of a cell.
     *
     * @see getLocalFaceData()
     * @see getLocalFacesCount()
     * @see buildFacets()
     ***************************************************************************************************************************************/
    inline int Mesh::findLocalFaceIndex(CfemInt cell_id, const LocalFace &target) const
    {
        const int n = getLocalFacesCount(m_cells[cell_id].type);
        for (int i = 0; i < n; ++i)
        {
            LocalFace other = getLocalFaceData(cell_id, i);
            if (other.size != target.size)
                continue;

            bool match = true;
            for (int j = 0; j < target.size; ++j)
            {
                if (!other.contains(target.v[j]))
                {
                    match = false;
                    break;
                }
            }
            if (match)
                return i;
        }
        return -1;
    }

    /** ************************************************************************************************************************************
     * @brief Returns the midpoint vertex of an edge, creating it if necessary.
     *
     * This function retrieves or creates the midpoint vertex associated with the edge defined by vertices `v1` and `v2`. If the midpoint has already
     * been created earlier during the quadratic upgrade process, the stored vertex ID is returned directly from the edge map.
     *
     * Otherwise, a new vertex is created at the geometric midpoint of the two vertices and stored in the edge map to ensure uniqueness.
     *
     * The map guarantees that each edge produces **exactly one midpoint node**, even when the edge is shared by multiple cells.
     *
     * ---
     *
     * Midpoint coordinates are computed as:
     *
     * \f[
     * p_{mid} = \frac{p_1 + p_2}{2}
     * \f]
     *
     * ---
     *
     * @param v1 First vertex index of the edge.
     * @param v2 Second vertex index of the edge.
     * @param edge_map Hash map storing previously created edge midpoints.
     *
     * @return The vertex ID of the midpoint node.
     *
     * @note
     * Used during mesh upgrade to quadratic elements (e.g. TRI3 -> TRI6, TET4 -> TET10, HEX8 -> HEX27).
     *
     * @warning
     * The order of `(v1, v2)` must be handled consistently by the map's hashing and equality operators to avoid duplicate midpoints.
     *
     * @see upgradeToQuadratic()
     ***************************************************************************************************************************************/
    CfemInt Mesh::getMidpoint(CfemInt v1, CfemInt v2, auto &edge_map)
    {
        auto it = edge_map.find({v1, v2});
        if (it != edge_map.end())
            return it->second;

        const auto &n1 = m_vertices[v1];
        const auto &n2 = m_vertices[v2];

        // Add new vertex and return its ID
        CfemInt newId = addVertex(
            (n1.coords[0] + n2.coords[0]) * 0.5,
            (n1.coords[1] + n2.coords[1]) * 0.5,
            (n1.coords[2] + n2.coords[2]) * 0.5);
        edge_map[{v1, v2}] = newId;
        return newId;
    }

    /** ************************************************************************************************************************************
     * @brief Returns the center vertex of a quadrilateral face, creating it if necessary.
     *
     * This function retrieves or creates the center node associated with the quadrilateral face defined by the four vertices `(v1, v2, v3, v4)`.
     *
     * If the center node already exists in the provided map, its vertex ID is returned immediately. Otherwise, a new vertex is created at the
     * geometric center of the four vertices.
     *
     * The map ensures that each quadrilateral face produces **only one center node**, even when the face is shared between multiple cells.
     *
     * ---
     *
     * Face center coordinates are computed as:
     *
     * \f[
     * p_{center} = \frac{p_1 + p_2 + p_3 + p_4}{4}
     * \f]
     *
     * ---
     *
     * @param v1 First vertex of the quadrilateral face.
     * @param v2 Second vertex of the quadrilateral face.
     * @param v3 Third vertex of the quadrilateral face.
     * @param v4 Fourth vertex of the quadrilateral face.
     * @param quadface_map Hash map storing previously created face-center nodes.
     *
     * @return The vertex ID of the face-center node.
     *
     * @note
     * This function is used when upgrading elements such as:
     * - HEX8 -> HEX27
     * - PRISM6 -> PRISM18
     *
     * @see upgradeToQuadratic()
     ***************************************************************************************************************************************/
    CfemInt Mesh::getQuadFaceCenter(CfemInt v1, CfemInt v2, CfemInt v3, CfemInt v4, auto &quadface_map)
    {
        auto it = quadface_map.find({v1, v2, v3, v4});
        if (it != quadface_map.end())
            return it->second;

        const auto &n1 = m_vertices[v1];
        const auto &n2 = m_vertices[v2];
        const auto &n3 = m_vertices[v3];
        const auto &n4 = m_vertices[v4];

        // Add new vertex and return its ID
        CfemInt newId = addVertex( (n1.coords[0] + n2.coords[0] + n3.coords[0] + n4.coords[0]) * 0.25,
                                   (n1.coords[1] + n2.coords[1] + n3.coords[1] + n4.coords[1]) * 0.25,
                                   (n1.coords[2] + n2.coords[2] + n3.coords[2] + n4.coords[2]) * 0.25 );
        quadface_map[{v1, v2, v3, v4}] = newId;
        return newId;
    }

    /** ************************************************************************************************************************************
     * @brief Populate the node list for quadratic elements during mesh upgrade.
     *
     * These helper functions construct the connectivity of quadratic finite elements from the vertex connectivity of their linear counterparts.
     *
     * Each routine:
     *  - receives the vertex IDs of a linear cell
     *  - creates the required mid-edge, face-center, or volume-center nodes
     *  - writes the resulting quadratic node ordering into `out`
     *
     * Mid-edge and face-center nodes are cached using hash maps to ensure that nodes shared by adjacent cells are created only once.
     *
     * The resulting node ordering follows the CFEM internal conventions.
     *
     * Typical upgrades:
     *
     * LINE2  -> LINE3  
     * TRI3   -> TRI6  
     * QUAD4  -> QUAD9  
     * TET4   -> TET10  
     * PRISM6 -> PRISM18  
     * HEX8   -> HEX27
     *
     * @param v
     * Vertex indices of the original linear element.
     *
     * @param edge_to_midpoint
     * Hash map storing mid-edge nodes to guarantee uniqueness.
     *
     * @param quadface_to_midpoint_map
     * Hash map storing quadrilateral face-center nodes (used for prism/hex).
     *
     * @param out
     * Output vector that will contain the quadratic connectivity.
     *
     * @note
     * The output vector is cleared and filled with the quadratic node ordering.
     *
     * @see upgradeToQuadratic()
     * @see getMidpoint()
     * @see getQuadFaceCenter()
     ***************************************************************************************************************************************/
    inline void Mesh::fillLine3Nodes(std::span<const CfemInt> v, auto &edge_to_midpoint, std::vector<CfemInt> &out)
    {
        out.clear();
        // 01
        CfemInt n2 = getMidpoint(v[0], v[1], edge_to_midpoint);
        out = {v[0], v[1], n2};
    }

    /** ************************************************************************************************************************************
     * @brief Construct TRI6 connectivity from a TRI3 element.
     *
     * Generates the six-node quadratic triangle by inserting mid-edge nodes on each edge of the linear triangle.
     *
     * Node ordering:
     *
     * v0 -- n3 -- v1
     *  \         /
     *   n5     n4
     *    \    /
     *      v2
     *
     * Resulting connectivity:
     *
     * [ v0, v1, v2, n3, n4, n5 ]
     *
     * where:
     *  n3 = midpoint(v0, v1)
     *  n4 = midpoint(v1, v2)
     *  n5 = midpoint(v2, v0)
     *
     * @param v Linear triangle vertices (size = 3)
     * @param edge_to_midpoint Map storing mid-edge nodes
     * @param out Output quadratic node list (size = 6)
     ***************************************************************************************************************************************/
    inline void Mesh::fillTri6Nodes(std::span<const CfemInt> v, auto &edge_to_midpoint, std::vector<CfemInt> &out)
    {
        out.clear();
        // 01, 12, 20
        CfemInt n3 = getMidpoint(v[0], v[1], edge_to_midpoint);
        CfemInt n4 = getMidpoint(v[1], v[2], edge_to_midpoint);
        CfemInt n5 = getMidpoint(v[2], v[0], edge_to_midpoint);
        out = {v[0], v[1], v[2], n3, n4, n5};
    }

    /** ************************************************************************************************************************************
     * @brief Construct the node connectivity for a QUAD9 element from a QUAD4 cell.
     *
     * This function upgrades a linear quadrilateral element (QUAD4) to a quadratic quadrilateral (QUAD9) by adding mid-edge nodes and a
     * center node.
     *
     * The resulting quadratic element contains:
     * - 4 corner vertices
     * - 4 mid-edge vertices
     * - 1 face-center vertex
     *
     * Node layout:
     *
     * v3 ---- n6 ---- v2
     * |               |
     * n7      n8      n5
     * |               |
     * v0 ---- n4 ---- v1
     *
     * Resulting connectivity ordering:
     *
     * [ v0, v1, v2, v3, n4, n5, n6, n7, n8 ]
     *
     * where:
     *  n4 = midpoint(v0, v1)
     *  n5 = midpoint(v1, v2)
     *  n6 = midpoint(v2, v3)
     *  n7 = midpoint(v3, v0)
     *  n8 = geometric center of the quadrilateral
     *
     * Mid-edge nodes are reused across neighboring cells using the `edge_to_midpoint` map to avoid duplicate vertices.
     *
     * @param v Vertex indices of the QUAD4 element (size = 4)
     * @param edge_to_midpoint Hash map storing previously created edge midpoints
     * @param out Output vector that will contain the QUAD9 node ordering
     *
     * @note
     * The center node is always created locally since it belongs uniquely to the current cell.
     *
     * @see upgradeToQuadratic()
     * @see fillHex27Nodes()
     * @see fillPrism18Nodes()
     ***************************************************************************************************************************************/
    inline void Mesh::fillQuad9Nodes(std::span<const CfemInt> v, auto &edge_to_midpoint, std::vector<CfemInt> &out)
    {
        out.clear();
        CfemInt n4 = getMidpoint(v[0], v[1], edge_to_midpoint);
        CfemInt n5 = getMidpoint(v[1], v[2], edge_to_midpoint);
        CfemInt n6 = getMidpoint(v[2], v[3], edge_to_midpoint);
        CfemInt n7 = getMidpoint(v[3], v[0], edge_to_midpoint);

        cfem::Vector3D center = (m_vertices[v[0]].coords + m_vertices[v[1]].coords +
                                 m_vertices[v[2]].coords + m_vertices[v[3]].coords) *
                                0.25;
        CfemInt n8 = addVertex(center.x, center.y, center.z);

        out = {v[0], v[1], v[2], v[3], n4, n5, n6, n7, n8};
    }

    /** ************************************************************************************************************************************
     * @brief Construct the node connectivity for a TET10 element from a TET4 cell.
     *
     * This function upgrades a linear tetrahedral element (TET4) to a
     * quadratic tetrahedron (TET10) by inserting mid-edge nodes on each
     * of the six edges of the tetrahedron.
     *
     * The resulting quadratic element contains:
     * - 4 corner vertices
     * - 6 mid-edge vertices
     *
     * Edge layout:
     *
     *        v3
     *       /|\
     *     n7 | n9
     *     /  |  \
     *   v0---n6---v2
     *    \  n8   /
     *     n4 | n5
     *       \|/
     *        v1
     *
     * Resulting connectivity ordering:
     *
     * [ v0, v1, v2, v3, n4, n5, n6, n7, n8, n9 ]
     *
     * where:
     *  n4 = midpoint(v0, v1)
     *  n5 = midpoint(v1, v2)
     *  n6 = midpoint(v2, v0)
     *  n7 = midpoint(v0, v3)
     *  n8 = midpoint(v1, v3)
     *  n9 = midpoint(v2, v3)
     *
     * Mid-edge nodes are reused across neighboring cells through the `edge_to_midpoint` map to guarantee unique vertex creation for
     * shared edges.
     *
     * @param v Vertex indices of the linear TET4 element (size = 4)
     * @param edge_to_midpoint Hash map storing previously created edge midpoints
     * @param out Output vector that will contain the TET10 node ordering
     *
     * @note
     * TET10 elements do not contain face-center or volume-center nodes, only mid-edge nodes.
     *
     * @see upgradeToQuadratic()
     * @see fillTri6Nodes()
     * @see fillHex27Nodes()
     ***************************************************************************************************************************************/
    inline void Mesh::fillTet10Nodes(std::span<const CfemInt> v, auto &edge_to_midpoint, std::vector<CfemInt> &out)
    {
        out.clear();
        CfemInt n4 = getMidpoint(v[0], v[1], edge_to_midpoint);
        CfemInt n5 = getMidpoint(v[1], v[2], edge_to_midpoint);
        CfemInt n6 = getMidpoint(v[2], v[0], edge_to_midpoint);
        CfemInt n7 = getMidpoint(v[0], v[3], edge_to_midpoint);
        CfemInt n8 = getMidpoint(v[1], v[3], edge_to_midpoint);
        CfemInt n9 = getMidpoint(v[2], v[3], edge_to_midpoint);

        out = {v[0], v[1], v[2], v[3], n4, n5, n6, n7, n8, n9};
    }

    /** ************************************************************************************************************************************
     * @brief Construct the node connectivity for a PRISM18 element from a PRISM6 cell.
     *
     * This function upgrades a linear triangular prism (PRISM6) to its quadratic counterpart (PRISM18). The quadratic prism contains:
     *
     * - 6 corner vertices
     * - 9 mid-edge vertices
     * - 3 quadrilateral face-center vertices
     *
     * No center node exists for the triangular faces.
     *
     * The node layout follows the CFEM quadratic prism convention:
     *
     * Corners:
     *   v0, v1, v2  -> bottom triangle
     *   v3, v4, v5  -> top triangle
     *
     * Mid-edge nodes:
     *   n6  = midpoint(v0, v1)
     *   n7  = midpoint(v1, v2)
     *   n8  = midpoint(v2, v0)
     *
     *   n9  = midpoint(v3, v4)
     *   n10 = midpoint(v4, v5)
     *   n11 = midpoint(v5, v3)
     *
     *   n12 = midpoint(v0, v3)
     *   n13 = midpoint(v1, v4)
     *   n14 = midpoint(v2, v5)
     *
     * Quadrilateral face centers:
     *   n15 = center(v0, v1, v4, v3)
     *   n16 = center(v1, v2, v5, v4)
     *   n17 = center(v2, v0, v3, v5)
     *
     * Resulting connectivity ordering:
     *
     * [ v0, v1, v2, v3, v4, v5,
     *   n6, n7, n8, n9, n10, n11, n12, n13, n14,
     *   n15, n16, n17 ]
     *
     * Mid-edge nodes are cached using `edge_to_midpoint` so that edges shared by multiple cells produce a single unique midpoint node.
     *
     * Quadrilateral face centers are cached using `quadface_to_midpoint_map` to ensure shared faces reuse the same node.
     *
     * @param v Vertex indices of the linear PRISM6 element (size = 6)
     * @param edge_to_midpoint Hash map storing previously created edge midpoints
     * @param quadface_to_midpoint_map Hash map storing quadrilateral face-center nodes
     * @param out Output vector that will contain the PRISM18 node ordering
     *
     * @note
     * Triangular faces of the prism do not receive center nodes. Only quadrilateral faces introduce additional nodes.
     *
     * @see upgradeToQuadratic()
     * @see fillTri6Nodes()
     * @see fillHex27Nodes()
     ***************************************************************************************************************************************/
    inline void Mesh::fillPrism18Nodes(std::span<const CfemInt> v, auto &edge_to_midpoint, auto &quadface_to_midpoint_map, std::vector<CfemInt> &out)
    {
        out.clear();
        out.reserve(18);

        // 1. Bottom Triangle Edges (0-1, 1-2, 2-0)
        CfemInt n6 = getMidpoint(v[0], v[1], edge_to_midpoint);
        CfemInt n7 = getMidpoint(v[1], v[2], edge_to_midpoint);
        CfemInt n8 = getMidpoint(v[2], v[0], edge_to_midpoint);

        // 2. Top Triangle Edges (3-4, 4-5, 5-3)
        CfemInt n9 = getMidpoint(v[3], v[4], edge_to_midpoint);
        CfemInt n10 = getMidpoint(v[4], v[5], edge_to_midpoint);
        CfemInt n11 = getMidpoint(v[5], v[3], edge_to_midpoint);

        // 3. Vertical Edges (0-3, 1-4, 2-5)
        CfemInt n12 = getMidpoint(v[0], v[3], edge_to_midpoint);
        CfemInt n13 = getMidpoint(v[1], v[4], edge_to_midpoint);
        CfemInt n14 = getMidpoint(v[2], v[5], edge_to_midpoint);

        // 4. Quad Face Centers (3 faces)
        // We use a lambda that directly adds the new vertex to the Mesh
        CfemInt n15 = getQuadFaceCenter(v[0], v[1], v[4], v[3], quadface_to_midpoint_map);
        CfemInt n16 = getQuadFaceCenter(v[1], v[2], v[5], v[4], quadface_to_midpoint_map); // Face 1-2-5-4
        CfemInt n17 = getQuadFaceCenter(v[2], v[0], v[3], v[5], quadface_to_midpoint_map); // Face 2-0-3-5

        // 5. Fill the buffer in the standard order
        out = {
            v[0], v[1], v[2], v[3], v[4], v[5],      // Corners (0-5)
            n6, n7, n8, n9, n10, n11, n12, n13, n14, // Mid-edges (6-14)
            n15, n16, n17                            // Face centers (15-17)
        };
    }

    /** ************************************************************************************************************************************
     * @brief Construct HEX27 connectivity from a HEX8 element.
     *
     * Generates the 27-node quadratic hexahedron consisting of:
     *
     * - 8 corner nodes
     * - 12 mid-edge nodes
     * - 6 face-center nodes
     * - 1 volume-center node
     *
     * Node ordering:
     *
     * Corners:        0-7
     * Edge midpoints: 8-19
     * Face centers:   20-25
     * Volume center:  26
     *
     * The mid-edge and face-center nodes are reused across neighboring cells using hash maps to avoid duplication.
     *
     * The final node ordering follows the CFEM quadratic hexahedron convention.
     *
     * @param v Linear HEX8 vertex list (size = 8)
     * @param edge_to_midpoint Map storing mid-edge nodes
     * @param quadface_to_midpoint_map Map storing quad face centers
     * @param out Output quadratic node list (size = 27)
     ***************************************************************************************************************************************/
    inline void Mesh::fillHex27Nodes(std::span<const CfemInt> v, auto &edge_to_midpoint, auto &quadface_to_midpoint_map, std::vector<CfemInt> &out)
    {
        out.clear();
        out.reserve(27);

        // Mid-edges (12)
        CfemInt n8 = getMidpoint(v[0], v[1], edge_to_midpoint);
        CfemInt n9 = getMidpoint(v[1], v[2], edge_to_midpoint);
        CfemInt n10 = getMidpoint(v[2], v[3], edge_to_midpoint);
        CfemInt n11 = getMidpoint(v[3], v[0], edge_to_midpoint);
        CfemInt n12 = getMidpoint(v[4], v[5], edge_to_midpoint);
        CfemInt n13 = getMidpoint(v[5], v[6], edge_to_midpoint);
        CfemInt n14 = getMidpoint(v[6], v[7], edge_to_midpoint);
        CfemInt n15 = getMidpoint(v[7], v[4], edge_to_midpoint);
        CfemInt n16 = getMidpoint(v[0], v[4], edge_to_midpoint);
        CfemInt n17 = getMidpoint(v[1], v[5], edge_to_midpoint);
        CfemInt n18 = getMidpoint(v[2], v[6], edge_to_midpoint);
        CfemInt n19 = getMidpoint(v[3], v[7], edge_to_midpoint);

        CfemInt n20 = getQuadFaceCenter(v[0], v[4], v[7], v[3], quadface_to_midpoint_map);
        CfemInt n21 = getQuadFaceCenter(v[5], v[1], v[2], v[6], quadface_to_midpoint_map);
        CfemInt n22 = getQuadFaceCenter(v[0], v[4], v[5], v[1], quadface_to_midpoint_map);
        CfemInt n23 = getQuadFaceCenter(v[3], v[2], v[6], v[7], quadface_to_midpoint_map);
        CfemInt n24 = getQuadFaceCenter(v[0], v[1], v[2], v[3], quadface_to_midpoint_map);
        CfemInt n25 = getQuadFaceCenter(v[4], v[5], v[6], v[7], quadface_to_midpoint_map);

        cfem::Vector3D bSum({0, 0, 0});
        for (auto id : v)
            bSum = bSum + m_vertices[id].coords;
        CfemInt n26 = addVertex(bSum.x * 0.125, bSum.y * 0.125, bSum.z * 0.125);

        out = {v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7],
               n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19,
               n20, n21, n22, n23, n24, n25, n26};
    }

    /** ************************************************************************************************************************************
     * @brief Retrieves the physical coordinates of all nodes belonging to an element.
     * @param elementId The index of the element.
     * @param localCoords [out] A vector that will be filled with the cfem::Vector3D coordinates.
     ***************************************************************************************************************************************/
    inline void Mesh::getElementCoordinates(CfemInt elementId, cfem::DynamicVector<cfem::Vector3D> &localCoords) const
    {
        auto nodeIndices = getCellConnectivity(elementId);

        const int n = static_cast<int>(nodeIndices.size());

        localCoords.resize(n);

        for (int i = 0; i < n; ++i)
        {
            CfemInt globalIdx = nodeIndices[i];
            localCoords[i] = m_vertices[globalIdx].coords;
        }
    }

    /** ************************************************************************************************************************************
     * @brief Axis-aligned bounding box defining a planar quadrilateral region.
     *
     * QuadBounds represents a rectangular region aligned with one of the principal coordinate planes (XY, XZ, or YZ). It is typically used
     * to generate structured quadrilateral meshes on planar surfaces.
     *
     * The constructor validates that the bounding box:
     * - lies entirely within one coordinate plane
     * - is not degenerate (non-zero area)
     *
     * Exactly one coordinate direction must remain constant:
     *
     * XY plane -> z is constant  
     * XZ plane -> y is constant  
     * YZ plane -> x is constant
     *
     * @param p_min Lower corner of the bounding box.
     * @param p_max Upper corner of the bounding box.
     *
     * @throws CFEM_ERROR if:
     * - the bounds do not define a planar region
     * - the region is degenerate
     *
     * @note
     * This structure is typically used when generating structured surface meshes.
     *
     * @see QuadCorners
     ***************************************************************************************************************************************/
    struct Mesh::QuadBounds
    {
        Vector3D min;
        Vector3D max;

        QuadBounds(const Vector3D &p_min, const Vector3D &p_max) : min(p_min), max(p_max)
        {
            constexpr CfemReal eps = cfem::ZERO_TOLERANCE;

            bool x_equal = std::abs(min.x - max.x) < eps;
            bool y_equal = std::abs(min.y - max.y) < eps;
            bool z_equal = std::abs(min.z - max.z) < eps;

            CfemInt equal_count = x_equal + y_equal + z_equal;

            if (equal_count != 1)
                CFEM_ERROR << "QuadBounds: must lie in XY, YZ or XZ plane.";

            if (std::abs((max - min).norm()) < eps)
                CFEM_ERROR << "QuadBounds: degenerate geometry.";
        }
    };

    /** ************************************************************************************************************************************
     * @brief Axis-aligned bounding box describing a 3D block domain.
     *
     * BlockBounds represents a rectangular volume aligned with the Cartesian axes. It is typically used to generate structured
     * hexahedral meshes.
     *
     * The constructor automatically normalizes the bounds so that:
     *
     * min = component-wise minimum  
     * max = component-wise maximum
     *
     * The block must have a non-zero size in all three directions.
     *
     * @param p_min One corner of the block.
     * @param p_max Opposite corner of the block.
     *
     * @throws CFEM_ERROR if the block is degenerate (zero thickness along any axis).
     *
     * @note
     * This structure is commonly used for structured mesh generation of rectangular domains.
     *
     * @see BlockCorners
     ***************************************************************************************************************************************/
    struct Mesh::BlockBounds
    {
        Vector3D min;
        Vector3D max;

        BlockBounds(const Vector3D &p_min, const Vector3D &p_max)
        {
            constexpr CfemReal eps = cfem::ZERO_TOLERANCE;

            min = {std::min(p_min.x, p_max.x), std::min(p_min.y, p_max.y), std::min(p_min.z, p_max.z)};
            max = {std::max(p_min.x, p_max.x), std::max(p_min.y, p_max.y), std::max(p_min.z, p_max.z)};

            if (std::abs(max.x - min.x) < eps ||
                std::abs(max.y - min.y) < eps ||
                std::abs(max.z - min.z) < eps)
            {
                CFEM_ERROR << "BlockBounds: Degenerate block.";
            }
        }
    };

    /** ************************************************************************************************************************************
     * @brief Geometric representation of a quadrilateral defined by its corners.
     *
     * QuadCorners stores the four vertices defining a quadrilateral surface. Unlike QuadBounds, the quadrilateral is not required to be axis-aligned
     * and may lie on any plane in 3D space.
     *
     * The constructor verifies that:
     * - the vertices are not collinear
     * - the geometry is not degenerate
     *
     * This structure is typically used for generating meshes on arbitrarily oriented quadrilateral surfaces.
     *
     * Vertex ordering must follow the mesh convention for consistent face orientation.
     *
     * @param verts Array containing the four quadrilateral vertices.
     *
     * @throws CFEM_ERROR if the vertices define a degenerate geometry.
     *
     * @note
     * Only basic geometric validation is performed. Full planarity checks are assumed to be handled by mesh construction routines.
     ***************************************************************************************************************************************/
    struct Mesh::QuadCorners
    {
        std::array<Vector3D, 4> vertices;

        QuadCorners(const std::array<Vector3D, 4> &verts) : vertices(verts)
        {
            // Check coplanarity
            Vector3D v0v1 = vertices[1] - vertices[0];
            Vector3D v0v2 = vertices[2] - vertices[0];
            Vector3D normal_0 = v0v1.cross(v0v2);

            Vector3D v3v1 = vertices[1] - vertices[3];
            Vector3D v3v2 = vertices[2] - vertices[3];
            Vector3D normal_3 = v3v1.cross(v3v2);

            if ((normal_0.norm() < cfem::ZERO_TOLERANCE) || (normal_3.norm() < cfem::ZERO_TOLERANCE))
            {
                CFEM_ERROR << "QuadCorners: vertices are degenerate or collinear.";
            }
        }

        Vector3D operator[](size_t i) const
        {
            if (!(i >= 0 && i < 4))
            {
                CFEM_ERROR << "Index out of bound.";
            }

            return vertices[i];
        }
    };

    /** ************************************************************************************************************************************
     * @brief Geometric representation of a hexahedral block defined by its corners.
     *
     * BlockCorners stores the eight vertices defining a hexahedral volume. This allows the creation of structured meshes on arbitrary hexahedral
     * domains, including skewed or warped blocks.
     *
     * The constructor performs geometric validation:
     *
     * - The hexahedron must have a positive non-zero volume.
     * - Each face must be approximately coplanar.
     *
     * Volume is computed by decomposing the hexahedron into tetrahedra.
     *
     * Face coplanarity is verified using the distance of the fourth vertex to the plane defined by the first three vertices.
     *
     * @param verts Array containing the eight vertices of the hexahedron.
     *
     * @throws CFEM_ERROR if:
     * - the hexahedron has zero or negative volume
     * - any face is not approximately planar
     *
     * @note
     * This structure allows generation of structured meshes on general hexahedral geometries rather than only axis-aligned blocks.
     *
     * @warning
     * Vertex ordering must follow the internal hexahedral convention used by the mesh generator.
     ***************************************************************************************************************************************/
    struct Mesh::BlockCorners
    {
        std::array<Vector3D, 8> vertices;

        BlockCorners(const std::array<Vector3D, 8> &verts) : vertices(verts)
        {
            double vol = compute_volume();
            if (vol < cfem::ZERO_TOLERANCE)
                CFEM_ERROR << "BlockCorners: degenerate hexahedron with zero or negative volume.";

            // Check each face is roughly coplanar
            check_faces_coplanar();
        }

        Vector3D operator[](size_t i) const { return vertices[i]; }

    private:
        // Compute signed volume of tetrahedron defined by points a,b,c,d
        static double tetrahedron_volume(const Vector3D &a,
                                         const Vector3D &b,
                                         const Vector3D &c,
                                         const Vector3D &d)
        {
            return ((b - a).cross(c - a)).dot(d - a) / 6.0;
        }

        double compute_volume() const
        {
            const auto &v = vertices;
            double vol = 0.0;
            vol += tetrahedron_volume(v[0], v[1], v[3], v[4]);
            vol += tetrahedron_volume(v[1], v[2], v[3], v[6]);
            vol += tetrahedron_volume(v[1], v[5], v[4], v[6]);
            vol += tetrahedron_volume(v[3], v[4], v[7], v[6]);
            vol += tetrahedron_volume(v[1], v[3], v[4], v[6]);
            return std::abs(vol);
        }

        // Check that 4 nodes of a face are approximately coplanar
        static void check_coplanar(const Vector3D &a,
                                   const Vector3D &b,
                                   const Vector3D &c,
                                   const Vector3D &d,
                                   double eps = 1e-12)
        {
            Vector3D normal = (b - a).cross(c - a);
            double dist = std::abs((d - a).dot(normal)) / normal.norm();
            if (dist > eps)
                CFEM_ERROR << "Corners<3>: face is not coplanar";
        }

        void check_faces_coplanar() const
        {
            const auto &v = vertices;
            // 6 faces of a hexahedron, each defined by 4 vertices
            // Face order is arbitrary but must match your mesh convention

            // Bottom face: 0,1,2,3
            check_coplanar(v[0], v[1], v[2], v[3]);

            // Top face: 4,5,6,7
            check_coplanar(v[4], v[5], v[6], v[7]);

            // Front face: 0,1,5,4
            check_coplanar(v[0], v[1], v[5], v[4]);

            // Back face: 2,3,7,6
            check_coplanar(v[2], v[3], v[7], v[6]);

            // Left face: 0,3,7,4
            check_coplanar(v[0], v[3], v[7], v[4]);

            // Right face: 1,2,6,5
            check_coplanar(v[1], v[2], v[6], v[5]);
        }
    };


    /** ************************************************************************************************************************************
     * @brief Compute the measure of a cell.
     *
     * Computes the geometric measure of a cell:
     *
     * | Dimension | Measure |
     * |-----------|---------|
     * | 1D        | Length  |
     * | 2D        | Area    |
     * | 3D        | Volume  |
     *
     * @param type Cell type.
     * @param v_ids Vertex indices of the element.
     *
     * @return The geometric measure of the cell.
     ***************************************************************************************************************************************/
    template<typename T>
    inline CfemReal Mesh::computeCellMeasure(CellType type, const T &v_ids) const
    {
        switch (type)
        {
        case CellType::LINE2:
            return compute2NodesLineLength(v_ids[0], v_ids[1]);

        case CellType::TRI3:
        case CellType::TRI6:
        {
            return computeSigned3NodesTriangleArea(v_ids[0], v_ids[1], v_ids[2]);
        }

        case CellType::QUAD4:
        case CellType::QUAD9:  // use this temporarery
        {
            cfem::Vector3D lReferenceNormale{0.0, 0.0, 0.0};
            CfemReal A1 = computeSigned3NodesTriangleArea(v_ids[0], v_ids[1], v_ids[2]);
            CfemReal A2 = computeSigned3NodesTriangleArea(v_ids[0], v_ids[2], v_ids[3]);

            if ((std::abs(A1) <= cfem::ZERO_TOLERANCE) || (std::abs(A2) <= cfem::ZERO_TOLERANCE))
            {
                CFEM_WARN << "Sub-triangle in QUAD4 has zero area. Check if nodes are colinear.";
            }

            if ((A1 > cfem::ZERO_TOLERANCE && A2 < -cfem::ZERO_TOLERANCE) ||
                (A1 < -cfem::ZERO_TOLERANCE && A2 > cfem::ZERO_TOLERANCE))
            {
                CFEM_WARN << "Poor quality QUAD4 cell: Self-intersecting geometry detected.";
            }

            CfemReal totalArea = A1 + A2;
            return totalArea;
        }

        case CellType::TET4:
        case CellType::TET10:
            return computeSigned4NodesTetraVolume(v_ids[0], v_ids[1], v_ids[2], v_ids[3]);

        case CellType::PRISM6:
        case CellType::PRISM18:
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

                if (tetra_vol > 0)
                    positiveTets++;
                else if (tetra_vol < 0)
                    negativeTets++;

                lTotalVolume += tetra_vol;
            }

            if (positiveTets > 0 && negativeTets > 0)
            {
                CFEM_WARN << "Invalid PRISM6: Twisted/Self-intersecting geometry detected.";
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
                    CFEM_WARN << "Invalid HEX8: Sub-tetra in HEX8 has zero volume. Check if nodes are coplanar.";
                }

                if (tetra_vol > 0)
                    positiveTets++;
                else if (tetra_vol < 0)
                    negativeTets++;

                lTotalVolume += tetra_vol;
            }

            if (positiveTets > 0 && negativeTets > 0)
            {
                CFEM_WARN << "Invalid HEX8: Twisted/Self-intersecting geometry detected.";
            }

            return lTotalVolume;
        }

        default:
            CFEM_ERROR << "Measure calculation not implemented for " << cellTypeToString(type);
            return 0;
        }
    }

    /** ************************************************************************************************************************************
     * @brief Compute the Euclidean distance between two mesh vertices.
     *
     * This function returns the physical length of the segment connecting two vertices in the mesh.
     *
     * The distance is computed using the standard Euclidean norm:
     *
     * |x₁ − x₀| = sqrt((dx)² + (dy)² + (dz)²)
     *
     * where:
     * dx = x₁ - x₀
     *
     * @param p_idx0 Index of the first vertex.
     * @param p_idx1 Index of the second vertex.
     *
     * @return Length of the line segment between the two vertices.
     *
     * @note
     * This function is used internally by several geometric operations such as element quality metrics and mesh validation.
     ***************************************************************************************************************************************/
    inline CfemReal Mesh::compute2NodesLineLength(CfemInt p_idx0, CfemInt p_idx1) const
    {
        const auto &x0 = m_vertices[p_idx0].coords;
        const auto &x1 = m_vertices[p_idx1].coords;

        const CfemReal dx = x1[0] - x0[0];
        const CfemReal dy = x1[1] - x0[1];
        const CfemReal dz = x1[2] - x0[2];

        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

    /** ************************************************************************************************************************************
     * @brief Compute the signed area of a triangle defined by three vertices.
     *
     * The triangle area is computed using the magnitude of the cross product:
     *
     * A = 0.5 * | (B − A) × (C − A) |
     *
     * If a global reference normal is provided (`m_global_reference_normal`), the function returns a **signed area** based on the orientation of
     * the triangle relative to that reference direction:
     *
     * A_signed = 0.5 * (n_ref · ((B − A) × (C − A)))
     *
     * Otherwise, the function returns the **positive geometric area**.
     *
     * @param p_idx0 Index of vertex A.
     * @param p_idx1 Index of vertex B.
     * @param p_idx2 Index of vertex C.
     *
     * @return
     * - Signed triangle area if a reference normal exists.
     * - Positive triangle area otherwise.
     *
     * @note
     * If the triangle is degenerate (area below ZERO_TOLERANCE),
     * the function returns 0.
     *
     * @warning
     * The sign of the area depends on the vertex ordering when a reference normal is used.
     ***************************************************************************************************************************************/
    inline CfemReal Mesh::computeSigned3NodesTriangleArea(CfemInt p_idx0,
                                                          CfemInt p_idx1,
                                                          CfemInt p_idx2) const
    {
        const auto &A = m_vertices[p_idx0].coords;
        const auto &B = m_vertices[p_idx1].coords;
        const auto &C = m_vertices[p_idx2].coords;

        cfem::Vector3D AB = {B[0] - A[0], B[1] - A[1], B[2] - A[2]};
        cfem::Vector3D AC = {C[0] - A[0], C[1] - A[1], C[2] - A[2]};

        // Cross product
        cfem::Vector3D cp = {
            AB[1] * AC[2] - AB[2] * AC[1],
            AB[2] * AC[0] - AB[0] * AC[2],
            AB[0] * AC[1] - AB[1] * AC[0]
        };

        CfemReal double_area_mag = std::sqrt(cp[0] * cp[0] + cp[1] * cp[1] + cp[2] * cp[2]);

        if (double_area_mag < cfem::ZERO_TOLERANCE) return 0.0;

        // If global reference is set, use it for the sign
        if (m_global_reference_normal.has_value()) {
            const auto& ref = *m_global_reference_normal;
            return 0.5 * (ref[0] * cp[0] + ref[1] * cp[1] + ref[2] * cp[2]);
        }

        // No reference? Return the physical magnitude (always positive)
        return 0.5 * double_area_mag;
    }

    /** ************************************************************************************************************************************
     * @brief Compute the signed volume of a tetrahedron defined by four vertices.
     *
     * The volume is calculated using the scalar triple product:
     *
     * V = ((p1 − p0) · ((p2 − p0) × (p3 − p0))) / 6
     *
     * The result is **signed**:
     *
     * - Positive volume indicates a right-handed orientation.
     * - Negative volume indicates an inverted tetrahedron.
     * - Zero volume indicates a degenerate configuration.
     *
     * @param p_idx0 Index of vertex p0.
     * @param p_idx1 Index of vertex p1.
     * @param p_idx2 Index of vertex p2.
     * @param p_idx3 Index of vertex p3.
     *
     * @return Signed tetrahedral volume.
     *
     * @note
     * This function is primarily used internally for:
     * - element orientation checks
     * - mesh validity tests
     * - geometric measure computations.
     *
     * @warning
     * A negative result indicates an inverted element.
     ***************************************************************************************************************************************/
    inline CfemReal Mesh::computeSigned4NodesTetraVolume(CfemInt p_idx0, CfemInt p_idx1, CfemInt p_idx2, CfemInt p_idx3) const
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


    /** ************************************************************************************************************************************
     * @brief Project the vertices of a boundary entity onto a sphere surface.
     *
     * This function moves all vertices of the specified boundary entity so that they lie on the surface of a sphere defined by a center and radius.
     *
     * The projection is done along the vector from the sphere center to the vertex:
     * 
     *   v_new = center + (v_old - center) / ||v_old - center|| * radius
     *
     * Vertices that are already extremely close to the center (distance < ZERO_TOLERANCE) are left unchanged to avoid division by zero.
     *
     * @param mesh Reference to the mesh containing the boundary.
     * @param center Sphere center as a 3D vector.
     * @param radius Sphere radius.
     *
     * @note
     * - `cfem::SYS_BOUNDARY_ID` is assumed to identify the boundary entity.
     * - This function only affects vertices in the boundary entity.
     * - Useful for generating curved meshes or projecting mid-edge nodes onto
     *   a spherical surface.
     *
     * @warning
     * Ensure that the `boundary` entity exists in the mesh. If not, the function silently returns without modifying any vertices.
     ***************************************************************************************************************************************/
    inline void Mesh::projectBoundaryMidNodes(Mesh &mesh, const cfem::Vector3D &center, CfemReal radius)
    {
        MeshEntity *boundary = mesh.getEntity(cfem::SYS_BOUNDARY_ID);
        if (!boundary)
            return;

        const auto &boundaryVids = boundary->getVertexIds();

        for (CfemInt vid : boundaryVids)
        {
            Vertex &v = mesh.getVertex(vid);

            cfem::Vector3D vec = v.coords - center;
            CfemReal currentDist = vec.norm();

            if (currentDist > cfem::ZERO_TOLERANCE)
            {
                v.coords = center + (vec / currentDist) * radius;
            }
        }
    }

    /** ************************************************************************************************************************************
     * @brief Smooth/relax the mesh using Laplacian relaxation.
     *
     * This function iteratively moves interior vertices toward the centroid of their neighbors, effectively "relaxing" the mesh to improve element quality.
     *
     * Boundary vertices are fixed and do not move.
     *
     * @param iterations Number of relaxation iterations to perform.
     * @param lambda Relaxation factor (step size) in range [0.1, 0.9].
     *
     * @details
     * Algorithm steps:
     * 1. Identify boundary vertices (using `cfem::SYS_BOUNDARY_ID`) and mark them as fixed.
     * 2. For each iteration:
     *    - For every interior vertex:
     *      - Compute the average position of all neighboring vertices (vertices that share at least one cell with this vertex).
     *      - Update the vertex position as:
     *          v_new = v_old + lambda * (avg_neighbor - v_old)
     *      - The factor lambda controls how aggressively the vertex moves toward its neighbors. Smaller lambda = slower, more stable movement.
     *    - All vertex updates are stored in a temporary buffer and applied after each iteration to avoid order-dependence.
     *
     * @note
     * - The method uses a simple Laplacian smoothing approach (a.k.a. Jacobi relaxation).
     * - It is parallelized with OpenMP.
     * - Boundary vertices remain fixed throughout the relaxation.
     * - The neighbor search uses all cells sharing the vertex, so it works for arbitrary unstructured meshes.
     ***************************************************************************************************************************************/
    inline void Mesh::relaxMesh(CfemInt iterations, CfemReal lambda)
    {
        PROFILE_SCOPE("Mesh::relaxMesh");

        const size_t numVertices = getNumVertices();
        if (numVertices == 0 || m_cells.empty())
            return;

        // 1. Identify boundary vertices (O(N))
        std::vector<bool> isBoundary(numVertices, false);
        if (MeshEntity *boundary = getEntity(cfem::SYS_BOUNDARY_ID))
        {
            for (CfemInt id : boundary->getVertexIds())
            {
                if (id >= 0 && (size_t)id < numVertices)
                    isBoundary[id] = true;
            }
        }

        // 2. Buffers for updates
        std::vector<cfem::Vector3D> nextCoords(numVertices);
        lambda = std::clamp(lambda, 0.1, 0.9);

        // 3. Iterations
        for (int iter = 0; iter < iterations; ++iter)
        {
#pragma omp parallel for schedule(static)
            for (long v_id = 0; v_id < (long)numVertices; ++v_id)
            {
                if (isBoundary[v_id])
                {
                    nextCoords[v_id] = m_vertices[v_id].coords;
                    continue;
                }

                const auto &vertex = m_vertices[v_id];
                cfem::Vector3D sum{0, 0, 0};
                int neighborCount = 0;

                // --- THE CSR APPROACH ---
                // Instead of a set, look at all cells sharing this vertex
                for (CfemInt i = 0; i < vertex.num_connected_cells; ++i)
                {
                    CfemInt c_id = vertex.connected_cells_ptr[i];
                    const auto &cell_vids = getCellConnectivity(c_id);

                    // In a cell, neighbors are the nodes adjacent to v_id
                    // For simplicity/robustness in any cell type:
                    // average ALL other nodes in the sharing cells
                    for (CfemInt other_v : cell_vids)
                    {
                        if (other_v != v_id)
                        {
                            sum = sum + m_vertices[other_v].coords;
                            neighborCount++;
                        }
                    }
                }

                if (neighborCount > 0)
                {
                    cfem::Vector3D avg = sum * (1.0 / neighborCount);
                    nextCoords[v_id] = m_vertices[v_id].coords + lambda * (avg - m_vertices[v_id].coords);
                }
                else
                {
                    nextCoords[v_id] = m_vertices[v_id].coords;
                }
            }

// 4. Update coordinates
#pragma omp parallel for schedule(static)
            for (long i = 0; i < (long)numVertices; ++i)
            {
                if (!isBoundary[i])
                    m_vertices[i].coords = nextCoords[i];
            }
        }
    }

} // namespace