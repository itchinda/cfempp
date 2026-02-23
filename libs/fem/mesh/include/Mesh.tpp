#ifndef CFEM_MESH_TPP
#define CFEM_MESH_TPP
#include "Mesh.h"
#include "DynamicVector.h"
#endif

namespace cfem
{

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

    inline CfemInt getNodeCountPerCell(CellType type)
    {
        switch (type)
        {
        case CellType::LINE2:
            return 2;
        case CellType::TRI3:
            return 3;
        case CellType::QUAD4:
            return 4;
        case CellType::TET4:
            return 4;
        case CellType::PRISM6:
            return 6;
        case CellType::HEX8:
            return 8;
        case CellType::LINE3:
            return 3;
        case CellType::TRI6:
            return 6;
        case CellType::QUAD9:
            return 9;
        case CellType::TET10:
            return 10;
        case CellType::PRISM18:
            return 18;
        case CellType::HEX27:
            return 27;
        default:
            return 0;
        }
    }

    /**
     * @brief Internal helper to define face-to-vertex topology based on CellType.
     */
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
                lf.v[0] = v[0];
                lf.v[1] = v[1];
            }
            else if (face_idx == 1)
            {
                lf.v[0] = v[1];
                lf.v[1] = v[2];
            }
            else
            {
                lf.v[0] = v[2];
                lf.v[1] = v[0];
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
                CFEM_ERROR("face_idx out of bound.");
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
            else
            {
                lf.v[0] = v[3];
                lf.v[1] = v[0];
                lf.v[2] = v[4];
                lf.v[3] = v[7];
                lf.v[4] = v[11];
                lf.v[5] = v[16];
                lf.v[6] = v[15];
                lf.v[7] = v[19];
                lf.v[8] = v[25];
            }
            break;

        default:
            break;
        }
        return lf;
    }

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
                CFEM_ERROR("Unsupported cell type."); 
                return -1;
        }
    }

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

    inline void Mesh::fillLine3Nodes(std::span<const CfemInt> v, auto &edge_to_midpoint, std::vector<CfemInt> &out)
    {
        out.clear();
        // 01
        CfemInt n2 = getMidpoint(v[0], v[1], edge_to_midpoint);
        out = {v[0], v[1], n2};
    }

    inline void Mesh::fillTri6Nodes(std::span<const CfemInt> v, auto &edge_to_midpoint, std::vector<CfemInt> &out)
    {
        out.clear();
        // 01, 12, 20
        CfemInt n3 = getMidpoint(v[0], v[1], edge_to_midpoint);
        CfemInt n4 = getMidpoint(v[1], v[2], edge_to_midpoint);
        CfemInt n5 = getMidpoint(v[2], v[0], edge_to_midpoint);
        out = {v[0], v[1], v[2], n3, n4, n5};
    }

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

    /**
     * @brief Retrieves the physical coordinates of all nodes belonging to an element.
     * @param elementId The index of the element.
     * @param localCoords [out] A vector that will be filled with the cfem::Vector3D coordinates.
     */
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
                CFEM_ERROR("QuadBounds: must lie in XY, YZ or XZ plane.");

            if (std::abs((max - min).norm()) < eps)
                CFEM_ERROR("QuadBounds: degenerate geometry.");
        }
    };

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
                CFEM_ERROR("BlockBounds: Degenerate block.");
            }
        }
    };

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
                CFEM_ERROR("QuadCorners: vertices are degenerate or collinear.");
            }
        }

        Vector3D operator[](size_t i) const
        {
            if (!(i >= 0 && i < 4))
            {
                CFEM_ERROR("Index out of bound.");
            }

            return vertices[i];
        }
    };

    struct Mesh::BlockCorners
    {
        std::array<Vector3D, 8> vertices;

        BlockCorners(const std::array<Vector3D, 8> &verts) : vertices(verts)
        {
            double vol = compute_volume();
            if (vol < cfem::ZERO_TOLERANCE)
                CFEM_ERROR("BlockCorners: degenerate hexahedron with zero or negative volume.");

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
                CFEM_ERROR("Corners<3>: face is not coplanar");
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
                CFEM_WARNING("Sub-triangle in QUAD4 has zero area. Check if nodes are colinear.");
            }

            if ((A1 > cfem::ZERO_TOLERANCE && A2 < -cfem::ZERO_TOLERANCE) ||
                (A1 < -cfem::ZERO_TOLERANCE && A2 > cfem::ZERO_TOLERANCE))
            {
                CFEM_WARNING("Poor quality QUAD4 cell: Self-intersecting geometry detected.");
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

            return lTotalVolume;
        }

        default:
            CFEM_ERROR("Measure calculation not implemented for " + cellTypeToString(type));
        }
    }

    inline CfemReal Mesh::compute2NodesLineLength(CfemInt p_idx0, CfemInt p_idx1) const
    {
        const auto &x0 = m_vertices[p_idx0].coords;
        const auto &x1 = m_vertices[p_idx1].coords;

        const CfemReal dx = x1[0] - x0[0];
        const CfemReal dy = x1[1] - x0[1];
        const CfemReal dz = x1[2] - x0[2];

        return std::sqrt(dx * dx + dy * dy + dz * dz);
    }

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

    // Internal helper for measure calculation
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



    /**
     * @brief
     *
     * @param mesh
     * @param center
     * @param radius
     */
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

    /**
     * @brief
     *
     * @param iterations
     * @param lambda
     */
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