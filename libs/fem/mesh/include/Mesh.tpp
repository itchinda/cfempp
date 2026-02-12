#ifndef CFEM_MESH_TPP
#define CFEM_MESH_TPP
#include "Mesh.h"
#include "DynamicVector.h"
#endif

namespace cfem {

inline void Mesh::getCellVerticesCoords(CfemInt cell_id, DynamicVector<Vector3D>& v_coords) const{
    CFEM_ASSERT(cell_id >=0 && cell_id < getNumCells());
    const Cell& cell = getCell(cell_id);
    size_t numVertices {cell.vertex_ids.size()};
    v_coords.resize(numVertices);

    for (size_t i = 0; i < numVertices; ++i){
        v_coords[i] = getVertex(cell.vertex_ids[i]).coords;
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
inline std::vector< std::pair< std::vector<CfemInt>, CellType>> getLocalFaces(const Cell &cell)
{
    const auto &v = cell.vertex_ids;
    CFEM_ASSERT(v.size() == getNodeCountPerCell(cell.type));
    switch (cell.type)
    {
    case CellType::LINE2:
    case CellType::LINE3:
        return { { {v[0]}, CellType::NODE}, {{v[1]}, CellType::NODE} } ;
    case CellType::TRI3:
        return { { {v[0], v[1]}, CellType::LINE2}, {{v[1], v[2]} , CellType::LINE2}, {{v[2], v[0]} , CellType::LINE2}};
    case CellType::QUAD4:
        return {{{v[0], v[1]}, CellType::LINE2}, {{v[1], v[2]}, CellType::LINE2}, {{v[2], v[3]}, CellType::LINE2}, {{v[3], v[0]}, CellType::LINE2}};
    case CellType::TET4:
        return {{{v[0], v[1], v[2]}, CellType::TRI3},
                {{v[0], v[1], v[3]}, CellType::TRI3},
                {{v[1], v[2], v[3]}, CellType::TRI3},
                {{v[2], v[0], v[3]}, CellType::TRI3}};
    case CellType::PRISM6:
        return {{{v[0], v[2], v[1]}, CellType::TRI3},
                {{v[3], v[4], v[5]}, CellType::TRI3},
                {{v[0], v[1], v[4], v[3]}, CellType::QUAD4},
                {{v[1], v[2], v[5], v[4]}, CellType::QUAD4},
                {{v[2], v[0], v[3], v[5]}, CellType::QUAD4}};
    case CellType::HEX8:
        return {{{v[0], v[3], v[2], v[1]}, CellType::QUAD4},
                {{v[4], v[5], v[6], v[7]}, CellType::QUAD4},
                {{v[0], v[1], v[5], v[4]}, CellType::QUAD4},
                {{v[1], v[2], v[6], v[5]}, CellType::QUAD4},
                {{v[2], v[3], v[7], v[6]}, CellType::QUAD4},
                {{v[3], v[0], v[4], v[7]}, CellType::QUAD4}};
    // Quad cells?
    case CellType::TRI6:
        return {{{v[0], v[1], v[3]}, CellType::LINE3}, {{v[1], v[2], v[4]}, CellType::LINE3}, {{v[2], v[0], v[5]}, CellType::LINE3}};
    case CellType::QUAD9:
        return {{{v[0], v[1], v[4]}, CellType::LINE3}, {{v[1], v[2], v[5]}, CellType::LINE3}, 
                {{v[2], v[3], v[6]}, CellType::LINE3}, {{v[3], v[0], v[7]}, CellType::LINE3}};
    case CellType::TET10:
        return {{{v[0], v[1], v[2], v[4], v[5], v[6]}, CellType::TRI6},
                {{v[0], v[1], v[3], v[4], v[8], v[7]}, CellType::TRI6},
                {{v[1], v[2], v[3], v[5], v[9], v[8]}, CellType::TRI6},
                {{v[2], v[0], v[3], v[6], v[7], v[9]}, CellType::TRI6}};
    case CellType::PRISM18:
        return {{{v[0], v[2], v[1], v[8], v[7], v[6]}, CellType::TRI6},
                {{v[3], v[4], v[5], v[9], v[10], v[11]}, CellType::TRI6},
                {{v[0], v[1], v[4], v[3], v[6], v[13], v[9], v[12], v[15]}, CellType::QUAD9},
                {{v[1], v[2], v[5], v[4], v[7], v[14], v[10], v[13], v[16]}, CellType::QUAD9},
                {{v[2], v[0], v[3], v[5], v[8], v[12], v[11], v[14], v[17]}, CellType::QUAD9}};
    case CellType::HEX27:
        return {{{v[0], v[3], v[2], v[1], v[11], v[10], v[9], v[8], v[24]}, CellType::QUAD9},
                {{v[4], v[5], v[6], v[7], v[12], v[13], v[14], v[15], v[25]}, CellType::QUAD9},
                {{v[0], v[1], v[5], v[4], v[8], v[17], v[12], v[16], v[22]}, CellType::QUAD9},
                {{v[1], v[2], v[6], v[5], v[9], v[18], v[13], v[17], v[21]}, CellType::QUAD9},
                {{v[2], v[3], v[7], v[6], v[10], v[19], v[14], v[18], v[23]}, CellType::QUAD9},
                {{v[3], v[0], v[4], v[7], v[11], v[16], v[15], v[19], v[25]}, CellType::QUAD9}};

    default:
        CFEM_ERROR("UNKNOWN CELL TYPE.");
        return {};
    }
}

inline CfemInt Mesh::getMidpoint(CfemInt p_edgeFirstVertexId,
                                 CfemInt p_edgeSecondVertexId,
                                 std::map<std::pair<CfemInt, CfemInt>, CfemInt> &edge_to_midpoint)
{

    auto edge_key = std::minmax(p_edgeFirstVertexId, p_edgeSecondVertexId);
    auto it = edge_to_midpoint.find(edge_key);
    if (it != edge_to_midpoint.end())
    {
        return it->second;
    }
    // Not found: Calculate physical midpoint and add vertex
    cfem::Vector3D p0 = m_vertices[p_edgeFirstVertexId].coords;
    cfem::Vector3D p1 = m_vertices[p_edgeSecondVertexId].coords;
    cfem::Vector3D mid = (p0 + p1) * 0.5;
    CfemInt new_id = addVertex(mid.x, mid.y, mid.z);
    edge_to_midpoint[edge_key] = new_id;

    return new_id;
}

inline void Mesh::upgradeLine2ToLine3(Cell &pLineCell,
                                      std::map<std::pair<CfemInt, CfemInt>, CfemInt> &edge_to_midpoint)
{
    std::vector<CfemInt> v = pLineCell.vertex_ids;
    pLineCell.vertex_ids.push_back(getMidpoint(v[0], v[1], edge_to_midpoint));
    pLineCell.type = CellType::LINE3;
}

inline void Mesh::upgradeTri3ToTri6(Cell &pTriCell,
                                    std::map<std::pair<CfemInt, CfemInt>, CfemInt> &edge_to_midpoint)
{
    std::vector<CfemInt> v = pTriCell.vertex_ids; // Copy original linear nodes
    // New nodes at mid-edges 01, 12, 20
    CfemInt n3 = getMidpoint(v[0], v[1], edge_to_midpoint);
    CfemInt n4 = getMidpoint(v[1], v[2], edge_to_midpoint);
    CfemInt n5 = getMidpoint(v[2], v[0], edge_to_midpoint);
    pTriCell.vertex_ids = {v[0], v[1], v[2], n3, n4, n5};
    pTriCell.type = CellType::TRI6;
}

inline void Mesh::upgradeQuad4ToQuad9(Cell &pQuadCell,
                                      std::map<std::pair<CfemInt, CfemInt>, CfemInt> &edge_to_midpoint)
{
    std::vector<CfemInt> v = pQuadCell.vertex_ids;
    // Mid-edges 01, 12, 23, 30
    CfemInt n4 = getMidpoint(v[0], v[1], edge_to_midpoint);
    CfemInt n5 = getMidpoint(v[1], v[2], edge_to_midpoint);
    CfemInt n6 = getMidpoint(v[2], v[3], edge_to_midpoint);
    CfemInt n7 = getMidpoint(v[3], v[0], edge_to_midpoint);

    // Center node (average of all 4 vertices)
    cfem::Vector3D center = (m_vertices[v[0]].coords + m_vertices[v[1]].coords +
                       m_vertices[v[2]].coords + m_vertices[v[3]].coords) *
                      0.25;
    CfemInt n8 = addVertex(center.x, center.y, center.z);

    pQuadCell.vertex_ids = {v[0], v[1], v[2], v[3], n4, n5, n6, n7, n8};
    pQuadCell.type = CellType::QUAD9;
}

inline void Mesh::upgradeTetra4ToTetra10(Cell &pTetraCell,
                                         std::map<std::pair<CfemInt, CfemInt>, CfemInt> &edge_to_midpoint)
{
    std::vector<CfemInt> v = pTetraCell.vertex_ids;
    CfemInt n4 = getMidpoint(v[0], v[1], edge_to_midpoint);
    CfemInt n5 = getMidpoint(v[1], v[2], edge_to_midpoint);
    CfemInt n6 = getMidpoint(v[2], v[0], edge_to_midpoint);
    CfemInt n7 = getMidpoint(v[0], v[3], edge_to_midpoint);
    CfemInt n8 = getMidpoint(v[1], v[3], edge_to_midpoint);
    CfemInt n9 = getMidpoint(v[2], v[3], edge_to_midpoint);

    pTetraCell.vertex_ids = {v[0], v[1], v[2], v[3], n4, n5, n6, n7, n8, n9};
    pTetraCell.type = CellType::TET10;
}

inline void Mesh::upgradePrism6ToPrism18(Cell &pPrismCell,
                                         std::map<std::pair<CfemInt, CfemInt>, CfemInt> &edge_to_midpoint)
{
    std::vector<CfemInt> v = pPrismCell.vertex_ids;

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
    auto getFaceCenter = [&](const std::vector<CfemInt> &f_v)
    {
        cfem::Vector3D sum({0, 0, 0});
        for (auto id : f_v)
            sum = sum + m_vertices[id].coords;
        cfem::Vector3D mid = sum * 0.25;
        return addVertex(mid.x, mid.y, mid.z);
    };

    CfemInt n15 = getFaceCenter({v[0], v[1], v[4], v[3]}); // Face 0-1-4-3
    CfemInt n16 = getFaceCenter({v[1], v[2], v[5], v[4]}); // Face 1-2-5-4
    CfemInt n17 = getFaceCenter({v[2], v[0], v[3], v[5]}); // Face 2-0-3-5

    pPrismCell.vertex_ids = {
        v[0], v[1], v[2], v[3], v[4], v[5],      // Corners
        n6, n7, n8, n9, n10, n11, n12, n13, n14, // Mid-edges
        n15, n16, n17                            // Face centers
    };
    pPrismCell.type = CellType::PRISM18;
}

inline void Mesh::upgradeHex8ToHex27(Cell &pHexCell,
                                     std::map<std::pair<CfemInt, CfemInt>, CfemInt> &edge_to_midpoint)
{
    std::vector<CfemInt> v = pHexCell.vertex_ids;
    // Mid-edges (12 nodes)
    // Bottom edges
    CfemInt n8 = getMidpoint(v[0], v[1], edge_to_midpoint);
    CfemInt n9 = getMidpoint(v[1], v[2], edge_to_midpoint);
    CfemInt n10 = getMidpoint(v[2], v[3], edge_to_midpoint);
    CfemInt n11 = getMidpoint(v[3], v[0], edge_to_midpoint);
    // Top edges
    CfemInt n12 = getMidpoint(v[4], v[5], edge_to_midpoint);
    CfemInt n13 = getMidpoint(v[5], v[6], edge_to_midpoint);
    CfemInt n14 = getMidpoint(v[6], v[7], edge_to_midpoint);
    CfemInt n15 = getMidpoint(v[7], v[4], edge_to_midpoint);
    // Vertical edges
    CfemInt n16 = getMidpoint(v[0], v[4], edge_to_midpoint);
    CfemInt n17 = getMidpoint(v[1], v[5], edge_to_midpoint);
    CfemInt n18 = getMidpoint(v[2], v[6], edge_to_midpoint);
    CfemInt n19 = getMidpoint(v[3], v[7], edge_to_midpoint);

    // Face Centers (6 nodes)
    auto getFaceCenter = [&](const std::vector<CfemInt> &f_v)
    {
        cfem::Vector3D sum({0, 0, 0});
        for (auto id : f_v)
            sum = sum + m_vertices[id].coords;
        cfem::Vector3D mid = sum * 0.25;
        return addVertex(mid.x, mid.y, mid.z);
    };

    CfemInt n20 = getFaceCenter({v[0], v[4], v[7], v[3]}); // Bottom
    CfemInt n21 = getFaceCenter({v[5], v[1], v[2], v[6]}); // Top
    CfemInt n22 = getFaceCenter({v[0], v[4], v[5], v[1]}); // Front
    CfemInt n23 = getFaceCenter({v[3], v[2], v[6], v[7]}); // Right
    CfemInt n24 = getFaceCenter({v[0], v[1], v[2], v[3]}); // Back
    CfemInt n25 = getFaceCenter({v[4], v[5], v[6], v[7]}); // Left

    // Body Center (1 node)
    cfem::Vector3D bSum({0, 0, 0});
    for (auto id : v)
        bSum = bSum + m_vertices[id].coords;
    CfemInt n26 = addVertex(bSum.x * 0.125, bSum.y * 0.125, bSum.z * 0.125);

    pHexCell.vertex_ids = {
        v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7],           // 0-7: Corners
        n8, n9, n10, n11, n12, n13, n14, n15, n16, n17, n18, n19, // 8-19: Edges
        n20, n21, n22, n23, n24, n25,                             // 20-25: Faces
        n26                                                       // 26: Body
    };

    pHexCell.type = CellType::HEX27;
}

/**
 * @brief Retrieves the physical coordinates of all nodes belonging to an element.
 * @param elementId The index of the element.
 * @param localCoords [out] A vector that will be filled with the cfem::Vector3D coordinates.
 */
inline void Mesh::getElementCoordinates(CfemInt elementId, cfem::DynamicVector<cfem::Vector3D>& localCoords) const 
{
    const std::vector<int>& nodeIndices = m_cells[elementId].vertex_ids;

    const int n = nodeIndices.size();

    localCoords.resize(n);

    for (int i = 0; i < n; ++i) {
        int globalIdx = nodeIndices[i];
        localCoords[i] = m_vertices[globalIdx].coords;
    }
}

/**
 * @brief
 *
 * @param mesh
 * @param nx
 * @param ny
 * @param corners
 * @param subDiv
 * @param names
 */
inline void QuadMeshBuilder::build(Mesh &mesh, 
                                   MeshOrder order,
                                   const Subdivisions& subDiv,
                                   const Corners &corners,
                                   ElementType elemType,
                                   const BoundaryNames &names)
{
    mesh.clear();
    mesh.setDimension(2);
    mesh.setOrder(MeshOrder::Linear);
    int nodeRows = subDiv.Nx + 1;
    int nodeCols = subDiv.Ny + 1;
    std::vector<CfemInt> gridIds(nodeRows * nodeCols);
    // 1. Vertices
    for (int j = 0; j < nodeCols; ++j)
    {
        double v = static_cast<double>(j) / subDiv.Ny;
        for (int i = 0; i < nodeRows; ++i)
        {
            double u = static_cast<double>(i) / subDiv.Nx;
            cfem::Vector3D p = corners.P0 * ((1 - u) * (1 - v)) + corners.P1 * (u * (1 - v)) +
                         corners.P2 * (u * v) + corners.P3 * ((1 - u) * v);
            gridIds[j * nodeRows + i] = mesh.addVertex(p.x, p.y, p.z);
        }
    }
    // 2. Cells Generation (Fixed the break in FirstDiagonal)
    for (int j = 0; j < subDiv.Ny; ++j)
    {
        for (int i = 0; i < subDiv.Nx; ++i)
        {
            CfemInt v0 = gridIds[j * nodeRows + i];
            CfemInt v1 = gridIds[j * nodeRows + (i + 1)];
            CfemInt v2 = gridIds[(j + 1) * nodeRows + (i + 1)];
            CfemInt v3 = gridIds[(j + 1) * nodeRows + i];
            switch (elemType)
            {
            case ElementType::None:
                mesh.addCell(CellType::QUAD4, {v0, v1, v2, v3});
                break;
            case ElementType::FirstDiagonal:
                mesh.addCell(CellType::TRI3, {v0, v1, v2});
                mesh.addCell(CellType::TRI3, {v0, v2, v3});
                break;
            case ElementType::SecondDiagonal:
                mesh.addCell(CellType::TRI3, {v1, v2, v3});
                mesh.addCell(CellType::TRI3, {v1, v3, v0});
                break;
            case ElementType::CrissCross:
            {
                cfem::Vector3D midP = (mesh.getVertex(v0).coords + mesh.getVertex(v2).coords) * 0.5;
                CfemInt vMid = mesh.addVertex(midP.x, midP.y, midP.z);
                mesh.addCell(CellType::TRI3, {v0, v1, vMid});
                mesh.addCell(CellType::TRI3, {v1, v2, vMid});
                mesh.addCell(CellType::TRI3, {v2, v3, vMid});
                mesh.addCell(CellType::TRI3, {v3, v0, vMid});
                break;
            }
            }
        }
    }
    // 3. Finalize Topology
    mesh.buildMesh();
    
    if (order == MeshOrder::Quadratic)
    {
        mesh.upgradeToQuadratic();
    }
    // 4. Tag Boundaries
    applyBoundaryTags(mesh, subDiv.Nx, subDiv.Ny, gridIds, names);
}

/**
 * @brief
 *
 * @param mesh
 * @param nx
 * @param ny
 * @param grid
 * @param names
 */
inline void QuadMeshBuilder::applyBoundaryTags(Mesh &mesh, int nx, int ny,
                                               const std::vector<CfemInt> &grid,
                                               const BoundaryNames &names)
{
    // Create Mesh Entities (Boundaries)
    CfemInt bottomId = mesh.createEntity(names.bottom, 1); // Dim 1 for edges
    CfemInt rightId = mesh.createEntity(names.right, 1);
    CfemInt topId = mesh.createEntity(names.top, 1);
    CfemInt leftId = mesh.createEntity(names.left, 1);

    int nodeRows = nx + 1;

    // Iterate over all facets created by buildMesh()
    //for (CfemInt f_idx = 0; f_idx < mesh.getNumFacets(); ++f_idx)
    for (auto f_idx : mesh.getEntity(SYS_BOUNDARY_ID)->getFacetIds())
    {
        const Facet &facet = mesh.getFacet(f_idx);

        if (facet.isBoundary())
        {
            // Get the vertices of this facet
            CfemInt fv[2] = {facet.vertex_ids[0], facet.vertex_ids[1]};
            
            // Logical checks using coordinates or indices
            bool isBottom = true, isTop = true, isLeft = true, isRight = true;

            for (CfemInt vid : fv)
            {
                // Find i, j from the grid
                auto it = std::find(grid.begin(), grid.end(), vid);
                int idx = std::distance(grid.begin(), it);
                int i = idx % nodeRows;
                int j = idx / nodeRows;
                if (j != 0)
                    isBottom = false;
                if (j != ny)
                    isTop = false;
                if (i != 0)
                    isLeft = false;
                if (i != nx)
                    isRight = false;
            }

            if (isBottom) {
                mesh.getEntity(names.bottom)->addFacet(f_idx);
                for (auto v_id : facet.vertex_ids) mesh.getEntity(names.bottom)->addVertex(v_id);
            }
            else if (isTop){
                mesh.getEntity(names.top)->addFacet(f_idx);
                for (auto v_id : facet.vertex_ids) mesh.getEntity(names.top)->addVertex(v_id);
            }
            else if (isLeft){
                mesh.getEntity(names.left)->addFacet(f_idx);
                for (auto v_id : facet.vertex_ids) mesh.getEntity(names.left)->addVertex(v_id);
            }
            else if (isRight){
                mesh.getEntity(names.right)->addFacet(f_idx);
                for (auto v_id : facet.vertex_ids) mesh.getEntity(names.right)->addVertex(v_id);
            }
        }
    }
    mesh.getEntity(bottomId)->finalize();
    mesh.getEntity(topId)->finalize();
    mesh.getEntity(leftId)->finalize();
    mesh.getEntity(rightId)->finalize();
}

/**
 * @brief
 *
 * @param mesh
 * @param nx
 * @param ny
 * @param nz
 * @param corners
 * @param elemType
 */
inline void CubicMeshBuilder::build(Mesh& mesh, 
                                    MeshOrder order, 
                                    const Subdivisions& subDiv, 
                                    const Corners& corners,
                                    ElementType elemType)
{
    mesh.clear();
    mesh.setOrder(MeshOrder::Linear);
    mesh.setDimension(3);
    CfemInt nx = subDiv.Nx;
    CfemInt ny = subDiv.Ny;
    CfemInt nz = subDiv.Nz;

    auto gridIdx = [&](int i, int j, int k)    
    {
        return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;
    };

    // 1. Vertices (Trilinear Interpolation)
    for (int k = 0; k <= nz; ++k)
    {
        double w = (double)k / nz;

        for (int j = 0; j <= ny; ++j)
        {
            double v = (double)j / ny;

            for (int i = 0; i <= nx; ++i)
            {
                double u = (double)i / nx;
                cfem::Vector3D p =
                    corners.P0 * (1 - u) * (1 - v) * (1 - w) + corners.P1 * u * (1 - v) * (1 - w) +
                    corners.P2 * u * v * (1 - w) + corners.P3 * (1 - u) * v * (1 - w) +
                    corners.P4 * (1 - u) * (1 - v) * w + corners.P5 * u * (1 - v) * w +
                    corners.P6 * u * v * w + corners.P7 * (1 - u) * v * w;
                mesh.addVertex(p.x, p.y, p.z);
            }
        }
    }

    // 2. Cell Generation
    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {
                // Local indices for the current "Voxel"
                CfemInt v0 = gridIdx(i, j, k);
                CfemInt v1 = gridIdx(i + 1, j, k);
                CfemInt v2 = gridIdx(i + 1, j + 1, k);
                CfemInt v3 = gridIdx(i, j + 1, k);
                CfemInt v4 = gridIdx(i, j, k + 1);
                CfemInt v5 = gridIdx(i + 1, j, k + 1);
                CfemInt v6 = gridIdx(i + 1, j + 1, k + 1);
                CfemInt v7 = gridIdx(i, j + 1, k + 1);

                if (elemType == ElementType::None)
                {
                    mesh.addCell(CellType::HEX8, {v0, v1, v2, v3, v4, v5, v6, v7});
                }
                else if (elemType == ElementType::SixTet)
                {
                    mesh.addCell(CellType::TET4, {v0, v1, v3, v4});
                    mesh.addCell(CellType::TET4, {v1, v3, v4, v5});
                    mesh.addCell(CellType::TET4, {v4, v5, v3, v7});
                    mesh.addCell(CellType::TET4, {v1, v3, v5, v2});
                    mesh.addCell(CellType::TET4, {v3, v7, v5, v2}); 
                    mesh.addCell(CellType::TET4, {v5, v7, v6, v2});
                }
            }
        }
    }

    mesh.buildMesh();
    if (order == MeshOrder::Quadratic)
        mesh.upgradeToQuadratic();
}

inline void DiskMeshBuilder::build(Mesh &mesh,
                                   MeshOrder order,
                                   CfemReal target_mesh_size,
                                   const cfem::Vector3D &center,
                                   double radius,
                                   ElementType type)
{
    mesh.clear();
    mesh.setOrder(MeshOrder::Linear);

    CfemInt n = static_cast<int>(std::ceil((2.0 * radius) / target_mesh_size));
    if (n < 1)
        n = 1;
    int nx = n;
    int ny = n;

    std::vector<CfemInt> ids((nx + 1) * (ny + 1));

    for (int j = 0; j <= ny; ++j)
    {
        double eta = -1.0 + 2.0 * j / ny;

        for (int i = 0; i <= nx; ++i)
        {
            double xi = -1.0 + 2.0 * i / nx;

            double x = 0.0, y = 0.0;

            if (std::abs(xi) < cfem::ZERO_TOLERANCE && std::abs(eta) < cfem::ZERO_TOLERANCE)
            {
                x = 0.0;
                y = 0.0;
            }
            else if (std::abs(xi) > std::abs(eta))
            {
                double r = xi;
                double theta = (M_PI / 4.0) * (eta / xi);

                x = r * std::cos(theta);
                y = r * std::sin(theta);
            }
            else
            {
                double r = eta;
                double theta = (M_PI / 4.0) * (xi / eta);

                x = r * std::sin(theta);
                y = r * std::cos(theta);
            }

            cfem::Vector3D p{
                center.x + radius * x,
                center.y + radius * y,
                center.z};

            ids[j * (nx + 1) + i] =
                mesh.addVertex(p.x, p.y, p.z);
        }
    }

    // 2. Quads
    for (int j = 0; j < ny; ++j)
    {
        for (int i = 0; i < nx; ++i)
        {

            CfemInt v0 = ids[j * (nx + 1) + i];
            CfemInt v1 = ids[j * (nx + 1) + i + 1];
            CfemInt v2 = ids[(j + 1) * (nx + 1) + i + 1];
            CfemInt v3 = ids[(j + 1) * (nx + 1) + i];

            bool isCorner = ((i == 0 || i == nx - 1) && (j == 0 || j == ny - 1));

            if (type == ElementType::Quadrilateral)
            {
                mesh.addCell(CellType::QUAD4, {v0, v1, v2, v3});
            }
            else if (type == ElementType::Triangle)
            {

                if (!isCorner)
                {
                    cfem::Vector3D midP = (mesh.getVertex(v0).coords + mesh.getVertex(v1).coords +
                                     mesh.getVertex(v2).coords + mesh.getVertex(v3).coords) *
                                    0.25;
                    CfemInt vMid = mesh.addVertex(midP.x, midP.y, midP.z);

                    mesh.addCell(CellType::TRI3, {v0, v1, vMid});
                    mesh.addCell(CellType::TRI3, {v1, v2, vMid});
                    mesh.addCell(CellType::TRI3, {v2, v3, vMid});
                    mesh.addCell(CellType::TRI3, {v3, v0, vMid});
                }
                else
                {
                    if ((i == 0 && j == 0) || (i == nx - 1 && j == ny - 1))
                    {
                        mesh.addCell(CellType::TRI3, {v0, v1, v2});
                        mesh.addCell(CellType::TRI3, {v0, v2, v3});
                    }
                    else if ((i == nx - 1 && j == 0) || (i == 0 && j == ny - 1))
                    {
                        mesh.addCell(CellType::TRI3, {v0, v1, v3});
                        mesh.addCell(CellType::TRI3, {v1, v2, v3});
                    }
                }
            }
        }
    }

    mesh.buildMesh();

    mesh.relaxMesh(5);

    if (order == MeshOrder::Quadratic)
    {
        mesh.upgradeToQuadratic();
        projectBoundaryMidNodes(mesh, center, radius);
    }
}

inline void DiskMeshBuilder::projectBoundaryMidNodes(Mesh &mesh, const cfem::Vector3D &center, double radius)
{
    MeshEntity *boundary = mesh.getEntity(cfem::SYS_BOUNDARY_ID);
    if (!boundary)
        return;

    const auto &boundaryVids = boundary->getVertexIds();

    for (CfemInt vid : boundaryVids)
    {
        Vertex &v = mesh.getVertex(vid);

        cfem::Vector3D vec = v.coords - center;
        double currentDist = vec.norm();

        if (currentDist > cfem::ZERO_TOLERANCE)
        {
            v.coords = center + (vec / currentDist) * radius;
        }
    }
}

inline void SphereMeshBuilder::build(Mesh &mesh,
                                     MeshOrder order,
                                     CfemReal h,
                                     const cfem::Vector3D &center,
                                     double radius,
                                     ElementType type)
{
    mesh.clear();
    mesh.setOrder(MeshOrder::Linear);
    mesh.setDimension(3);

    CfemInt n = static_cast<CfemInt>(std::ceil((2.0 * radius) / h));
    if (n < 1)
        n = 1;

    int nx = n;
    int ny = n;
    int nz = n;

    std::vector<CfemInt> ids((nx + 1) * (ny + 1) * (nz + 1));

    auto idx = [&](int i, int j, int k)
    {
        return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;
    };

    /*----------------------------
      1. Cube → Ball mapping
    -----------------------------*/
    for (int k = 0; k <= nz; ++k)
    {
        double zeta = -1.0 + 2.0 * k / nz;

        for (int j = 0; j <= ny; ++j)
        {
            double eta = -1.0 + 2.0 * j / ny;

            for (int i = 0; i <= nx; ++i)
            {
                double xi = -1.0 + 2.0 * i / nx;

                double x = xi * std::sqrt(
                                    1.0 - 0.5 * eta * eta - 0.5 * zeta * zeta + eta * eta * zeta * zeta / 3.0);

                double y = eta * std::sqrt(
                                     1.0 - 0.5 * xi * xi - 0.5 * zeta * zeta + xi * xi * zeta * zeta / 3.0);

                double z = zeta * std::sqrt(
                                      1.0 - 0.5 * xi * xi - 0.5 * eta * eta + xi * xi * eta * eta / 3.0);

                cfem::Vector3D p{
                    center.x + radius * x,
                    center.y + radius * y,
                    center.z + radius * z};

                ids[idx(i, j, k)] =
                    mesh.addVertex(p.x, p.y, p.z);
            }
        }
    }

    /*----------------------------
      2. Cells
    -----------------------------*/
    for (int k = 0; k < nz; ++k)
    {
        for (int j = 0; j < ny; ++j)
        {
            for (int i = 0; i < nx; ++i)
            {

                CfemInt v0 = ids[idx(i, j, k)];
                CfemInt v1 = ids[idx(i + 1, j, k)];
                CfemInt v2 = ids[idx(i + 1, j + 1, k)];
                CfemInt v3 = ids[idx(i, j + 1, k)];

                CfemInt v4 = ids[idx(i, j, k + 1)];
                CfemInt v5 = ids[idx(i + 1, j, k + 1)];
                CfemInt v6 = ids[idx(i + 1, j + 1, k + 1)];
                CfemInt v7 = ids[idx(i, j + 1, k + 1)];

                if (type == ElementType::Hexahedron)
                {

                    mesh.addCell(CellType::HEX8,
                                 {v0, v1, v2, v3, v4, v5, v6, v7});
                }
                else
                {
                    cfem::Vector3D midP = (mesh.getVertex(v0).coords + mesh.getVertex(v1).coords +
                                     mesh.getVertex(v2).coords + mesh.getVertex(v3).coords + 
                                     mesh.getVertex(v4).coords + mesh.getVertex(v5).coords +
                                     mesh.getVertex(v6).coords + mesh.getVertex(v7).coords) *
                                    0.125;
                    CfemInt vMid = mesh.addVertex(midP.x, midP.y, midP.z);
                    // 6-tet decomposition
                    mesh.addCell(CellType::TET4, {v0, v1, v3, v4});
                    mesh.addCell(CellType::TET4, {v1, v3, v4, v5});
                    mesh.addCell(CellType::TET4, {v4, v5, v3, v7});
                    mesh.addCell(CellType::TET4, {v1, v3, v5, v2});
                    mesh.addCell(CellType::TET4, {v3, v7, v5, v2}); 
                    mesh.addCell(CellType::TET4, {v5, v7, v6, v2});

                }
            }
        }
    }

    mesh.buildMesh();

    mesh.relaxMesh(10);

    if (order == MeshOrder::Quadratic)
    {
        mesh.upgradeToQuadratic();
        projectBoundaryNodes(mesh, center, radius);
    }
}

inline void SphereMeshBuilder::projectBoundaryNodes(
    Mesh &mesh,
    const cfem::Vector3D &center,
    double radius)
{
    MeshEntity *boundary = mesh.getEntity("Boundary");
    if (!boundary)
        return;

    for (CfemInt vid : boundary->getVertexIds())
    {

        Vertex &v = mesh.getVertex(vid);

        cfem::Vector3D d = v.coords - center;
        double r = d.norm();

        if (r > cfem::ZERO_TOLERANCE)
        {

            v.coords = center + d * (radius / r);
        }
    }
}


inline void Mesh::relaxMesh(int iterations,
                           double lambda)
{
    const size_t numVertices = getNumVertices();
    const size_t numCells    = getNumCells();

    if (numVertices == 0 || numCells == 0)
        return;

    // --------------------------------------------------
    // 0. Make sure topology exists
    // --------------------------------------------------
#ifdef NDEBUG
    CFEM_ASSERT(!m_cells.empty());
#endif

    // --------------------------------------------------
    // 1. Identify boundary vertices (O(1))
    // --------------------------------------------------
    std::vector<bool> isBoundary(numVertices, false);

    MeshEntity* boundary = getEntity(cfem::SYS_BOUNDARY_ID);
    if (boundary)
    {
        for (CfemInt id : boundary->getVertexIds())
        {
            if (id >= 0 && (size_t)id < numVertices)
                isBoundary[id] = true;
        }
    }

    // --------------------------------------------------
    // 2. Build vertex adjacency (once)
    // --------------------------------------------------
    std::vector<std::set<CfemInt>> adj(numVertices);

    for (size_t c = 0; c < numCells; ++c)
    {
        const Cell& cell = getCell(c);
        const auto& v    = cell.vertex_ids;
        const size_t n   = v.size();

        for (size_t i = 0; i < n; ++i)
        {
            CfemInt v0 = v[i];
            CfemInt v1 = v[(i + 1) % n];

            if (v0 >= 0 && v1 >= 0)
            {
                adj[v0].insert(v1);
                adj[v1].insert(v0);
            }
        }
    }

    // --------------------------------------------------
    // 3. Flatten adjacency (cache-friendly)
    // --------------------------------------------------
    std::vector<std::vector<CfemInt>> neighbors(numVertices);

    for (size_t i = 0; i < numVertices; ++i)
    {
        neighbors[i].assign(adj[i].begin(), adj[i].end());
    }

    adj.clear();
    adj.shrink_to_fit();

    // --------------------------------------------------
    // 4. Buffers
    // --------------------------------------------------
    std::vector<cfem::Vector3D> nextCoords(numVertices);

    // Clamp lambda
    if (lambda <= 0.0) lambda = 0.3;
    if (lambda >= 1.0) lambda = 0.8;

    // --------------------------------------------------
    // 5. Iterations
    // --------------------------------------------------
    for (int iter = 0; iter < iterations; ++iter)
    {

#pragma omp parallel for schedule(static)
        for (long i = 0; i < (long)numVertices; ++i)
        {
            // Boundary / isolated nodes
            if (isBoundary[i] || neighbors[i].empty())
            {
                nextCoords[i] = getVertex(i).coords;
                continue;
            }

            // Laplacian average
            cfem::Vector3D sum{0, 0, 0};

            for (CfemInt j : neighbors[i])
            {
                sum = sum + getVertex(j).coords;
            }

            cfem::Vector3D avg = sum * (1.0 / neighbors[i].size());
            cfem::Vector3D old = getVertex(i).coords;

            // Damped update (Taubin-style lambda step)
            nextCoords[i] = old + lambda * (avg - old);
        }

        // --------------------------------------------------
        // 6. Commit update
        // --------------------------------------------------
        for (size_t i = 0; i < numVertices; ++i)
        {
            if (!isBoundary[i])
            {
                getVertex(i).coords = nextCoords[i];
            }
        }
    }
}

}//namespace