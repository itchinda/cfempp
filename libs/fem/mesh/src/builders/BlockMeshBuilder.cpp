#include "BlockMeshBuilder.h"
#include "omp.h"

namespace cfem
{

    std::shared_ptr<Mesh> BlockMeshBuilder::create(
                                                    MeshOrder order,
                                                    const BlockSubdivisions &subDiv,
                                                    const BlockCorners &corners,
                                                    BlockElementType elemType,
                                                    const BlockBoundaryNames &names
                                                )
    {
        auto mesh = Mesh::create();
        BlockMeshBuilder builder(ConstructorToken{});
        builder.build(*mesh, order, subDiv, corners, elemType, names);
        return mesh;
    }

    /**
     * @details
     * The mesh construction proceeds in three stages:
     *
     * 1. Vertex generation:
     *    Vertices are created using trilinear interpolation of the
     *    eight block corner coordinates.
     *
     * 2. Cell generation:
     *    Each structured voxel is subdivided into either:
     *      - One HEX8 cell
     *      - Six TET4 cells
     *
     * 3. Post-processing:
     *    - Mesh topology is finalized via buildMesh().
     *    - If requested, the mesh is upgraded to quadratic order.
     *    - Boundary entities are tagged.
     *
     * The resulting mesh is fully connected and ready for FE space construction.
     */
    inline void BlockMeshBuilder::build(
                                            Mesh &mesh,
                                            MeshOrder order,
                                            const BlockSubdivisions &subDiv,
                                            const BlockCorners &corners,
                                            BlockElementType elemType,
                                            const BlockBoundaryNames &names
                                        )
    {
        mesh.clear();
        mesh.setOrder(MeshOrder::Linear);
        mesh.setDimension(3);
        CfemInt nx = subDiv.nx;
        CfemInt ny = subDiv.ny;
        CfemInt nz = subDiv.nz;
        CfemInt totalVertices = (nx + 1) * (ny + 1) * (nz + 1);

        mesh.m_vertices.reserve(totalVertices);
        std::vector<CfemInt> gridIds(totalVertices);

        int cellMultiplier = (elemType == BlockElementType::Hexahedron) ? 1 : 6;
        if (elemType == BlockElementType::Edge) cellMultiplier = 2;

        mesh.m_cells.reserve(nx * ny * nz * cellMultiplier);

        auto gridIdx = [&](CfemInt i, CfemInt j, CfemInt k)
        {
            return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;
        };

        // 1. Vertices (Trilinear Interpolation)
        for (CfemInt k = 0; k <= nz; ++k)
        {
            CfemReal w = (CfemReal)k / nz;

            for (CfemInt j = 0; j <= ny; ++j)
            {
                CfemReal v = (CfemReal)j / ny;

                for (CfemInt i = 0; i <= nx; ++i)
                {
                    CfemReal u = (CfemReal)i / nx;
                    cfem::Vector3D p =
                        corners.vertices[0] * (1 - u) * (1 - v) * (1 - w) + corners.vertices[1] * u * (1 - v) * (1 - w) +
                        corners.vertices[2] * u * v * (1 - w) + corners.vertices[3] * (1 - u) * v * (1 - w) +
                        corners.vertices[4] * (1 - u) * (1 - v) * w + corners.vertices[5] * u * (1 - v) * w +
                        corners.vertices[6] * u * v * w + corners.vertices[7] * (1 - u) * v * w;
                    gridIds[gridIdx(i, j, k)] = mesh.addVertex(p.x, p.y, p.z);
                }
            }
        }

        // 2. Cell Generation
        for (CfemInt k = 0; k < nz; ++k)
        {
            for (CfemInt j = 0; j < ny; ++j)
            {
                for (CfemInt i = 0; i < nx; ++i)
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

                    if (elemType == BlockElementType::Hexahedron)
                    {
                        mesh.addCell(CellType::HEX8, {v0, v1, v2, v3, v4, v5, v6, v7});
                    }
                    else if (elemType == BlockElementType::Edge)
                    {
                        mesh.addCell(CellType::PRISM6, {v0, v1, v3, v4, v5, v7});
                        mesh.addCell(CellType::PRISM6, {v1, v2, v3, v5, v6, v7});
                    }
                    else if (elemType == BlockElementType::Tetrahedron)
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

        applyBoundaryTags(mesh, subDiv.nx, subDiv.ny, subDiv.nz, gridIds, names);
    }

    /**
     * @details
     * Boundary detection is performed using structured indexing.
     *
     * Each boundary facet is classified by checking whether all its
     * vertices lie on one of the planes:
     *
     *   i = 0        -> xmin
     *   i = nx       -> xmax
     *   j = 0        -> ymin
     *   j = ny       -> ymax
     *   k = 0        -> zmin
     *   k = nz       -> zmax
     *
     * Facets are then assigned to the corresponding mesh entity.
     *
     * Duplicate vertices are removed during entity finalization.
     */
    void BlockMeshBuilder::applyBoundaryTags(
                                                Mesh &mesh,
                                                CfemInt nx,
                                                CfemInt ny,
                                                CfemInt nz,
                                                const std::vector<CfemInt> &grid,
                                                const BlockBoundaryNames &names
                                            )
    {
        // Create/Retrieve Mesh Entities for each boundary
        CfemInt xminId = mesh.createEntity(names.xmin, 2);
        CfemInt xmaxId = mesh.createEntity(names.xmax, 2);
        CfemInt yminId = mesh.createEntity(names.ymin, 2);
        CfemInt ymaxId = mesh.createEntity(names.ymax, 2);
        CfemInt zminId = mesh.createEntity(names.zmin, 2);
        CfemInt zmaxId = mesh.createEntity(names.zmax, 2);

        // Compute strides for structured indexing
        int nx1 = nx + 1;
        int ny1 = ny + 1;

        // Build inverse map: vertex_id -> linear index in structured grid
        std::vector<int> invGrid(mesh.getNumVertices(), -1);
        for (CfemInt k = 0; k < grid.size(); ++k)
            invGrid[grid[k]] = k;

        // Retrieve system boundary facets
        auto *sysBoundary = mesh.getEntity(cfem::SYS_BOUNDARY_ID);
        if (!sysBoundary)
            return;

        for (auto f_idx : sysBoundary->getFacetIds())
        {
            const Facet &facet = mesh.getFacet(f_idx);
            auto fv = mesh.getFacetConnectivity(f_idx);

            bool onXmin = true, onXmax = true;
            bool onYmin = true, onYmax = true;
            bool onZmin = true, onZmax = true;

            for (CfemInt vid : fv)
            {
                int idx = invGrid[vid];
                if (idx == -1)
                    continue; // Should not happen

                int i = idx % nx1;
                int j = (idx / nx1) % ny1;
                int k = idx / (nx1 * ny1);

                if (i != 0)
                    onXmin = false;
                if (i != nx)
                    onXmax = false;
                if (j != 0)
                    onYmin = false;
                if (j != ny)
                    onYmax = false;
                if (k != 0)
                    onZmin = false;
                if (k != nz)
                    onZmax = false;
            }

            // Lambda to assign facet to entity
            auto assignTo = [&](CfemInt entityId)
            {
                auto *entity = mesh.getEntity(entityId);
                entity->addFacet(f_idx);
                for (auto v_id : fv)
                    entity->addVertex(v_id);
            };

            if (onXmin)
                assignTo(xminId);
            else if (onXmax)
                assignTo(xmaxId);
            else if (onYmin)
                assignTo(yminId);
            else if (onYmax)
                assignTo(ymaxId);
            else if (onZmin)
                assignTo(zminId);
            else if (onZmax)
                assignTo(zmaxId);
        }

        // Finalize all entities to remove duplicate vertices
        mesh.getEntity(xminId)->finalize();
        mesh.getEntity(xmaxId)->finalize();
        mesh.getEntity(yminId)->finalize();
        mesh.getEntity(ymaxId)->finalize();
        mesh.getEntity(zminId)->finalize();
        mesh.getEntity(zmaxId)->finalize();

    }




inline void BlockMeshBuilder::build_omp(
    Mesh &mesh,
    MeshOrder order,
    const BlockSubdivisions &subDiv,
    const BlockCorners &corners,
    BlockElementType elemType,
    const BlockBoundaryNames &names)
{
    mesh.clear();
    mesh.setOrder(MeshOrder::Linear);
    mesh.setDimension(3);

    CfemInt nx = subDiv.nx;
    CfemInt ny = subDiv.ny;
    CfemInt nz = subDiv.nz;
    CfemInt totalVertices = (nx + 1) * (ny + 1) * (nz + 1);

    // -------------------------------
    // 1. Vertices
    // -------------------------------
    std::vector<cfem::Vector3D> vertexCoords(totalVertices);

#pragma omp parallel for collapse(3) schedule(static)
    for (CfemInt k = 0; k <= nz; ++k)
        for (CfemInt j = 0; j <= ny; ++j)
            for (CfemInt i = 0; i <= nx; ++i)
            {
                CfemReal u = (CfemReal)i / nx;
                CfemReal v = (CfemReal)j / ny;
                CfemReal w = (CfemReal)k / nz;
                vertexCoords[k * (nx + 1) * (ny + 1) + j * (nx + 1) + i] =
                    corners.vertices[0] * (1 - u) * (1 - v) * (1 - w) +
                    corners.vertices[1] * u * (1 - v) * (1 - w) +
                    corners.vertices[2] * u * v * (1 - w) +
                    corners.vertices[3] * (1 - u) * v * (1 - w) +
                    corners.vertices[4] * (1 - u) * (1 - v) * w +
                    corners.vertices[5] * u * (1 - v) * w +
                    corners.vertices[6] * u * v * w +
                    corners.vertices[7] * (1 - u) * v * w;
            }

    std::vector<CfemInt> gridIds(totalVertices);
    for (CfemInt idx = 0; idx < totalVertices; ++idx)
        gridIds[idx] = mesh.addVertex(vertexCoords[idx].x, vertexCoords[idx].y, vertexCoords[idx].z);

    // -------------------------------
    // 2. Cells
    // -------------------------------
    int cellMultiplier = (elemType == BlockElementType::Hexahedron) ? 1 :
                         (elemType == BlockElementType::Edge) ? 2 : 6;

    CfemInt totalCells = nx * ny * nz * cellMultiplier;
    std::vector<std::vector<CfemInt>> cells(totalCells);

    auto gridIdx = [&](CfemInt i, CfemInt j, CfemInt k) {
        return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;
    };

#pragma omp parallel for collapse(3) schedule(static)
    for (CfemInt k = 0; k < nz; ++k)
        for (CfemInt j = 0; j < ny; ++j)
            for (CfemInt i = 0; i < nx; ++i)
            {
                int baseId = (k * ny * nx + j * nx + i) * cellMultiplier;
                CfemInt v0 = gridIdx(i, j, k);
                CfemInt v1 = gridIdx(i + 1, j, k);
                CfemInt v2 = gridIdx(i + 1, j + 1, k);
                CfemInt v3 = gridIdx(i, j + 1, k);
                CfemInt v4 = gridIdx(i, j, k + 1);
                CfemInt v5 = gridIdx(i + 1, j, k + 1);
                CfemInt v6 = gridIdx(i + 1, j + 1, k + 1);
                CfemInt v7 = gridIdx(i, j + 1, k + 1);

                if (elemType == BlockElementType::Hexahedron)
                    cells[baseId] = {v0, v1, v2, v3, v4, v5, v6, v7};
                else if (elemType == BlockElementType::Edge)
                {
                    cells[baseId + 0] = {v0, v1, v3, v4, v5, v7};
                    cells[baseId + 1] = {v1, v2, v3, v5, v6, v7};
                }
                else if (elemType == BlockElementType::Tetrahedron)
                {
                    cells[baseId + 0] = {v0, v1, v3, v4};
                    cells[baseId + 1] = {v1, v3, v4, v5};
                    cells[baseId + 2] = {v4, v5, v3, v7};
                    cells[baseId + 3] = {v1, v3, v5, v2};
                    cells[baseId + 4] = {v3, v7, v5, v2};
                    cells[baseId + 5] = {v5, v7, v6, v2};
                }
            }

    for (auto &c : cells)
        mesh.addCell(c.size() == 4 ? CellType::TET4 :
                     c.size() == 6 ? CellType::PRISM6 : CellType::HEX8, c);

    mesh.buildMesh();
    if (order == MeshOrder::Quadratic)
        mesh.upgradeToQuadratic();

    // -------------------------------
    // 3. Boundary tagging (thread-local)
    // -------------------------------
    applyBoundaryTags_omp(mesh, nx, ny, nz, gridIds, names);
}

// -------------------------------
// Thread-safe boundary tagging
// -------------------------------
void BlockMeshBuilder::applyBoundaryTags_omp(
    Mesh &mesh,
    CfemInt nx, CfemInt ny, CfemInt nz,
    const std::vector<CfemInt> &grid,
    const BlockBoundaryNames &names)
{
    CfemInt nx1 = nx + 1, ny1 = ny + 1;
    std::vector<int> invGrid(mesh.getNumVertices(), -1);
    for (CfemInt k = 0; k < grid.size(); ++k)
        invGrid[grid[k]] = k;

    CfemInt xminId = mesh.createEntity(names.xmin, 2);
    CfemInt xmaxId = mesh.createEntity(names.xmax, 2);
    CfemInt yminId = mesh.createEntity(names.ymin, 2);
    CfemInt ymaxId = mesh.createEntity(names.ymax, 2);
    CfemInt zminId = mesh.createEntity(names.zmin, 2);
    CfemInt zmaxId = mesh.createEntity(names.zmax, 2);

    auto *sysBoundary = mesh.getEntity(cfem::SYS_BOUNDARY_ID);
    if (!sysBoundary)
        return;

    // Thread-local assignments
#pragma omp parallel
    {
        std::vector<CfemInt> xminFacets, xmaxFacets, yminFacets, ymaxFacets, zminFacets, zmaxFacets;
        std::vector<std::vector<CfemInt>> xminVerts, xmaxVerts, yminVerts, ymaxVerts, zminVerts, zmaxVerts;

#pragma omp for schedule(static)
        for (size_t f_idx = 0; f_idx < sysBoundary->getFacetIds().size(); ++f_idx)
        {
            CfemInt idx = sysBoundary->getFacetIds()[f_idx];
            const Facet &facet = mesh.getFacet(idx);
            auto fv = mesh.getFacetConnectivity(idx);

            bool onXmin = true, onXmax = true, onYmin = true, onYmax = true, onZmin = true, onZmax = true;
            for (CfemInt vid : fv)
            {
                int vIdx = invGrid[vid];
                if (vIdx == -1) continue;
                int i = vIdx % nx1;
                int j = (vIdx / nx1) % ny1;
                int k = vIdx / (nx1 * ny1);

                if (i != 0) onXmin = false;
                if (i != nx) onXmax = false;
                if (j != 0) onYmin = false;
                if (j != ny) onYmax = false;
                if (k != 0) onZmin = false;
                if (k != nz) onZmax = false;
            }

            auto assignTo = [&](std::vector<CfemInt> &facetList,
                    std::vector<std::vector<CfemInt>> &vertsList)
                {
                    facetList.push_back(idx);
                    vertsList.emplace_back(fv.begin(), fv.end());
                };

            if (onXmin) assignTo(xminFacets, xminVerts);
            else if (onXmax) assignTo(xmaxFacets, xmaxVerts);
            else if (onYmin) assignTo(yminFacets, yminVerts);
            else if (onYmax) assignTo(ymaxFacets, ymaxVerts);
            else if (onZmin) assignTo(zminFacets, zminVerts);
            else if (onZmax) assignTo(zmaxFacets, zmaxVerts);
        }

#pragma omp critical
        {
            auto addEntityData = [&](CfemInt entityId, const std::vector<CfemInt> &facets, const std::vector<std::vector<CfemInt>> &verts)
            {
                auto *entity = mesh.getEntity(entityId);
                for (size_t k = 0; k < facets.size(); ++k)
                {
                    entity->addFacet(facets[k]);
                    for (CfemInt vid : verts[k])
                        entity->addVertex(vid);
                }
            };

            addEntityData(xminId, xminFacets, xminVerts);
            addEntityData(xmaxId, xmaxFacets, xmaxVerts);
            addEntityData(yminId, yminFacets, yminVerts);
            addEntityData(ymaxId, ymaxFacets, ymaxVerts);
            addEntityData(zminId, zminFacets, zminVerts);
            addEntityData(zmaxId, zmaxFacets, zmaxVerts);
        }
    }

    // Finalize all entities
    for (auto eid : {xminId, xmaxId, yminId, ymaxId, zminId, zmaxId})
        mesh.getEntity(eid)->finalize();
}

} // namespace cfem