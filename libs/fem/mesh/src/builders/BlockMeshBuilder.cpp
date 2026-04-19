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
        builder.build_omp(*mesh, order, subDiv, corners, elemType, names);
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

        // int cellMultiplier = (elemType == BlockElementType::Hexahedron) ? 1 : 6;
        int cellMultiplier = (elemType == BlockElementType::Hexahedron) ? 1 : 5;
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
                        // using six tets

                        // mesh.addCell(CellType::TET4, {v0, v1, v7, v4});
                        // mesh.addCell(CellType::TET4, {v0, v1, v2, v7});
                        // mesh.addCell(CellType::TET4, {v0, v2, v3, v7});
                        // mesh.addCell(CellType::TET4, {v1, v6, v4, v5});
                        // mesh.addCell(CellType::TET4, {v1, v2, v7, v6});
                        // mesh.addCell(CellType::TET4, {v1, v6, v7, v4});
                        
                        // using five tets
                        bool isEven = (i + j + k) % 2 == 0;
                        // if (isEven) {
                        //     mesh.addCell(CellType::TET4, {v0, v1, v3, v4});
                        //     mesh.addCell(CellType::TET4, {v1, v2, v3, v6});
                        //     mesh.addCell(CellType::TET4, {v4, v5, v6, v1});
                        //     mesh.addCell(CellType::TET4, {v4, v6, v7, v3});
                        //     mesh.addCell(CellType::TET4, {v1, v3, v4, v6}); // Central
                        // } else {
                        //     mesh.addCell(CellType::TET4, {v1, v0, v2, v5});
                        //     mesh.addCell(CellType::TET4, {v3, v0, v2, v7});
                        //     mesh.addCell(CellType::TET4, {v5, v4, v7, v0});
                        //     mesh.addCell(CellType::TET4, {v7, v6, v5, v2});
                        //     mesh.addCell(CellType::TET4, {v0, v2, v5, v7}); // Central
                        // }

                        if ((i + j + k) % 2 == 0) {
                            // Version A (Toutes les orientations corrigées pour detJ > 0)
                            mesh.addCell(CellType::TET4, {v0, v1, v3, v4});
                            mesh.addCell(CellType::TET4, {v1, v2, v3, v6});
                            mesh.addCell(CellType::TET4, {v1, v5, v6, v4}); // Corrigé (était v4, v5, v6, v1)
                            mesh.addCell(CellType::TET4, {v4, v7, v6, v3}); // Corrigé (était v4, v6, v7, v3)
                            mesh.addCell(CellType::TET4, {v1, v3, v4, v6}); // Central : Déjà bon
                        } else {
                            // Version B (Toutes les orientations corrigées pour detJ > 0)
                            mesh.addCell(CellType::TET4, {v0, v1, v2, v5}); // Corrigé (était v1, v0, v2, v5)
                            mesh.addCell(CellType::TET4, {v0, v2, v3, v7}); // Corrigé (était v3, v0, v2, v7)
                            mesh.addCell(CellType::TET4, {v4, v7, v5, v0}); // Corrigé (était v5, v4, v7, v0)
                            mesh.addCell(CellType::TET4, {v7, v6, v5, v2}); // Déjà bon
                            mesh.addCell(CellType::TET4, {v0, v5, v2, v7}); // Central : Corrigé (était v0, v2, v5, v7)
                        }

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
    // 1. Vertices (Calcul en parallèle, insertion séquentielle rapide)
    // -------------------------------
    std::vector<cfem::Vector3D> vertexCoords(totalVertices);

    #pragma omp parallel for collapse(3) schedule(static)
    for (CfemInt k = 0; k <= nz; ++k) {
        for (CfemInt j = 0; j <= ny; ++j) {
            for (CfemInt i = 0; i <= nx; ++i) {
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
        }
    }

    std::vector<CfemInt> gridIds(totalVertices);
    // Insertion séquentielle ultra-rapide (la mémoire de m_vertices a été pré-allouée)
    for (CfemInt idx = 0; idx < totalVertices; ++idx) {
        gridIds[idx] = mesh.addVertex(vertexCoords[idx].x, vertexCoords[idx].y, vertexCoords[idx].z);
    }

    // -------------------------------
    // 2. Cells (Data-Oriented Design)
    // -------------------------------
    int cellMultiplier = (elemType == BlockElementType::Hexahedron) ? 1 :
                         (elemType == BlockElementType::Edge) ? 2 : 5;
    int nodesPerCell = (elemType == BlockElementType::Hexahedron) ? 8 :
                       (elemType == BlockElementType::Edge) ? 6 : 4;

    CfemInt totalCells = nx * ny * nz * cellMultiplier;
    
    // LE SECRET EST ICI : Un seul bloc de mémoire plat contigu
    std::vector<CfemInt> bigCellBlock(totalCells * nodesPerCell);

    auto gridIdx = [&](CfemInt i, CfemInt j, CfemInt k) {
        return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;
    };

    #pragma omp parallel for collapse(3) schedule(static)
    for (CfemInt k = 0; k < nz; ++k) {
        for (CfemInt j = 0; j < ny; ++j) {
            for (CfemInt i = 0; i < nx; ++i) {
                int baseCellId = (k * ny * nx + j * nx + i) * cellMultiplier;
                int bIdx = baseCellId * nodesPerCell;

                CfemInt v0 = gridIdx(i, j, k);
                CfemInt v1 = gridIdx(i + 1, j, k);
                CfemInt v2 = gridIdx(i + 1, j + 1, k);
                CfemInt v3 = gridIdx(i, j + 1, k);
                CfemInt v4 = gridIdx(i, j, k + 1);
                CfemInt v5 = gridIdx(i + 1, j, k + 1);
                CfemInt v6 = gridIdx(i + 1, j + 1, k + 1);
                CfemInt v7 = gridIdx(i, j + 1, k + 1);

                if (elemType == BlockElementType::Hexahedron) {
                    bigCellBlock[bIdx+0]=v0; bigCellBlock[bIdx+1]=v1; bigCellBlock[bIdx+2]=v2; bigCellBlock[bIdx+3]=v3;
                    bigCellBlock[bIdx+4]=v4; bigCellBlock[bIdx+5]=v5; bigCellBlock[bIdx+6]=v6; bigCellBlock[bIdx+7]=v7;
                }
                else if (elemType == BlockElementType::Edge) {
                    bigCellBlock[bIdx+0]=v0; bigCellBlock[bIdx+1]=v1; bigCellBlock[bIdx+2]=v3; bigCellBlock[bIdx+3]=v4; bigCellBlock[bIdx+4]=v5; bigCellBlock[bIdx+5]=v7;
                    bIdx += 6;
                    bigCellBlock[bIdx+0]=v1; bigCellBlock[bIdx+1]=v2; bigCellBlock[bIdx+2]=v3; bigCellBlock[bIdx+3]=v5; bigCellBlock[bIdx+4]=v6; bigCellBlock[bIdx+5]=v7;
                }
                else if (elemType == BlockElementType::Tetrahedron) {
                    // using six tets
                    // bigCellBlock[bIdx+0]=v0; bigCellBlock[bIdx+1]=v1; bigCellBlock[bIdx+2]=v7; bigCellBlock[bIdx+3]=v4; bIdx+=4;
                    // bigCellBlock[bIdx+0]=v0; bigCellBlock[bIdx+1]=v1; bigCellBlock[bIdx+2]=v2; bigCellBlock[bIdx+3]=v7; bIdx+=4;
                    // bigCellBlock[bIdx+0]=v0; bigCellBlock[bIdx+1]=v2; bigCellBlock[bIdx+2]=v3; bigCellBlock[bIdx+3]=v7; bIdx+=4;
                    // bigCellBlock[bIdx+0]=v1; bigCellBlock[bIdx+1]=v6; bigCellBlock[bIdx+2]=v4; bigCellBlock[bIdx+3]=v5; bIdx+=4;
                    // bigCellBlock[bIdx+0]=v1; bigCellBlock[bIdx+1]=v2; bigCellBlock[bIdx+2]=v7; bigCellBlock[bIdx+3]=v6; bIdx+=4;
                    // bigCellBlock[bIdx+0]=v1; bigCellBlock[bIdx+1]=v6; bigCellBlock[bIdx+2]=v7; bigCellBlock[bIdx+3]=v4;

                    // using 5 tets
                    bool isEven = (i + j + k) % 2 == 0;
                    if (isEven) {
                        bigCellBlock[bIdx+0]=v0; bigCellBlock[bIdx+1]=v1; bigCellBlock[bIdx+2]=v3; bigCellBlock[bIdx+3]=v4; bIdx+=4;
                        bigCellBlock[bIdx+0]=v1; bigCellBlock[bIdx+1]=v2; bigCellBlock[bIdx+2]=v3; bigCellBlock[bIdx+3]=v6; bIdx+=4;
                        bigCellBlock[bIdx+0]=v1; bigCellBlock[bIdx+1]=v5; bigCellBlock[bIdx+2]=v6; bigCellBlock[bIdx+3]=v4; bIdx+=4;
                        bigCellBlock[bIdx+0]=v4; bigCellBlock[bIdx+1]=v7; bigCellBlock[bIdx+2]=v6; bigCellBlock[bIdx+3]=v3; bIdx+=4;
                        bigCellBlock[bIdx+0]=v1; bigCellBlock[bIdx+1]=v3; bigCellBlock[bIdx+2]=v4; bigCellBlock[bIdx+3]=v6; // Central
                    } else { // Alterning: important for mesh comformity
                        bigCellBlock[bIdx+0]=v0; bigCellBlock[bIdx+1]=v1; bigCellBlock[bIdx+2]=v2; bigCellBlock[bIdx+3]=v5; bIdx+=4;
                        bigCellBlock[bIdx+0]=v0; bigCellBlock[bIdx+1]=v2; bigCellBlock[bIdx+2]=v3; bigCellBlock[bIdx+3]=v7; bIdx+=4;
                        bigCellBlock[bIdx+0]=v4; bigCellBlock[bIdx+1]=v7; bigCellBlock[bIdx+2]=v5; bigCellBlock[bIdx+3]=v0; bIdx+=4;
                        bigCellBlock[bIdx+0]=v7; bigCellBlock[bIdx+1]=v6; bigCellBlock[bIdx+2]=v5; bigCellBlock[bIdx+3]=v2; bIdx+=4;
                        bigCellBlock[bIdx+0]=v0; bigCellBlock[bIdx+1]=v5; bigCellBlock[bIdx+2]=v2; bigCellBlock[bIdx+3]=v7; // Central
                    }

                }
            }
        }
    }

    // Déduction du type
    CellType cType = (elemType == BlockElementType::Hexahedron) ? CellType::HEX8 :
                     (elemType == BlockElementType::Edge) ? CellType::PRISM6 : CellType::TET4;

    // Insertion séquentielle ultra-rapide avec un buffer réutilisé (zéro allocation dynamique ici !)
    std::vector<CfemInt> tempBuffer(nodesPerCell);
    for (CfemInt c = 0; c < totalCells; ++c) {
        for (int n = 0; n < nodesPerCell; ++n) {
            tempBuffer[n] = bigCellBlock[c * nodesPerCell + n];
        }
        mesh.addCell(cType, tempBuffer);
    }

    mesh.buildMesh();
    if (order == MeshOrder::Quadratic)
        mesh.upgradeToQuadratic();

    // -------------------------------
    // 3. Boundary tagging
    // -------------------------------
    applyBoundaryTags_omp(mesh, nx, ny, nz, gridIds, names);
}


// -------------------------------
// Thread-safe boundary tagging (Optimisé)
// -------------------------------
void BlockMeshBuilder::applyBoundaryTags_omp(
    Mesh &mesh,
    CfemInt nx, CfemInt ny, CfemInt nz,
    const std::vector<CfemInt> &grid,
    const BlockBoundaryNames &names)
{
    CfemInt nx1 = nx + 1, ny1 = ny + 1;
    std::vector<int> invGrid(mesh.getNumVertices(), -1);
    for (size_t k = 0; k < grid.size(); ++k)
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

    const auto& sysFacets = sysBoundary->getFacetIds();
    CfemInt numSysFacets = sysFacets.size();

    #pragma omp parallel
    {
        // Tableaux locaux au thread (Évite les contentions de verrou global)
        std::vector<CfemInt> t_xminFacets, t_xmaxFacets, t_yminFacets, t_ymaxFacets, t_zminFacets, t_zmaxFacets;
        std::vector<CfemInt> t_xminVerts, t_xmaxVerts, t_yminVerts, t_ymaxVerts, t_zminVerts, t_zmaxVerts;

        #pragma omp for schedule(static)
        for (CfemInt f_idx = 0; f_idx < numSysFacets; ++f_idx)
        {
            CfemInt idx = sysFacets[f_idx];
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

            // Ajout dans les tableaux LOCAUX (Thread-safe par nature, pas de lock !)
            auto assignTo = [&](std::vector<CfemInt> &facetList, std::vector<CfemInt> &vertsList) {
                facetList.push_back(idx);
                vertsList.insert(vertsList.end(), fv.begin(), fv.end());
            };

            if (onXmin) assignTo(t_xminFacets, t_xminVerts);
            else if (onXmax) assignTo(t_xmaxFacets, t_xmaxVerts);
            else if (onYmin) assignTo(t_yminFacets, t_yminVerts);
            else if (onYmax) assignTo(t_ymaxFacets, t_ymaxVerts);
            else if (onZmin) assignTo(t_zminFacets, t_zminVerts); // Fixed typo from previous
            else if (onZmax) assignTo(t_zmaxFacets, t_zmaxVerts);
        }

        // Le verrou critique n'est appelé qu'UNE SEULE FOIS par thread (énorme gain de temps)
        #pragma omp critical
        {
            auto pushToEntity = [&](CfemInt entityId, const std::vector<CfemInt>& facets, const std::vector<CfemInt>& verts) {
                auto *entity = mesh.getEntity(entityId);
                for(auto f : facets) entity->addFacet(f);
                for(auto v : verts)  entity->addVertex(v);
            };

            pushToEntity(xminId, t_xminFacets, t_xminVerts);
            pushToEntity(xmaxId, t_xmaxFacets, t_xmaxVerts);
            pushToEntity(yminId, t_yminFacets, t_yminVerts);
            pushToEntity(ymaxId, t_ymaxFacets, t_ymaxVerts);
            pushToEntity(zminId, t_zminFacets, t_zminVerts);
            pushToEntity(zmaxId, t_zmaxFacets, t_zmaxVerts);
        }
    }

    // Finalize all entities
    for (auto eid : {xminId, xmaxId, yminId, ymaxId, zminId, zmaxId}) {
        mesh.getEntity(eid)->finalize();
    }
}
} // namespace cfem