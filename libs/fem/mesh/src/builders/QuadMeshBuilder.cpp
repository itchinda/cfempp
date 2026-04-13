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

#include "QuadMeshBuilder.h"

namespace cfem
{
    /**
     * @brief Create a structured quadrilateral mesh with boundary naming.
     *
     * @param order Polynomial order of the mesh.
     * @param subDiv Number of subdivisions in x (nx) and y (ny).
     * @param corners Geometric corners of the quadrilateral domain.
     * @param elemType Type of elements to generate.
     * @param names Optional boundary names.
     *
     * @return Shared pointer to the generated Mesh.
     *
     * The domain is mapped using bilinear interpolation between the four corners.
     *
     * Boundary tags are assigned consistently with logical grid directions.
     */
    std::shared_ptr<Mesh> QuadMeshBuilder::create(MeshOrder                order,
                                                  QuadSubdivisions         subDiv,
                                                  const QuadCorners       &corners,
                                                  QuadElementType          elemType,
                                                  const QuadBoundaryNames &names)
    {
        auto mesh = Mesh::create();

        QuadMeshBuilder builder(ConstructorToken{});

        builder.build(*mesh, order, subDiv, corners, elemType, names);

        return mesh;
    }

    /**
     * @brief Internal mesh construction routine.
     *
     * Generates:
     *  - Structured grid nodes
     *  - Element connectivity
     *  - Boundary tags
     *
     * @param mesh Mesh object to populate.
     * @param order Polynomial order.
     * @param subDiv Structured subdivisions.
     * @param corners Domain geometry.
     * @param elemType Element topology.
     * @param names Boundary names.
     * 
     * Geometry mapping:
     *
     * The physical coordinates \f$ \displaystyle (x,y)\f$ are obtained from the reference
     * square \f$ \displaystyle (\zeta,eta) \in [0,1]^2\f$ using bilinear interpolation:
     *
     *     \f$ \displaystyle X(\zeta,\eta) = \sum N_i(\zeta,\eta) P_i \f$
     *
     * where \f$ \displaystyle P_i \f$ are the four corner points and Ni are bilinear
     * shape functions.
     *
     * This guarantees:
     *  - Affine mapping for parallelograms
     *  - Bilinear mapping for general quads
     */
    void QuadMeshBuilder::build(Mesh                     &mesh,
                                MeshOrder                order,
                                const QuadSubdivisions   &subDiv,
                                const QuadCorners        &corners,
                                QuadElementType           elemType,
                                const QuadBoundaryNames  &names)
    {
        PROFILE_SCOPE("QuadMeshBuilder::build");
        mesh.clear();
        mesh.setDimension(2);
        mesh.setOrder(MeshOrder::Linear);

        const int nodeRows = subDiv.nx + 1;
        const int nodeCols = subDiv.ny + 1;
        const int numQuads = subDiv.nx * subDiv.ny;

        int totalVertices = nodeRows * nodeCols;
        if (elemType == QuadElementType::TriangleCrissCross)
        {
            totalVertices += numQuads; // Add one vertex for the center of every quad
        }

        // Pre-allocate everything to avoid reallocations and memory spikes
        std::vector<CfemInt> gridIds(nodeRows * nodeCols);
        mesh.m_vertices.reserve(totalVertices);

        int cellMultiplier = (elemType == QuadElementType::TriangleCrissCross) ? 4 : 2;
        if (elemType == QuadElementType::Quadrilateral)
            cellMultiplier = 1;
        mesh.m_cells.reserve(numQuads * cellMultiplier);

        // Optimized Vertex Generation
        // We pre-calculate row vectors to turn bilinear interpolation into simple additions
        for (int j = 0; j < nodeCols; ++j)
        {
            double v = static_cast<double>(j) / subDiv.ny;
            // Points along the side edges at height v
            cfem::Vector3D leftP = corners.vertices[0] * (1.0 - v) + corners.vertices[3] * v;
            cfem::Vector3D rightP = corners.vertices[1] * (1.0 - v) + corners.vertices[2] * v;
            cfem::Vector3D rowVec = rightP - leftP;

            for (int i = 0; i < nodeRows; ++i)
            {
                double u = static_cast<double>(i) / subDiv.nx;
                cfem::Vector3D p = leftP + rowVec * u;
                gridIds[j * nodeRows + i] = mesh.addVertex(p.x, p.y, p.z);
            }
        }

        // Cell Generation: Move switch OUTSIDE the loop
        // This allows the CPU to stay in a tight, predictable loop for 8 million iterations
        auto generateCells = [&](auto func)
        {
            for (int j = 0; j < subDiv.ny; ++j)
            {
                for (int i = 0; i < subDiv.nx; ++i)
                {
                    CfemInt v0 = gridIds[j * nodeRows + i];
                    CfemInt v1 = gridIds[j * nodeRows + (i + 1)];
                    CfemInt v2 = gridIds[(j + 1) * nodeRows + (i + 1)];
                    CfemInt v3 = gridIds[(j + 1) * nodeRows + i];
                    func(v0, v1, v2, v3);
                }
            }
        };

        switch (elemType)
        {
        case QuadElementType::Quadrilateral:
            generateCells([&](CfemInt v0, CfemInt v1, CfemInt v2, CfemInt v3)
                          { mesh.addCell(CellType::QUAD4, {v0, v1, v2, v3}); });
            break;

        case QuadElementType::TriangleFirstDiag:
            generateCells([&](CfemInt v0, CfemInt v1, CfemInt v2, CfemInt v3)
                          {
            mesh.addCell(CellType::TRI3, {v0, v1, v2});
            mesh.addCell(CellType::TRI3, {v0, v2, v3}); });
            break;

        case QuadElementType::TriangleSecondDiag:
            generateCells([&](CfemInt v0, CfemInt v1, CfemInt v2, CfemInt v3)
                          {
            mesh.addCell(CellType::TRI3, {v1, v2, v3});
            mesh.addCell(CellType::TRI3, {v1, v3, v0}); });
            break;

        case QuadElementType::TriangleCrissCross:
            generateCells([&](CfemInt v0, CfemInt v1, CfemInt v2, CfemInt v3)
                          {
            // Midpoint calculation is still expensive, but necessary for CrissCross
            cfem::Vector3D p0 = mesh.getVertex(v0).coords;
            cfem::Vector3D p2 = mesh.getVertex(v2).coords;
            CfemInt vMid = mesh.addVertex((p0.x + p2.x)*0.5, (p0.y + p2.y)*0.5, (p0.z + p2.z)*0.5);
            mesh.addCell(CellType::TRI3, {v0, v1, vMid});
            mesh.addCell(CellType::TRI3, {v1, v2, vMid});
            mesh.addCell(CellType::TRI3, {v2, v3, vMid});
            mesh.addCell(CellType::TRI3, {v3, v0, vMid}); });
            break;
        }

        // 4. Finalize
        mesh.buildMesh();

        if (order == MeshOrder::Quadratic)
        {
            // Now that buildFacets is fast, this will also be significantly faster
            mesh.upgradeToQuadratic();
        }

        applyBoundaryTags(mesh, subDiv.nx, subDiv.ny, gridIds, names);
    }

    /**
     * @brief Assign boundary tags to exterior facets.
     *
     * Boundaries are identified from the structured grid topology:
     *  - xmin  -> left boundary
     *  - xmax  -> right boundary
     *  - ymin  -> bottom boundary
     *  - ymax  -> top boundary
     *
     * @param mesh Mesh to modify.
     * @param nx Number of subdivisions in x.
     * @param ny Number of subdivisions in y.
     * @param grid Structured vertex indexing.
     * @param names Boundary names.
     */
    void QuadMeshBuilder::applyBoundaryTags(Mesh &mesh, 
                                            CfemInt nx, 
                                            CfemInt ny,
                                            const std::vector<CfemInt> &grid,
                                            const QuadBoundaryNames &names)
    {
        // Create/Retrieve Mesh Entities
        CfemInt yminId = mesh.createEntity(names.ymin, 1);
        CfemInt xmaxId = mesh.createEntity(names.xmax, 1);
        CfemInt ymaxId = mesh.createEntity(names.ymax, 1);
        CfemInt xminId = mesh.createEntity(names.xmin, 1);

        int nodeRows = nx + 1;

        // OPTIMIZATION: Create an inverse map once to avoid std::find
        // vertex_id -> linear_index_in_grid
        // If 'grid' is just a permutation of 0..N, this is very fast.
        std::vector<int> invGrid(mesh.getNumVertices(), -1);
        for (CfemInt k = 0; k < grid.size(); ++k)
        {
            invGrid[grid[k]] = k;
        }

        // Iterate only over the System Boundary facets
        auto *sysBoundary = mesh.getEntity(cfem::SYS_BOUNDARY_ID);
        if (!sysBoundary)
            return;

        for (auto f_idx : sysBoundary->getFacetIds())
        {
            const Facet &facet = mesh.getFacet(f_idx);

            auto fv = mesh.getFacetConnectivity(f_idx);

            bool isBottom = true, isTop = true, isLeft = true, isRight = true;

            for (CfemInt vid : fv)
            {
                int idx = invGrid[vid];
                if (idx == -1) continue; // Should not happen in a valid grid

                CfemInt i = idx % nodeRows;
                CfemInt j = idx / nodeRows;

                if (j != 0) isBottom = false;
                if (j != ny) isTop = false;
                if (i != 0) isLeft = false;
                if (i != nx) isRight = false;
            }

            // Assignment Logic
            auto assignTo = [&](CfemInt entityId)
            {
                auto *entity = mesh.getEntity(entityId);
                entity->addFacet(f_idx);
                for (auto v_id : fv)
                    entity->addVertex(v_id);
            };

            if (isBottom) assignTo(yminId);
            else if (isTop) assignTo(ymaxId);
            else if (isLeft) assignTo(xminId);
            else if (isRight) assignTo(xmaxId);
        }

        // Finalize (removes duplicate vertices added during the loop)
        mesh.getEntity(yminId)->finalize();
        mesh.getEntity(ymaxId)->finalize();
        mesh.getEntity(xminId)->finalize();
        mesh.getEntity(xmaxId)->finalize();
    }

}