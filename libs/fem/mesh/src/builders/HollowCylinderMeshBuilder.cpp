#include "HollowCylinderMeshBuilder.h"
#include "numeric_policy.h"

namespace cfem
{
    std::shared_ptr<Mesh> HollowCylinderMeshBuilder::create(
        MeshOrder order, CfemReal target_mesh_size,
        const cfem::Vector3D &center, CfemReal r_in, CfemReal r_out, CfemReal height,
        const cfem::Vector3D &axis, CylinderElementType type,
        const HollowCylinderBoundaryNames &names)
    {
        auto mesh = Mesh::create();
        HollowCylinderMeshBuilder builder(ConstructorToken{});
        builder.build(*mesh, order, target_mesh_size, center, r_in, r_out, height, axis, type, names);
        return mesh;
    }

    inline void HollowCylinderMeshBuilder::build(
        Mesh &mesh, MeshOrder order, CfemReal h,
        const cfem::Vector3D &center, CfemReal r_in, CfemReal r_out, CfemReal height,
        const cfem::Vector3D &axis, CylinderElementType type,
        const HollowCylinderBoundaryNames &names)
    {
        if (r_in < ZERO_TOLERANCE)
            CFEM_ERROR << "HollowCylinderMeshBuilder requires r_in > 0.";
        mesh.clear();
        mesh.setOrder(MeshOrder::Linear);
        mesh.setDimension(3);

        // Normalize axis
        CfemReal axis_len = std::sqrt(axis.x * axis.x + axis.y * axis.y + axis.z * axis.z);
        if (axis_len < ZERO_TOLERANCE)
            CFEM_ERROR << "Cylinder axis has zero length.";
        CfemReal ux = axis.x / axis_len;
        CfemReal uy = axis.y / axis_len;
        CfemReal uz = axis.z / axis_len;

        // 1. Layer Structure Calculation
        CfemReal thickness = r_out - r_in;
        int num_rings = std::max(1, static_cast<int>(std::round(thickness / h)));
        CfemReal dr = thickness / num_rings;

        int nz = std::max(1, static_cast<int>(std::round(height / h)));
        CfemReal dz = height / nz;

        // BUILD AT ORIGIN: Center is assumed (0,0,0) for generation
        CfemReal z_min = -height * 0.5;

        std::vector<std::vector<std::vector<CfemInt>>> grid_nodes(nz + 1);

        // 2. Vertex Generation (At Origin)
        for (int k = 0; k <= nz; ++k)
        {
            grid_nodes[k].resize(num_rings + 1);
            CfemReal z = z_min + k * dz;

            for (int r = 0; r <= num_rings; ++r)
            {
                CfemReal current_r = r_in + r * dr;
                int num_points = std::max(6, static_cast<int>(std::round(2.0 * cfem::pi * current_r / h)));

                for (int i = 0; i < num_points; ++i)
                {
                    CfemReal theta = (2.0 * cfem::pi * i) / num_points;
                    CfemReal x = current_r * std::cos(theta);
                    CfemReal y = current_r * std::sin(theta);
                    grid_nodes[k][r].push_back(mesh.addVertex(x, y, z));
                }
            }
        }

        // 3. Extruded Volume Helper (Omitted for brevity, unchanged from previous)
        auto addExtrudedVolume = [&](CfemInt b0, CfemInt b1, CfemInt b2, CfemInt t0, CfemInt t1, CfemInt t2)
        {
            const auto &p0 = mesh.getVertex(b0).coords;
            const auto &p1 = mesh.getVertex(b1).coords;
            const auto &p2 = mesh.getVertex(b2).coords;
            
            CfemReal dx1 = p1.x - p0.x;
            CfemReal dy1 = p1.y - p0.y;
            CfemReal dx2 = p2.x - p0.x;
            CfemReal dy2 = p2.y - p0.y;

            CfemReal signedArea = dx1*dy2 - dy1*dx2;

            CfemReal norm_scale = std::hypot(dx1, dy1) * std::hypot(dx2, dy2);
            if (signedArea < 0)
            {
                std::swap(b1, b2);
                std::swap(t1, t2);
            }
            else if (NumericPolicy::isZero(signedArea, norm_scale))
                return;

            if (type == CylinderElementType::Wedge)
            {
                mesh.addCell(CellType::PRISM6, {b0, b1, b2, t0, t1, t2});
            }
            else
            {
                // Standard Wedge-to-Tet split (3 Tets)
                mesh.addCell(CellType::TET4, {b0, b1, b2, t0});
                mesh.addCell(CellType::TET4, {t0, b1, b2, t1});
                mesh.addCell(CellType::TET4, {t0, t1, b2, t2});
            }
        };

        // 4. Element Connectivity Generation (Unchanged)
        for (int k = 0; k < nz; ++k)
        {
            for (int r = 0; r < num_rings; ++r)
            {
                const auto &in_b = grid_nodes[k][r];
                const auto &out_b = grid_nodes[k][r + 1];
                const auto &in_t = grid_nodes[k + 1][r];
                const auto &out_t = grid_nodes[k + 1][r + 1];
                size_t n_in = in_b.size(), n_out = out_b.size(), j = 0;
                for (size_t i = 0; i < n_in; ++i)
                {
                    size_t next_i = (i + 1) % n_in;
                    size_t target_j = static_cast<size_t>((static_cast<double>(next_i) / n_in) * n_out);
                    while (j != target_j || (next_i == 0 && j != 0))
                    {
                        size_t next_j = (j + 1) % n_out;
                        addExtrudedVolume(in_b[i], out_b[next_j], out_b[j], in_t[i], out_t[next_j], out_t[j]);
                        j = next_j;
                        if (j == 0 && next_i != 0)
                            break;
                    }
                    addExtrudedVolume(in_b[i], in_b[next_i], out_b[j], in_t[i], in_t[next_i], out_t[j]);
                }
            }
        }

        mesh.buildMesh();
        //mesh.relaxMesh(20);

        // 5. Boundary Tagging (At Origin)
        CfemInt inId = mesh.createEntity(names.inner, 2);
        CfemInt outId = mesh.createEntity(names.outer, 2);
        CfemInt botId = mesh.createEntity(names.bottom, 2);
        CfemInt topId = mesh.createEntity(names.top, 2);

        auto *sysBoundary = mesh.getEntity(cfem::SYS_BOUNDARY_ID);
        if (sysBoundary)
        {
            CfemReal tol = 10 * ZERO_TOLERANCE * std::max({r_out, height, CfemReal(1.0)});
            for (auto f_idx : sysBoundary->getFacetIds())
            {
                auto fv = mesh.getFacetConnectivity(f_idx);
                bool on_bot = true, on_top = true, on_in = true, on_out = true;

                for (CfemInt vid : fv)
                {
                    const auto &p = mesh.getVertex(vid).coords;
                    if (!NumericPolicy::areEqual(p.z, z_min, height) )
                        on_bot = false;
                    if (!NumericPolicy::areEqual(p.z, height * 0.5, height ))
                        on_top = false;

                    CfemReal dr_current = std::hypot(p.x, p.y);
                    if (!NumericPolicy::areEqual(dr_current, r_in, r_in))
                        on_in = false;
                    if (!NumericPolicy::areEqual(dr_current, r_out, r_out))
                        on_out = false;
                }

                auto assignTo = [&](CfemInt entityId)
                {
                    auto *entity = mesh.getEntity(entityId);
                    entity->addFacet(f_idx);
                    for (CfemInt vid : fv)
                        entity->addVertex(vid);
                };

                if (on_in)
                    assignTo(inId);
                else if (on_out)
                    assignTo(outId);
                else if (on_bot)
                    assignTo(botId);
                else if (on_top)
                    assignTo(topId);
            }
        }
        mesh.getEntity(inId)->finalize();
        mesh.getEntity(outId)->finalize();
        mesh.getEntity(botId)->finalize();
        mesh.getEntity(topId)->finalize();

        // 6. Quadratic Upgrade & Projection (Still at Origin!)
        if (order == MeshOrder::Quadratic)
        {
            mesh.upgradeToQuadratic();
            projectQuadraticSurfaces(mesh, r_in, r_out, names);
        }

        // ====================================================================
        // 7. RIGID BODY TRANSFORMATION (Rotation + Translation)
        // ====================================================================

        // Loop over all vertices in the mesh.
        // (Assuming getNumVertices() gets the total count. Adapt to your API if needed)
        for (CfemInt vid = 0; vid < mesh.getNumVertices(); ++vid)
        {
            auto &p = mesh.getVertex(vid).coords;
            CfemReal x = p.x, y = p.y, z = p.z;

            // Step A: Rotation
            if (NumericPolicy::areEqual(uz, -1.0))
            {
                // Special Case: axis is (0,0,-1). Rotate 180 deg around X-axis.
                y = -y;
                z = -z;
            }
            else if (uz < 1.0 - ZERO_TOLERANCE)
            {
                // General Case: Rodrigues' rotation matrix
                CfemReal h = 1.0 / (1.0 + uz);
                CfemReal vx = -uy;
                CfemReal vy = ux;

                CfemReal rx = (uz + h * vx * vx) * x + (h * vx * vy) * y + (vy)*z;
                CfemReal ry = (h * vx * vy) * x + (uz + h * vy * vy) * y + (-vx) * z;
                CfemReal rz = (-vy) * x + (vx)*y + (uz)*z;

                x = rx;
                y = ry;
                z = rz;
            }

            // Step B: Translation
            p.x = x + center.x;
            p.y = y + center.y;
            p.z = z + center.z;
        }
    }

    void HollowCylinderMeshBuilder::projectQuadraticSurfaces(
                                                                Mesh &mesh, 
                                                                CfemReal r_in, 
                                                                CfemReal r_out,
                                                                const HollowCylinderBoundaryNames &names
                                                            )
    {
        auto project = [&](const std::string &name, CfemReal target_radius)
        {
            auto entity = mesh.getEntity(name);
            if (!entity)
                return;

            entity->update();

            for (CfemInt vid : entity->getVertexIds())
            {
                auto &node = mesh.getVertex(vid);
                // Center is (0,0,0) here because projection happens before translation
                CfemReal current_r = std::hypot(node.coords.x, node.coords.y);

                if (current_r > ZERO_TOLERANCE)
                {
                    node.coords.x = (node.coords.x / current_r) * target_radius;
                    node.coords.y = (node.coords.y / current_r) * target_radius;
                }
            }
        };

        project(names.inner, r_in);
        project(names.outer, r_out);
    }
}