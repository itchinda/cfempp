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

#include "CylinderMeshBuilder.h"

namespace cfem
{
    std::shared_ptr<Mesh> CylinderMeshBuilder::create(
        MeshOrder order, CfemReal target_mesh_size,
        const cfem::Vector3D &center, CfemReal radius,
        CfemReal height, const cfem::Vector3D &axis,
        CylinderElementType type, const CylinderBoundaryNames &names)
    {
        auto mesh = Mesh::create();
        CylinderMeshBuilder builder(ConstructorToken{});
        builder.build(*mesh, order, target_mesh_size, center, radius, height, axis, type, names);
        return mesh;
    }

    inline void CylinderMeshBuilder::build(
        Mesh &mesh, MeshOrder order, CfemReal h,
        const cfem::Vector3D &center, CfemReal radius,
        CfemReal height, const cfem::Vector3D &axis,
        CylinderElementType type, const CylinderBoundaryNames &names)
    {
        mesh.clear();
        mesh.setOrder(MeshOrder::Linear);
        mesh.setDimension(3);

        // ------------------------
        // 1. Layer Structure
        // ------------------------
        int num_rings = std::max(1, static_cast<int>(std::round(radius / h)));
        CfemReal dr = radius / num_rings;
        int nz = std::max(1, static_cast<int>(std::round(height / h)));
        CfemReal dz = height / nz;
        CfemReal z_min = -height * 0.5; // generate around origin

        std::vector<std::vector<std::vector<CfemInt>>> grid_nodes(nz + 1);

        // ------------------------
        // 2. Vertex Generation
        // ------------------------
        for (int k = 0; k <= nz; ++k) {
            grid_nodes[k].resize(num_rings + 1);
            CfemReal z = z_min + k * dz;

            // Center node
            grid_nodes[k][0].push_back(mesh.addVertex(0.0, 0.0, z));

            for (int r = 1; r <= num_rings; ++r) {
                CfemReal current_r = r * dr;
                int num_points = std::max(3, static_cast<int>(std::round(2.0 * cfem::pi * current_r / h)));

                for (int i = 0; i < num_points; ++i) {
                    CfemReal theta = (2.0 * cfem::pi * i) / num_points;
                    CfemReal x = current_r * std::cos(theta);
                    CfemReal y = current_r * std::sin(theta);
                    grid_nodes[k][r].push_back(mesh.addVertex(x, y, z));
                }
            }
        }

        // ------------------------
        // 3. Helper: Build Elements
        // ------------------------
        auto addExtrudedVolume = [&](CfemInt b0, CfemInt b1, CfemInt b2,
                                     CfemInt t0, CfemInt t1, CfemInt t2) {
            const auto& p0 = mesh.getVertex(b0).coords;
            const auto& p1 = mesh.getVertex(b1).coords;
            const auto& p2 = mesh.getVertex(b2).coords;

            CfemReal signedArea = (p1.x - p0.x)*(p2.y - p0.y) - (p1.y - p0.y)*(p2.x - p0.x);

            if (signedArea < 0) { std::swap(b1,b2); std::swap(t1,t2); }
            else if (std::abs(signedArea) < 1e-15) return;

            if (type == CylinderElementType::Wedge)
                mesh.addCell(CellType::PRISM6, {b0,b1,b2,t0,t1,t2});
            else {
                mesh.addCell(CellType::TET4, {b0,b1,b2,t0});
                mesh.addCell(CellType::TET4, {t0,b1,b2,t1});
                mesh.addCell(CellType::TET4, {t0,t1,b2,t2});
            }
        };

        // ------------------------
        // 4. Connectivity
        // ------------------------
        for (int k = 0; k < nz; ++k) {
            const auto& bot_center = grid_nodes[k][0][0];
            const auto& top_center = grid_nodes[k+1][0][0];
            const auto& bot_r1 = grid_nodes[k][1];
            const auto& top_r1 = grid_nodes[k+1][1];

            // Inner triangles
            for (size_t i = 0; i < bot_r1.size(); ++i) {
                size_t next_i = (i+1)%bot_r1.size();
                addExtrudedVolume(bot_center, bot_r1[i], bot_r1[next_i],
                                  top_center, top_r1[i], top_r1[next_i]);
            }

            // Outer rings
            for (int r = 1; r < num_rings; ++r) {
                const auto& in_b = grid_nodes[k][r];   const auto& out_b = grid_nodes[k][r+1];
                const auto& in_t = grid_nodes[k+1][r]; const auto& out_t = grid_nodes[k+1][r+1];
                size_t n_in = in_b.size(), n_out = out_b.size(), j = 0;

                for (size_t i = 0; i < n_in; ++i) {
                    size_t next_i = (i+1)%n_in;
                    size_t target_j = static_cast<size_t>((static_cast<double>(next_i)/n_in)*n_out);
                    while (j != target_j || (next_i==0 && j!=0)) {
                        size_t next_j = (j+1)%n_out;
                        addExtrudedVolume(in_b[i], out_b[next_j], out_b[j],
                                          in_t[i], out_t[next_j], out_t[j]);
                        j = next_j;
                        if (j==0 && next_i!=0) break;
                    }
                    addExtrudedVolume(in_b[i], in_b[next_i], out_b[j],
                                      in_t[i], in_t[next_i], out_t[j]);
                }
            }
        }

        mesh.buildMesh();

        // ------------------------
        // 5. Boundary tagging
        // ------------------------
        CfemInt latId = mesh.createEntity(names.lateral, 2);
        CfemInt botId = mesh.createEntity(names.bottom, 2);
        CfemInt topId = mesh.createEntity(names.top, 2);

        auto* sysBoundary = mesh.getEntity(cfem::SYS_BOUNDARY_ID);
        if (sysBoundary) {
            for (auto f_idx : sysBoundary->getFacetIds()) {
                auto fv = mesh.getFacetConnectivity(f_idx);
                bool on_bot=true, on_top=true, on_lat=true;
                for (CfemInt vid : fv) {
                    const auto& p = mesh.getVertex(vid).coords;
                    if (std::abs(p.z - z_min) > 1e-6) on_bot=false;
                    if (std::abs(p.z - (z_min+height)) > 1e-6) on_top=false;
                    CfemReal dr_current = std::sqrt(p.x*p.x + p.y*p.y);
                    if (std::abs(dr_current - radius) > 1e-6) on_lat=false;
                }
                auto assignTo = [&](CfemInt entityId){ auto* e=mesh.getEntity(entityId); e->addFacet(f_idx); for(auto vid : fv) e->addVertex(vid); };
                if (on_bot) assignTo(botId);
                else if (on_top) assignTo(topId);
                else if (on_lat) assignTo(latId);
            }
        }
        mesh.getEntity(latId)->finalize();
        mesh.getEntity(botId)->finalize();
        mesh.getEntity(topId)->finalize();

        // ------------------------
        // 6. Quadratic Upgrade
        // ------------------------
        if (order == MeshOrder::Quadratic) {
            mesh.upgradeToQuadratic();
            auto* latEntity = mesh.getEntity(latId);
            if (latEntity) {
                for (CfemInt vid : latEntity->getVertexIds()) {
                    auto& p = mesh.getVertex(vid).coords;
                    CfemReal current_r = std::sqrt(p.x*p.x + p.y*p.y);
                    if (current_r > 1e-12) { p.x = (p.x/current_r)*radius; p.y = (p.y/current_r)*radius; }
                }
            }
            mesh.updateAffineFlags();
        }

        // ------------------------
        // 7. Axis Alignment + Translation
        // ------------------------
        cfem::Vector3D target_axis = axis;
        CfemReal axis_len = std::sqrt(target_axis.x*target_axis.x + target_axis.y*target_axis.y + target_axis.z*target_axis.z);
        if (axis_len < ZERO_TOLERANCE) CFEM_ERROR << "Cylinder axis has zero length.";
        target_axis.x/=axis_len; target_axis.y/=axis_len; target_axis.z/=axis_len;

        CfemReal cos_theta = target_axis.z;
        CfemReal sin_theta = std::sqrt(1.0 - cos_theta*cos_theta);

        cfem::Vector3D k = {0,0,1};
        if (sin_theta > ZERO_TOLERANCE) {
            k.x = -target_axis.y; k.y = target_axis.x; k.z=0;
            CfemReal k_len = std::sqrt(k.x*k.x + k.y*k.y); k.x/=k_len; k.y/=k_len;
        }

        for (CfemInt vid=0; vid<mesh.getNumVertices(); ++vid) {
            auto &p = mesh.getVertex(vid).coords;
            CfemReal x = p.x, y = p.y, z = p.z;

            if (sin_theta > ZERO_TOLERANCE) {
                CfemReal v_dot_k = x*k.x + y*k.y + z*k.z;
                CfemReal cross_x = k.y*z - k.z*y;
                CfemReal cross_y = k.z*x - k.x*z;
                CfemReal cross_z = k.x*y - k.y*x;
                CfemReal xr = x*cos_theta + cross_x*sin_theta + k.x*v_dot_k*(1-cos_theta);
                CfemReal yr = y*cos_theta + cross_y*sin_theta + k.y*v_dot_k*(1-cos_theta);
                CfemReal zr = z*cos_theta + cross_z*sin_theta + k.z*v_dot_k*(1-cos_theta);
                x = xr; y = yr; z = zr;
            }

            // Translate to desired center
            p.x = x + center.x;
            p.y = y + center.y;
            p.z = z + center.z;
        }
    }
}