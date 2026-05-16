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

#include "DiskMeshBuilder.h"

namespace cfem
{
    std::shared_ptr<Mesh> DiskMeshBuilder::create(
        MeshOrder order,
        CfemReal target_mesh_size,
        const cfem::Vector3D &center,
        CfemReal radius,
        const cfem::Vector3D &normal,
        DiskElementType type)
    {
        auto mesh = Mesh::create();
        DiskMeshBuilder builder(ConstructorToken{});
        builder.build(*mesh, order, target_mesh_size, center, radius, normal, type);
        return mesh;
    }

    inline void DiskMeshBuilder::build(
        Mesh &mesh,
        MeshOrder order,
        CfemReal h,
        const cfem::Vector3D &center,
        CfemReal radius,
        const cfem::Vector3D &normal,
        DiskElementType type)
    {
        mesh.clear();
        mesh.setOrder(MeshOrder::Linear);
        mesh.setTopologicalDimension(2);

        // 1. Layer Structure
        CfemInt num_rings = static_cast<int>(std::round(radius / h));
        if (num_rings < 1) num_rings = 1;
        CfemReal dr = radius / num_rings;

        std::vector<std::vector<CfemInt>> ring_nodes(num_rings + 1);

        // 2. Compute orthonormal basis (u,v) in the plane of the disk
        cfem::Vector3D u, v;
        computeOrthonormalBasis(normal, u, v);

        // 3. Vertex generation
        ring_nodes[0].push_back(mesh.addVertex(center.x, center.y, center.z)); // center

        for (CfemInt r = 1; r <= num_rings; ++r) {
            CfemReal current_r = r * dr;
            CfemInt num_points = static_cast<int>(std::round(2.0 * M_PI * current_r / h));
            if (num_points < 3) num_points = 3;

            for (CfemInt i = 0; i < num_points; ++i) {
                CfemReal theta = (2.0 * M_PI * i) / num_points;
                cfem::Vector3D pos = center + u * (current_r * std::cos(theta)) + v * (current_r * std::sin(theta));
                ring_nodes[r].push_back(mesh.addVertex(pos.x, pos.y, pos.z));
            }
        }

        // 4. Helper to guarantee CCW orientation
        auto addCCWTriangle = [&](CfemInt v0, CfemInt v1, CfemInt v2) {
            const auto& p0 = mesh.getVertex(v0).coords;
            const auto& p1 = mesh.getVertex(v1).coords;
            const auto& p2 = mesh.getVertex(v2).coords;

            // Project onto plane defined by normal
            cfem::Vector3D n = normal.normalized();
            cfem::Vector3D proj1 = p1 - p0;
            cfem::Vector3D proj2 = p2 - p0;
            cfem::Vector3D cross = proj1.cross(proj2);
            CfemReal signedArea = cross.dot(n);

            if (signedArea < 0)
                mesh.addCell(CellType::TRI3, {v0, v2, v1});
            else if (std::abs(signedArea) > ZERO_TOLERANCE)
                mesh.addCell(CellType::TRI3, {v0, v1, v2});
        };

        // 5. Triangles for first ring (center)
        for (size_t i = 0; i < ring_nodes[1].size(); ++i) {
            size_t next_i = (i + 1) % ring_nodes[1].size();
            addCCWTriangle(ring_nodes[0][0], ring_nodes[1][i], ring_nodes[1][next_i]);
        }

        // 6. Tessellation for outer rings
        for (CfemInt r = 1; r < num_rings; ++r) {
            const auto& inner = ring_nodes[r];
            const auto& outer = ring_nodes[r+1];

            size_t n_in = inner.size();
            size_t n_out = outer.size();
            size_t j = 0;

            for (size_t i = 0; i < n_in; ++i) {
                size_t next_i = (i + 1) % n_in;
                size_t target_j = static_cast<size_t>((static_cast<double>(next_i)/n_in) * n_out);

                while (j != target_j || (next_i == 0 && j != 0)) {
                    size_t next_j = (j + 1) % n_out;
                    addCCWTriangle(inner[i], outer[next_j], outer[j]);
                    j = next_j;
                    if (j == 0 && next_i != 0) break;
                }
                addCCWTriangle(inner[i], inner[next_i], outer[j]);
            }
        }

        mesh.buildMesh();

        // 7. Relaxation (optional)
        mesh.relaxMesh(20);

        // 8. Quadratic upgrade
        if (order == MeshOrder::Quadratic) {
            mesh.upgradeToQuadratic();
            Mesh::projectBoundaryMidNodes(mesh, center, radius);
        }
        mesh.updateAffineFlags();
    }
}