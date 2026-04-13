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

#include "SphereMeshBuilder.h"

#include "Mesh.h"

namespace cfem {

    std::shared_ptr<Mesh> SphereMeshBuilder::create(
                                                MeshOrder order,
                                                CfemReal target_mesh_size,
                                                const cfem::Vector3D &center,
                                                CfemReal radius,
                                                SphereElementType type)
    {
        auto mesh = Mesh::create();

        SphereMeshBuilder builder(ConstructorToken{});

        builder.build(*mesh, order, target_mesh_size, center, radius, type);

        return mesh;
    }

    void SphereMeshBuilder::build(
                                Mesh &mesh,
                                MeshOrder order,
                                CfemReal h,
                                const cfem::Vector3D &center,
                                CfemReal radius,
                                SphereElementType type)
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
            CfemReal zeta = -1.0 + 2.0 * k / nz;

            for (int j = 0; j <= ny; ++j)
            {
                CfemReal eta = -1.0 + 2.0 * j / ny;

                for (int i = 0; i <= nx; ++i)
                {
                    CfemReal xi = -1.0 + 2.0 * i / nx;
                    CfemReal x = xi * std::sqrt(
                                        1.0 - 0.5 * eta * eta - 0.5 * zeta * zeta + eta * eta * zeta * zeta / 3.0);
                    CfemReal y = eta * std::sqrt(
                                        1.0 - 0.5 * xi * xi - 0.5 * zeta * zeta + xi * xi * zeta * zeta / 3.0);
                    CfemReal z = zeta * std::sqrt(
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

                    if (type == SphereElementType::Hexahedron)
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
            Mesh::projectBoundaryMidNodes(mesh, center, radius);
        }
    }


}