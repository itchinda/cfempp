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

#include "SolidSphereBuilder.h"

#include "Mesh.h"

namespace cfem {

    std::shared_ptr<Mesh> SolidSphereBuilder::create(
                                    MeshOrder order,
                                    CfemReal target_mesh_size,
                                    const cfem::Vector3D &center,
                                    CfemReal radius,
                                    SolidSphereElementType type)
    {
        auto mesh = Mesh::create();

        SolidSphereBuilder builder(ConstructorToken{});

        builder.build(*mesh, order, target_mesh_size, center, radius, type);

        return mesh;
    }

    void SolidSphereBuilder::build(
                                Mesh &mesh,
                                MeshOrder order,
                                CfemReal h,
                                const cfem::Vector3D &center,
                                CfemReal radius,
                                SolidSphereElementType type)
    {
        mesh.clear();
        mesh.setOrder(MeshOrder::Linear);
        mesh.setTopologicalDimension(3);

        CfemInt n = static_cast<CfemInt>(std::ceil((2.0 * radius) / h));
        if (n < 1)
            n = 1;

        CfemInt nx = n;
        CfemInt ny = n;
        CfemInt nz = n;

        std::vector<CfemInt> ids((nx + 1) * (ny + 1) * (nz + 1));

        auto idx = [&](CfemInt i, CfemInt j, CfemInt k)
        {
            return k * (nx + 1) * (ny + 1) + j * (nx + 1) + i;
        };

        /*----------------------------
        1. Cube → Ball mapping
        -----------------------------*/
        for (CfemInt k = 0; k <= nz; ++k)
        {
            CfemReal zeta = -1.0 + 2.0 * k / nz;

            for (CfemInt j = 0; j <= ny; ++j)
            {
                CfemReal eta = -1.0 + 2.0 * j / ny;

                for (CfemInt i = 0; i <= nx; ++i)
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
        for (CfemInt k = 0; k < nz; ++k)
        {
            for (CfemInt j = 0; j < ny; ++j)
            {
                for (CfemInt i = 0; i < nx; ++i)
                {

                    CfemInt v0 = ids[idx(i, j, k)];
                    CfemInt v1 = ids[idx(i + 1, j, k)];
                    CfemInt v2 = ids[idx(i + 1, j + 1, k)];
                    CfemInt v3 = ids[idx(i, j + 1, k)];

                    CfemInt v4 = ids[idx(i, j, k + 1)];
                    CfemInt v5 = ids[idx(i + 1, j, k + 1)];
                    CfemInt v6 = ids[idx(i + 1, j + 1, k + 1)];
                    CfemInt v7 = ids[idx(i, j + 1, k + 1)];

                    if (type == SolidSphereElementType::Hexahedron)
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
                        // CfemInt vMid = mesh.addVertex(midP.x, midP.y, midP.z);
                        // 6-tet decomposition
                        // mesh.addCell(CellType::TET4, {v0, v1, v3, v4});
                        // mesh.addCell(CellType::TET4, {v1, v3, v4, v5});
                        // mesh.addCell(CellType::TET4, {v4, v5, v3, v7});
                        // mesh.addCell(CellType::TET4, {v1, v3, v5, v2});
                        // mesh.addCell(CellType::TET4, {v3, v7, v5, v2}); 
                        // mesh.addCell(CellType::TET4, {v5, v7, v6, v2});

                        bool isEven = (i + j + k) % 2 == 0;
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

        mesh.relaxMesh(10);

        if (order == MeshOrder::Quadratic)
        {
            mesh.upgradeToQuadratic();
            Mesh::projectBoundaryMidNodes(mesh, center, radius);
            mesh.updateAffineFlags();
        }


    }

// void SolidSphereBuilder::build(
//     Mesh &mesh,
//     MeshOrder order,
//     CfemReal h,
//     const cfem::Vector3D &center,
//     CfemReal radius,
//     SolidSphereElementType type)
// {
//     mesh.clear();
//     mesh.setTopologicalDimension(3);
//     mesh.setOrder(MeshOrder::Linear);

//     const CfemInt N = std::max(2, (int)std::ceil(2.0 * radius / h));
//     const CfemReal dx = 2.0 * radius / N;

//     auto index = [&](CfemInt i,CfemInt j,CfemInt k){
//         return i + (N+1)*(j + (N+1)*k);
//     };

//     std::vector<CfemInt> vid((N+1)*(N+1)*(N+1), -1);

//     // --------------------------------------------------
//     // 1. Génération des sommets
//     // --------------------------------------------------
//     for (CfemInt k = 0; k <= N; ++k)
//     for (CfemInt j = 0; j <= N; ++j)
//     for (CfemInt i = 0; i <= N; ++i)
//     {
//         cfem::Vector3D p;
//         p.x = -radius + i * dx;
//         p.y = -radius + j * dx;
//         p.z = -radius + k * dx;

//         p += center;

//         CfemReal dist = (p - center).norm();

//         // ✔️ projection douce UNIQUEMENT près de la frontière
//         if (dist > 0.85 * radius)
//         {
//             if (dist > 1e-12)
//                 p = center + radius * (p - center) / dist;
//         }

//         vid[index(i,j,k)] = mesh.addVertex(p.x, p.y, p.z);
//     }

//     // --------------------------------------------------
//     // 2. Construction des cellules
//     // --------------------------------------------------
//     for (CfemInt k = 0; k < N; ++k)
//     for (CfemInt j = 0; j < N; ++j)
//     for (CfemInt i = 0; i < N; ++i)
//     {
//         CfemInt v[8] = {
//             vid[index(i  ,j  ,k  )],
//             vid[index(i+1,j  ,k  )],
//             vid[index(i+1,j+1,k  )],
//             vid[index(i  ,j+1,k  )],
//             vid[index(i  ,j  ,k+1)],
//             vid[index(i+1,j  ,k+1)],
//             vid[index(i+1,j+1,k+1)],
//             vid[index(i  ,j+1,k+1)]
//         };

//         // centre du cube
//         cfem::Vector3D c(0,0,0);
//         for(CfemInt m=0;m<8;m++)
//             c += mesh.getVertex(v[m]).coords;
//         c /= 8.0;

//         if ((c - center).norm() > radius)
//             continue;

//         if (type == SolidSphereElementType::Hexahedron)
//         {
//             mesh.addCell(CellType::HEX8, {
//                 v[0],v[1],v[2],v[3],
//                 v[4],v[5],v[6],v[7]
//             });
//         }
//         else
//         {
//             bool isEven = (i + j + k) % 2 == 0;
//             if (isEven) {
//                 // Ajouter les 5 tets de la Version A
//                 mesh.addCell(CellType::TET4, {v[0], v[1], v[3], v[4]});
//                 mesh.addCell(CellType::TET4, {v[1], v[2], v[3], v[6]});
//                 mesh.addCell(CellType::TET4, {v[4], v[5], v[6], v[1]});
//                 mesh.addCell(CellType::TET4, {v[4], v[6], v[7], v[3]});
//                 mesh.addCell(CellType::TET4, {v[1], v[3], v[4], v[6]}); // Central
//             } else {
//                 // Ajouter les 5 tets de la Version B
//                 mesh.addCell(CellType::TET4, {v[1], v[0], v[2], v[5]});
//                 mesh.addCell(CellType::TET4, {v[3], v[0], v[2], v[7]});
//                 mesh.addCell(CellType::TET4, {v[5], v[4], v[7], v[0]});
//                 mesh.addCell(CellType::TET4, {v[7], v[6], v[5], v[2]});
//                 mesh.addCell(CellType::TET4, {v[0], v[2], v[5], v[7]}); // Central
//             }
//         }
//     }

//     mesh.buildMesh();

//     // --------------------------------------------------
//     // 3. Relaxation (très important)
//     // --------------------------------------------------
//     mesh.relaxMesh(10);

//     // --------------------------------------------------
//     // 4. Quadratique
//     // --------------------------------------------------
//     if (order == MeshOrder::Quadratic)
//     {
//         mesh.upgradeToQuadratic();
//         Mesh::projectBoundaryMidNodes(mesh, center, radius);
//     }

//     mesh.updateAffineFlags();
// }
}