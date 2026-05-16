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


#include "SurfaceSphereBuilder.h"
#include <map>

namespace cfem {

    std::shared_ptr<Mesh> SurfaceSphereBuilder::create(
                                                MeshOrder order,
                                                CfemReal target_mesh_size,
                                                const cfem::Vector3D &center,
                                                CfemReal radius)
    {
        auto mesh = Mesh::create();

        SurfaceSphereBuilder builder(ConstructorToken{});

        builder.build(*mesh, order, target_mesh_size, center, radius);

        return mesh;
    }

    void SurfaceSphereBuilder::build(
        Mesh &mesh,
        MeshOrder order,
        CfemReal h,
        const cfem::Vector3D &center,
        CfemReal radius) 
    {
        mesh.clear();
        mesh.setOrder(MeshOrder::Linear);
        mesh.setTopologicalDimension(2); // Surface 2D dans l'espace 3D

        // Ajout des 12 sommets de base (normalisés et projetés)
        auto addProjectedVertex = [&](CfemReal x, CfemReal y, CfemReal z) {
            cfem::Vector3D p = (cfem::Vector3D(x, y, z)).normalized() * radius + center;
            return mesh.addVertex(p.x, p.y, p.z);
        };

        std::vector<CfemInt> v;
        const CfemReal t = (1.0 + std::sqrt(5.0)) / 2.0;
        v.push_back(addProjectedVertex(-1,  t,  0)); v.push_back(addProjectedVertex( 1,  t,  0));
        v.push_back(addProjectedVertex(-1, -t,  0)); v.push_back(addProjectedVertex( 1, -t,  0));
        v.push_back(addProjectedVertex( 0, -1,  t)); v.push_back(addProjectedVertex( 0,  1,  t));
        v.push_back(addProjectedVertex( 0, -1, -t)); v.push_back(addProjectedVertex( 0,  1, -t));
        v.push_back(addProjectedVertex( t,  0, -1)); v.push_back(addProjectedVertex( t,  0,  1));
        v.push_back(addProjectedVertex(-t,  0, -1)); v.push_back(addProjectedVertex(-t,  0,  1));

        // Faces initiales
        std::vector<std::vector<CfemInt>> faces = {
            {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
            {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
            {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
            {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}
        };

        // Calcul du nombre de subdivisions
        // L'arête initiale de l'icosaèdre est d'environ 1.05 * radius. 
        // Chaque subdivision divise l'arête par 2.
        CfemInt subdivisions = 0;
        CfemReal current_h = 1.05 * radius;
        while (current_h > h) { 
            current_h /= 2.0;
            subdivisions++;
        }

        std::unordered_map<uint64_t, CfemInt> midpoint_cache;

        auto getMidpoint = [&](CfemInt v1, CfemInt v2) {
            //uint64_t key = (v1 < v2) ? ((uint64_t)v1 << 32 | v2) : ((uint64_t)v2 << 32 | v1);
            uint64_t key = (uint64_t)std::min(v1, v2) << 32 | (uint64_t)std::max(v1, v2);
            auto it = midpoint_cache.find(key);
            if (it != midpoint_cache.end()) return it->second;

            const auto& p1 = mesh.getVertex(v1).coords;
            const auto& p2 = mesh.getVertex(v2).coords;
            cfem::Vector3D mid = center + ((p1 + p2) * 0.5 - center).normalized() * radius;
            
            CfemInt new_id = mesh.addVertex(mid.x, mid.y, mid.z);
            midpoint_cache[key] = new_id;
            return new_id;
        };

        for (CfemInt s = 0; s < subdivisions; ++s) {
            std::vector<std::vector<CfemInt>> next_faces;
            next_faces.reserve(faces.size() * 4); 
            midpoint_cache.reserve(faces.size() * 3);

            for (const auto& f : faces) {
                CfemInt a = getMidpoint(f[0], f[1]);
                CfemInt b = getMidpoint(f[1], f[2]);
                CfemInt c = getMidpoint(f[2], f[0]);

                next_faces.push_back({f[0], a, c});
                next_faces.push_back({f[1], b, a});
                next_faces.push_back({f[2], c, b});
                next_faces.push_back({a, b, c});
            }
            faces = std::move(next_faces);
            midpoint_cache.clear(); // On vide le cache car les arêtes changent d'ID à chaque niveau
        }

        for (const auto& f : faces) {
            mesh.addCell(CellType::TRI3, f);
        }

        mesh.buildMesh();

        if (order == MeshOrder::Quadratic) {
            mesh.upgradeToQuadratic();
            for (CfemInt i = 0; i < mesh.getNumVertices(); ++i) {
                auto& v = mesh.getVertex(i).coords;
                v = center + (v - center).normalized() * radius;
            }
        }
        
        mesh.updateAffineFlags();
    }
} // namespace cfem