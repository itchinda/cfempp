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
    mesh.setDimension(2); // Surface 2D dans l'espace 3D

    // 1. Calcul du nombre de subdivisions nécessaires pour atteindre h
    // La longueur d'arête d'un icosaèdre unité est d'environ 1.05.
    // Après n subdivisions, h_approx = (radius * 1.05) / 2^n
    int subdivisions = 0;
    CfemReal current_h = radius * 1.05;
    while (current_h > h && subdivisions < 8) {
        current_h /= 2.0;
        subdivisions++;
    }

    // 2. Génération des sommets de l'icosaèdre de base
    const CfemReal t = (1.0 + std::sqrt(5.0)) / 2.0;
    std::vector<cfem::Vector3D> base_coords = {
        {-1,  t,  0}, { 1,  t,  0}, {-1, -t,  0}, { 1, -t,  0},
        { 0, -1,  t}, { 0,  1,  t}, { 0, -1, -t}, { 0,  1, -t},
        { t,  0, -1}, { t,  0,  1}, {-t,  0, -1}, {-t,  0,  1}
    };

    for (auto &v : base_coords) {
        v = center + v.normalized() * radius;
        mesh.addVertex(v.x, v.y, v.z);
    }

    // Faces de l'icosaèdre (indices des sommets)
    std::vector<std::vector<CfemInt>> faces = {
        {0, 11, 5}, {0, 5, 1}, {0, 1, 7}, {0, 7, 10}, {0, 10, 11},
        {1, 5, 9}, {5, 11, 4}, {11, 10, 2}, {10, 7, 6}, {7, 1, 8},
        {3, 9, 4}, {3, 4, 2}, {3, 2, 6}, {3, 6, 8}, {3, 8, 9},
        {4, 9, 5}, {2, 4, 11}, {6, 2, 10}, {8, 6, 7}, {9, 8, 1}
    };

    // Utilisation d'un cache pour ne pas créer deux fois le même nœud milieu
    std::map<std::pair<CfemInt, CfemInt>, CfemInt> midpoint_cache;
    auto getMidpoint = [&](CfemInt v1, CfemInt v2) {
        if (v1 > v2) std::swap(v1, v2);
        auto edge = std::make_pair(v1, v2);
        if (midpoint_cache.count(edge)) return midpoint_cache[edge];

        const auto& p1 = mesh.getVertex(v1).coords;
        const auto& p2 = mesh.getVertex(v2).coords;
        cfem::Vector3D mid = center + ((p1 + p2) * 0.5 - center).normalized() * radius;
        
        CfemInt new_id = mesh.addVertex(mid.x, mid.y, mid.z);
        return midpoint_cache[edge] = new_id;
    };

    for (int s = 0; s < subdivisions; ++s) {
        std::vector<std::vector<CfemInt>> new_faces;
        for (const auto& f : faces) {
            CfemInt a = getMidpoint(f[0], f[1]);
            CfemInt b = getMidpoint(f[1], f[2]);
            CfemInt c = getMidpoint(f[2], f[0]);

            new_faces.push_back({f[0], a, c});
            new_faces.push_back({f[1], b, a});
            new_faces.push_back({f[2], c, b});
            new_faces.push_back({a, b, c});
        }
        faces = std::move(new_faces);
        midpoint_cache.clear();
    }

    // Ajout des cellules au maillage
    for (const auto& f : faces) {
        mesh.addCell(CellType::TRI3, f);
    }

    mesh.buildMesh();

    // mesh.relaxMesh(10);

    if (order == MeshOrder::Quadratic) {
        mesh.upgradeToQuadratic();
        // IMPORTANT : Projeter TOUS les nœuds sur la sphère, pas seulement le bord
        for (CfemInt i = 0; i < mesh.getNumVertices(); ++i) {
            auto& v = mesh.getVertex(i).coords;
            v = center + (v - center).normalized() * radius;
        }
    }
    
    mesh.updateAffineFlags();
}

} // namespace cfem