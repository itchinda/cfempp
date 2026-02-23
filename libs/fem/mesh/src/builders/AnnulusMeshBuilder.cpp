#include "AnnulusMeshBuilder.h"


namespace cfem
{


    void AnnulusMeshBuilder::build(Mesh &mesh,
                      MeshOrder order,
                      CfemReal h,
                      const Vector3D &center,
                      CfemReal r_in,
                      CfemReal r_out,
                      DiskElementType type,
                      const AnnulusBoundaryNames &names)
    {
        mesh.clear();
        mesh.setDimension(2);
        mesh.setOrder(MeshOrder::Linear);

        // 1. Calcul des subdivisions radiales
        CfemReal thickness = r_out - r_in;
        int num_rings = static_cast<int>(std::round(thickness / h));
        if (num_rings < 1) num_rings = 1;
        CfemReal dr = thickness / num_rings;

        std::vector<std::vector<CfemInt>> ring_nodes(num_rings + 1);

        // 2. Génération des Sommets
        for (int j = 0; j <= num_rings; ++j) {
            CfemReal current_r = r_in + j * dr;
            int num_points = static_cast<int>(std::round(2.0 * M_PI * current_r / h));
            if (num_points < 8) num_points = 8; 

            for (int i = 0; i < num_points; ++i) {
                CfemReal theta = (2.0 * M_PI * i) / num_points;
                CfemReal x = center.x + current_r * std::cos(theta);
                CfemReal y = center.y + current_r * std::sin(theta);
                ring_nodes[j].push_back(mesh.addVertex(x, y, center.z));
            }
        }

        // 3. Connectivité (Triangles CCW)
        auto addCCW = [&](CfemInt v0, CfemInt v1, CfemInt v2) {
            const auto &p0 = mesh.getVertex(v0).coords, &p1 = mesh.getVertex(v1).coords, &p2 = mesh.getVertex(v2).coords;
            if (((p1.x-p0.x)*(p2.y-p0.y)-(p1.y-p0.y)*(p2.x-p0.x)) < 0) {
                mesh.addCell(CellType::TRI3, {v0, v2, v1});
            }else {
                mesh.addCell(CellType::TRI3, {v0, v1, v2});
            }
        };

        for (int j = 0; j < num_rings; ++j) {
            const auto &inner = ring_nodes[j], &outer = ring_nodes[j+1];
            size_t n_in = inner.size(), n_out = outer.size(), k = 0;
            for (size_t i = 0; i < n_in; ++i) {
                size_t next_i = (i + 1) % n_in;
                size_t target_k = static_cast<size_t>((static_cast<double>(next_i)/n_in)*n_out);
                while (k != target_k || (next_i == 0 && k != 0)) {
                    size_t next_k = (k + 1) % n_out;
                    addCCW(inner[i], outer[next_k], outer[k]);
                    k = next_k;
                    if (k == 0 && next_i != 0) break;
                }
                addCCW(inner[i], inner[next_i], outer[k]);
            }
        }

        // 4. Finalisation du maillage P1
        mesh.buildMesh();

        // 5. Marquage des entités (Inspiré de QuadMeshBuilder)
        applyRadialBoundaryTags(mesh, ring_nodes, names);

        // 6. Passage au Quadratique et Projection Sélective
        if (order == MeshOrder::Quadratic) {
            mesh.upgradeToQuadratic();
            projectQuadraticAnnulus(mesh, center, r_in, r_out, names);
        }
    }

    void AnnulusMeshBuilder::applyRadialBoundaryTags(
                                                        Mesh &mesh, 
                                                        const std::vector<std::vector<CfemInt>> &ring_nodes,
                                                        const AnnulusBoundaryNames &names
                                                    ) 
    {
        CfemInt innerId = mesh.createEntity(names.inner, 1);
        CfemInt outerId = mesh.createEntity(names.outer, 1);

        // Map inverse pour identification rapide : vertex_id -> ring_index
        std::unordered_map<CfemInt, int> nodeToRing;
        for(int j=0; j < ring_nodes.size(); ++j)
            for(auto vid : ring_nodes[j]) nodeToRing[vid] = j;

        auto *sysBoundary = mesh.getEntity(cfem::SYS_BOUNDARY_ID);
        if (!sysBoundary) return;

        int lastRing = ring_nodes.size() - 1;

        for (auto f_idx : sysBoundary->getFacetIds()) {
            auto fv = mesh.getFacetConnectivity(f_idx);
            bool isInner = true, isOuter = true;

            for (auto vid : fv) {
                int ring = nodeToRing[vid];
                if (ring != 0) isInner = false;
                if (ring != lastRing) isOuter = false;
            }

            auto assign = [&](CfemInt entId) {
                auto *ent = mesh.getEntity(entId);
                ent->addFacet(f_idx);
                for (auto vid : fv) ent->addVertex(vid);
            };

            if (isInner) assign(innerId);
            else if (isOuter) assign(outerId);
        }
        mesh.getEntity(innerId)->finalize();
        mesh.getEntity(outerId)->finalize();
    }

    void AnnulusMeshBuilder::projectQuadraticAnnulus(
                                                        Mesh &mesh, 
                                                        const Vector3D &center, 
                                                        CfemReal r_in, 
                                                        CfemReal r_out,
                                                        const AnnulusBoundaryNames &names
                                                    ) 
    {
        auto project = [&](const std::string &name, CfemReal radius) {
            auto entity = mesh.getEntity(name);
            if (!entity) return;
            
            for (auto f_idx : entity->getFacetIds()) {
                const Facet &f = mesh.getFacet(f_idx);
                // On récupère le noeud milieu (index 2 sur une face de bord TRI3->6)
                CfemInt midNodeIdx = mesh.getFacetConnectivity(f_idx).back(); 
                auto &node = mesh.getVertex(midNodeIdx);
                
                CfemReal dx = node.coords.x - center.x;
                CfemReal dy = node.coords.y - center.y;
                CfemReal d = std::hypot(dx, dy);
                if (d > ZERO_TOLERANCE) {
                    node.coords.x = center.x + (dx/d) * radius;
                    node.coords.y = center.y + (dy/d) * radius;
                }
            }
        };
        project(names.inner, r_in);
        project(names.outer, r_out);
    }
}