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

#include "MeshConformityChecker.h"
#include <queue>
#include <algorithm>
#include <map>
#include <cmath>

namespace cfem {

std::vector<CfemInt> MeshConformityChecker::getCanonicalFaceKey(const LocalFace& local_face) {
    std::vector<CfemInt> key;
    key.reserve(local_face.size);
    for (CfemInt i = 0; i < local_face.size; ++i) {
        key.push_back(local_face.v[i]);
    }
    // Sorting allows us to match faces regardless of node ordering/orientation
    std::sort(key.begin(), key.end());
    return key;
}

ConformityReport MeshConformityChecker::check(const Mesh& mesh, bool check_geometry) {
    ConformityReport report;
    std::ostringstream detail_stream;
    
    const size_t num_cells = mesh.getNumCells();
    const size_t num_nodes = mesh.getNumVertices();

    if (num_cells == 0) {
        report.is_valid = false;
        report.details = " - [CRITICAL] Mesh contains zero cells.\n";
        return report;
    }

    // --- 1. REACHABILITY CHECK (BFS) ---
    // This detects if the mesh is topologically "zipped" together.
    std::vector<bool> visited(num_cells, false);
    std::queue<CfemInt> q;

    q.push(0); // Start walk from the first cell
    visited[0] = true;
    size_t reached_count = 0;

    while (!q.empty()) {
        CfemInt curr = q.front();
        q.pop();
        reached_count++;

        // Access the "Big Block" adjacency via the span-based accessor
        for (CfemInt neighbor : mesh.getCellNeighbors(curr)) {
            if (neighbor != -1 && !visited[neighbor]) {
                visited[neighbor] = true;
                q.push(neighbor);
            }
        }
    }

    if (reached_count < num_cells) {
        report.is_valid = false;
        report.unreachable_cells = static_cast<CfemInt>(num_cells - reached_count);
        detail_stream << " - [TOPOLOGY] Found " << report.unreachable_cells 
                      << " unreachable cells. Parts of the mesh are disconnected.\n";
    }

    // --- 2. ELEMENT AUDIT (Nodes, Faces, Geometry) ---
    std::vector<bool> node_used(num_nodes, false);
    std::map<std::vector<CfemInt>, int> face_usage_count;

    for (CfemInt c_id = 0; c_id < (CfemInt)num_cells; ++c_id) {
        const auto& cell = mesh.getCell(c_id);
        auto conn = mesh.getCellConnectivity(c_id);

        // Mark nodes as used
        for (CfemInt vid : conn) {
            if (vid >= 0 && vid < (CfemInt)num_nodes) {
                node_used[vid] = true;
            }
        }

        // Tally local faces for Manifold check
        CfemInt local_faces_count = Cell::getFaceCount(cell.type);
        for (CfemInt f_idx = 0; f_idx < local_faces_count; ++f_idx) {
            LocalFace lf = mesh.getLocalFaceData(c_id, f_idx);
            face_usage_count[getCanonicalFaceKey(lf)]++;
        }

        // Geometric quality check
        if (check_geometry) {
            CfemReal measure = mesh.computeCellMeasure(cell.type, conn);
            // Using 10^-12 as a robust epsilon for double precision
            if (measure < -1e-12) {
                report.inverted_elements++;
                report.is_valid = false;
                detail_stream << " - [GEOMETRY] Inverted cell detected at ID: " << c_id << "\n";
            } else if (std::abs(measure) <= 1e-12) {
                report.degenerate_elements++;
                report.is_valid = false;
                detail_stream << " - [GEOMETRY] Degenerate (zero-volume) cell at ID: " << c_id << "\n";
            }
        }
    }

    // --- 3. FINAL TALLYING ---
    
    // Check for orphan nodes
    for (CfemInt vid = 0; vid < (CfemInt)num_nodes; ++vid) {
        if (!node_used[vid]) {
            report.isolated_vertices++;
            report.is_valid = false;
        }
    }
    if (report.isolated_vertices > 0) {
        detail_stream << " - [NODAL] Found " << report.isolated_vertices << " isolated vertices.\n";
    }

    // Check for non-manifold edges/faces (shared by > 2 elements)
    for (const auto& [key, count] : face_usage_count) {
        if (count > 2) {
            report.non_manifold_faces++;
            report.is_valid = false;
        }
    }
    if (report.non_manifold_faces > 0) {
        detail_stream << " - [TOPOLOGY] Found " << report.non_manifold_faces 
                      << " non-manifold faces (shared by more than 2 cells).\n";
    }

    report.details = detail_stream.str();
    return report;
}

} // namespace cfem