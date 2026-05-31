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

/**
 * @file GmshParser.h
 * @author itchinda
 * @brief Gmsh 4.1 MSH format parser for the CFEM++ solver.
 * @version 0.2
 * @date 2026-01-31
 * * ### MESH ARCHITECTURE & ENTITY STRATEGY
 * This parser implements a strict decoupling between geometric data and computational 
 * semantics using a dual-layer entity system:
 * * #### 1. Mandatory "Universe" Domain (ID -1)
 * - The solver requires a Physical Group named **"Universe"** in the .msh file.
 * - This entity must match the maximum dimension of the mesh (e.g., 3D for volumes).
 * - **Filtering:** Elements not belonging to any Physical Group are ignored to 
 * ensure a clean computational domain.
 * * #### 2. Automated System Boundary (ID -2)
 * - After topological reconstruction (`buildMesh`), the parser identifies all 
 * facets without neighbors (exterior faces).
 * - These facets are automatically stored in the **"SYS_BOUNDARY"** entity.
 * - This ensures boundary coverage even if the user fails to tag surfaces in Gmsh.
 * * #### 3. User-Defined Entities (ID > 0)
 * - Other Physical Groups (e.g., "Inlet", "Steel_Region") are mapped to positive internal IDs.
 * - **Robust Matching:** Uses a "Corner-only" vertex intersection logic. A facet is 
 * assigned to a group only if its primary geometric vertices (corners) are tagged, 
 * making it immune to diagonal flips or mismatched mid-side nodes in quadratic meshes.
 * * ### OPTIMIZATION TECHNIQUES
 * - **Tag Interning:** Multiple physical tags per element are reduced to a single 
 * integer index (`tag_list_id`) pointing to a persistent table in the Mesh class.
 * - **Vertex Harvesting:** Each entity collects and unique-ifies its member vertices 
 * via a `finalize()` call, providing O(1) access to DoFs for Dirichlet conditions.
 * - **Memory Reservation:** Pre-allocates vectors for nodes and elements based on 
 * Gmsh block headers to minimize reallocations.
 * * @copyright Copyright (c) 2026
 */

#ifndef CFEM_GMSH_PARSER_H
#define CFEM_GMSH_PARSER_H

#include <algorithm>
#include <fstream>
#include <map>
#include <memory>
#include <optional>
#include <set>
#include <string>
#include <unordered_map>
#include <vector>

#include "cfemutils.h"
#include "Vector3D.h"

namespace cfem{
/**
 * @brief Hash utility for vertex vectors.
 * Enables unique identification of facets via spatial hashing in O(N).
 */
struct VectorHash {
    size_t operator()(const std::vector<CfemInt>& v) const {
        size_t seed = v.size();
        for (auto x : v) {
            seed ^= std::hash<CfemInt>{}(x) + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
        return seed;
    }
};

/**
 * @brief Permutation table for Gmsh HEX27 to CFEM internal ordering.
 */
static const CfemInt GMSH_TO_CFEM_HEX27[] = {
        0, 1, 2, 3, 4, 5, 6, 7, 8,   // 0-8 (No change)
        11, 13, 9,                   // 9-11 
        16, 18, 19, 17,              // 12-15 
        10, 12, 14, 15,              // 16-19 
        22, 23, 21, 24, 20, 25,  // 20-25
        26                           // 26 (No change)
    };

class Mesh;

class GmshParser {
public:
    /**
     * @brief Raw data structure for transient element storage during parsing.
     */
    struct GMSHElemData {
        CfemInt type;           ///< Gmsh element type ID
        CfemInt entityTag;      ///< Associated Gmsh entity tag
        std::vector<int> gmshNodeIds; ///< Raw Gmsh node IDs
    };

    /**
     * @brief Primary entry point for loading .msh files.
     * @param filename Path to the .msh file.
     * @return A shared pointer to the Mesh if successful, raise a CFEM_ERROR otherwise.
     */
    std::shared_ptr<Mesh> load(const std::string& filename, std::optional<Vector3D> pReferenceNormal = std::nullopt);

private:
    static void parsePhysicalNames(std::ifstream& f, std::map<std::pair<int, int>, std::string>& names);
    static void parseEntities(std::ifstream& f, std::map<int, std::map<int, std::vector<int>>>& entMap);
    void parseNodes(std::ifstream& f, const std::shared_ptr<Mesh>& mesh, std::map<int, int>& nodeMap);
    void parseElements(std::ifstream& f, std::vector<GMSHElemData>& raw, int& maxDim);
};

}
#endif