/**
 * @file GmshParser.h
 * @author itchinda
 * @brief CFEM Parser for gmsh (.msh) files.
 * @version 0.2
 * @date 2026-01-31
 * 
 * @copyright Copyright (c) 2026
 * 
 */



#ifndef CFEM_GMSH_PARSER_H
#define CFEM_GMSH_PARSER_H

#include "Mesh.h"
#include <fstream>
#include <map>
#include <unordered_map>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

/**
 * @brief Utilitaire pour le hachage de vecteurs (permet d'utiliser unordered_map avec des listes de sommets)
 * 
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
 * @class GmshParser
 * @brief CFEM parser for Gmsh (.msh) file format version 4.1.
 *   This class handles the translation of Gmsh geometric data into the internal CFEM Mesh structure.
 *   It is optimized for large-scale meshes by utilizing:
 *    - Tag Interning: Combining multiple physical tags into a single integer ID to reduce memory.
 *    - Spatial Hashing: Using an unordered_map with a custom vertex-based hash to match boundary elements
 *      in O(1) time, avoiding expensive O(N) geometric searches.
 *    - Entity Mapping: Automatic population of MeshEntity containers used to define Physical Groups 
 *      (e.g., Boundary Conditions on facets or Material Regions on cells).
 */

class GmshParser {
public:
    /**
     * @brief Raw data structure to temporarily hold Gmsh element info during parsing.
     */
    struct GMSHElemData {
        int type;
        int entityTag;
        std::vector<int> gmshNodeIds;
    };

    bool load(Mesh& mesh, const std::string& filename);

private:
    static void parsePhysicalNames(std::ifstream& f, std::map<int, std::string>& names);
    static void parseEntities(std::ifstream& f, std::map<int, std::map<int, std::vector<int>>>& entMap);
    void parseNodes(std::ifstream& f, Mesh& mesh, std::map<int, int>& nodeMap);
    void parseElements(std::ifstream& f, std::vector<GMSHElemData>& raw, int& maxDim);

};

#endif