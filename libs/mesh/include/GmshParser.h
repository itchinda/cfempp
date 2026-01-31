#ifndef CFEM_GMSH_PARSER_H
#define CFEM_GMSH_PARSER_H

#include "Mesh.h"
#include <fstream>
#include <map>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>

static inline void clean(std::string& s) {
    while(!s.empty() && (s.back()=='\r' || s.back()=='\n' || s.back()==' ' || s.back()=='\t'))
        s.pop_back();
}

class GmshParser {
public:
    struct GMSHElemData {
        int type;
        int entityTag;
        std::vector<int> gmshNodeIds;
    };

    static bool read(Mesh& mesh, const std::string& filename) {
        std::ifstream file(filename);
        if (!file.is_open()) {
            std::cerr << "Error: Could not open file " << filename << std::endl;
            return false;
        }

        mesh.clear();
        std::map<int, std::string> physicalNames; 
        std::map<int, std::map<int, std::vector<int>>> entityToPhysical; 
        std::map<int, int> nodeIndexMap; 
        std::vector<GMSHElemData> rawElements;
        int maxDim = 0;

        std::string line;
        while (std::getline(file, line)) {
            clean(line);
            if (line.empty()) continue;

            if (line == "$PhysicalNames") {
                parsePhysicalNames(file, physicalNames);
            } else if (line == "$Entities") {
                parseEntities(file, entityToPhysical);
            } else if (line == "$Nodes") {
                parseNodes(file, mesh, nodeIndexMap);
            } else if (line == "$Elements") {
                parseElements(file, rawElements, maxDim);
            }
        }

        // 1. Add Cells to the Mesh
        for (const auto& e : rawElements) {
            if (getDimFromGmshType(e.type) == maxDim) {
                std::vector<CfemInt> vids;
                for (int gid : e.gmshNodeIds) vids.push_back(nodeIndexMap.at(gid));
                mesh.addCell(mapGmshType(e.type), vids);
            }
        }

        if (mesh.getNumCells() == 0) {
            std::cerr << "[CFEM ERROR] No cells found. Check Physical Surface/Volume in GMSH." << std::endl;
            return false;
        }

        mesh.buildMesh();
        if (mesh.getOrder() == 2) mesh.upgradeToQuadratic();

        // // 2. Prepare MeshEntities
        // std::map<std::string, CfemInt> nameToEntityId;
        // for (auto const& [tag, name] : physicalNames) {
        //     // Find the dimension of this physical group from the physicalNames map
        //     // In GMSH 4.1, the physicalNames section already tells us the dim
        //     int dimOfGroup = -1;
        //     for(auto const& [t, n] : physicalNames) if(n == name) { /* find dim logic */ }
            
        //     // For simplicity, we create the entity; buildMesh handles the actual storage
        //     nameToEntityId[name] = mesh.createEntity(name, 0); // Dim will be updated below
        // }

        // // 3. Generalized Linker: Maps Elements to Physical Entities
        // for (const auto& e : rawElements) {
        //     int dim = getDimFromGmshType(e.type);
            
        //     if (entityToPhysical[dim].count(e.entityTag)) {
        //         std::vector<CfemInt> vids;
        //         for (int gid : e.gmshNodeIds) vids.push_back(nodeIndexMap.at(gid));

        //         for (int pTag : entityToPhysical[dim][e.entityTag]) {
        //             std::string pName = physicalNames[pTag];
        //             MeshEntity* entity = mesh.getEntity(nameToEntityId[pName]);
        //             entity->setDimension(dim);

        //             if (dim == maxDim) {
        //                 // Interior Entity (Cells)
        //                 CfemInt c_idx = findCellByVertices(mesh, vids);
        //                 if (c_idx != -1) entity->addCell(c_idx);
        //             } 
        //             else if (dim == maxDim - 1) {
        //                 // Boundary Entity (Facets)
        //                 CfemInt f_idx = findFacetByVertices(mesh, vids);
        //                 if (f_idx != -1) entity->addFacet(f_idx);
        //             }

        //             // Always add nodes to the entity
        //             for (CfemInt v : vids) entity->addVertex(v);
        //         }
        //     }
        // }
        return true;
    }

private:
    static void consumeFooter(std::ifstream& f, const std::string& footer) {
        std::string line;
        while (std::getline(f, line)) {
            clean(line);
            if (line == footer) break;
        }
    }

    static void parsePhysicalNames(std::ifstream& f, std::map<int, std::string>& names) {
        int count; f >> count;
        for (int i = 0; i < count; ++i) {
            int dim, tag; std::string name;
            f >> dim >> tag >> name;
            if (name.front() == '"') name = name.substr(1, name.size() - 2);
            names[tag] = name;
        }
        consumeFooter(f, "$EndPhysicalNames");
    }

    static void parseEntities(std::ifstream& f, std::map<int, std::map<int, std::vector<int>>>& entMap) {
        size_t counts[4]; // pts, curves, surfs, vols
        f >> counts[0] >> counts[1] >> counts[2] >> counts[3];
        for (int dim = 0; dim < 4; ++dim) {
            for (size_t i = 0; i < counts[dim]; ++i) {
                int tag; double x, y, z; f >> tag >> x >> y >> z;
                if (dim > 0) { double X, Y, Z; f >> X >> Y >> Z; }
                size_t nPhys; f >> nPhys;
                for (size_t p = 0; p < nPhys; ++p) {
                    int pTag; f >> pTag;
                    entMap[dim][tag].push_back(pTag);
                }
                if (dim > 0) { // Skip boundary tags
                    size_t nBound; f >> nBound;
                    for (size_t b = 0; b < nBound; ++b) { int bTag; f >> bTag; }
                }
            }
        }
        consumeFooter(f, "$EndEntities");
    }

    static void parseNodes(std::ifstream& f, Mesh& mesh, std::map<int, int>& nodeMap) {
        size_t nBlocks, totalNodes, minT, maxT;
        f >> nBlocks >> totalNodes >> minT >> maxT;
        int idx = mesh.getNumVertices();
        for (size_t b = 0; b < nBlocks; ++b) {
            int d, t, p; size_t n; f >> d >> t >> p >> n;
            std::vector<int> tags(n);
            for (size_t i = 0; i < n; ++i) f >> tags[i];
            for (size_t i = 0; i < n; ++i) {
                double x, y, z; f >> x >> y >> z;
                if (p > 0) { double dummy; for(int k=0; k<(d*p); ++k) f >> dummy; }
                mesh.addVertex(x, y, z);
                nodeMap[tags[i]] = idx++;
            }
        }
        consumeFooter(f, "$EndNodes");
    }

    static void parseElements(std::ifstream& f, std::vector<GMSHElemData>& raw, int& maxDim) {
        size_t nBlocks, totalElems, minT, maxT;
        f >> nBlocks >> totalElems >> minT >> maxT;
        for (size_t b = 0; b < nBlocks; ++b) {
            int d, t, type; size_t n; f >> d >> t >> type >> n;
            if (d > maxDim) maxDim = d;
            for (size_t i = 0; i < n; ++i) {
                int tag; f >> tag;
                int nN = getNodesPerType(type);
                std::vector<int> v(nN);
                for (int j = 0; j < nN; ++j) f >> v[j];
                raw.push_back({type, t, v});
            }
        }
        consumeFooter(f, "$EndElements");
    }

    static int getDimFromGmshType(int t) {
        if (t == 1 || t == 8) return 1;
        if (t == 2 || t == 3 || t == 9 || t == 10) return 2;
        if (t == 4 || t == 5 || t == 11 || t == 12) return 3;
        return 0;
    }

    static CellType mapGmshType(int t) {
        switch(t) {
            case 1: return CellType::LINE2; case 8:  return CellType::LINE3;
            case 2: return CellType::TRI3;  case 9:  return CellType::TRI6;
            case 3: return CellType::QUAD4; case 10: return CellType::QUAD9;
            case 4: return CellType::TET4;  case 11: return CellType::TET10;
            case 5: return CellType::HEX8;  case 12: return CellType::HEX27;
            default: return CellType::TRI3;
        }
    }

    static int getNodesPerType(int t) {
        switch(t) {
            case 1: return 2; case 8: return 3;
            case 2: return 3; case 9: return 6;
            case 3: return 4; case 10: return 9;
            case 4: return 4; case 11: return 10;
            case 5: return 8; case 12: return 27;
            default: return 0;
        }
    }

    static CfemInt findFacetByVertices(Mesh& mesh, const std::vector<CfemInt>& vids) {
        for (CfemInt i = 0; i < mesh.getNumFacets(); ++i) {
            const auto& f = mesh.getFacet(i);
            if (f.number_of_vertices < (CfemInt)vids.size()) continue;
            bool match = true;
            for (size_t v_req : vids) {
                bool found = false;
                for (int k = 0; k < f.number_of_vertices; ++k)
                    if (f.vertex_ids[k] == v_req) { found = true; break; }
                if (!found) { match = false; break; }
            }
            if (match) return i;
        }
        return -1;
    }

    static CfemInt findCellByVertices(Mesh& mesh, const std::vector<CfemInt>& vids) {
        for (CfemInt i = 0; i < mesh.getNumCells(); ++i) {
            const auto& c = mesh.getCell(i);
            size_t number_of_vertices {c.vertex_ids.size() };
            if (number_of_vertices < vids.size()) continue;
            bool match = true;
            for (size_t v_req : vids) {
                bool found = false;
                for (int k = 0; k < number_of_vertices; ++k)
                    if (c.vertex_ids[k] == v_req) { found = true; break; }
                if (!found) { match = false; break; }
            }
            if (match) return i;
        }
        return -1;
    }
};

#endif