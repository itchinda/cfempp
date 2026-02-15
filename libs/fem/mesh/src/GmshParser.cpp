#include "GmshParser.h"
#include <iostream>

namespace cfem{

/**
 * @brief Maps Gmsh element type IDs to CFEM internal CellType.
 */
inline static CellType mapGmshType(int t) {
    switch (t) {
        case 1: return CellType::LINE2; case 8:  return CellType::LINE3;
        case 2: return CellType::TRI3;  case 9:  return CellType::TRI6;
        case 3: return CellType::QUAD4; case 10: return CellType::QUAD9;
        case 4: return CellType::TET4;  case 11: return CellType::TET10;
        case 5: return CellType::HEX8;  case 12: return CellType::HEX27;
        default: return CellType::TRI3;
    }
}

/**
 * @brief Get the Nodes Per Type object
 * 
 * @param t 
 * @return int 
 */
inline static int getNodesPerType(int t) {
    switch (t) {
        case 1: return 2; case 8: return 3;
        case 2: return 3; case 9: return 6;
        case 3: return 4; case 10: return 9;
        case 4: return 4; case 11: return 10;
        case 5: return 8; case 12: return 27;
        default: return 0;
    }
}

/**
 * @brief Get the Dim From Gmsh Type object
 * 
 * @param t GMSH Type 
 * @return int 
 */
inline static int getDimFromGmshType(int t) {
    if (t == 1 || t == 8) return 1;
    if (t == 2 || t == 3 || t == 9 || t == 10) return 2;
    if (t == 4 || t == 5 || t == 11 || t == 12) return 3;
    return 0;
}


/**
 * @brief Returns the number of corner (geometric) nodes for a given type.
 * This is the core of the robust boundary matching logic.
 */
inline static int getNumMainCorners(CellType type) {
    switch (type) {
        
        case CellType::LINE2:
        case CellType::LINE3:
            return 2;

        case CellType::TRI3:
        case CellType::TRI6:
            return 3; 

        case CellType::QUAD4:
        case CellType::QUAD9:
            return 4;

        case CellType::TET4:
        case CellType::TET10:
            return 4;
        
        case CellType::PRISM6:
        case CellType::PRISM18:
            return 6;

        case CellType::HEX8:
        case CellType::HEX27:
            return 4;

        default:
            CFEM_ERROR("Non implemented Cell/Facet type.");
            return 0;
    }
}

inline static void clean(std::string &s) {
    while (!s.empty() && (s.back() == '\r' || s.back() == '\n' || s.back() == ' ' || s.back() == '\t'))
        s.pop_back();
}

inline static void consumeFooter(std::ifstream &f, const std::string &footer) {
    std::string line;
    while (std::getline(f, line)) {
        clean(line);
        if (line == footer) break;
    }
}

 /**
  * @brief Parses the $PhysicalNames section.
  * Maps (dimension, tag) pairs to their string names (e.g., {3, 1} -> "Universe").
  * * This is the first step in identifying user-defined regions.
  * 
  * @param f gmsh file
  * @param names (dimension, tag) pairs to their string names
  */
void GmshParser::parsePhysicalNames(std::ifstream& f, std::map<std::pair<int, int>, std::string>& names) {
    int count; f >> count;
    for (int i = 0; i < count; ++i) {
        int dim, tag; std::string name;
        f >> dim >> tag >> name;
        if (name.front() == '"') name = name.substr(1, name.size() - 2);
        names[{dim, tag}] = name;
    }
    consumeFooter(f, "$EndPhysicalNames");
}

/**
 * @brief Parses the $Entities section.
 * This section establishes the relationship between internal Gmsh 'tags' 
 * and user-defined 'Physical Groups'. 
 * * It populates a map where entMap[dimension][entityTag] returns a list of 
 * associated Physical Group IDs.
 * 
 * @param f gmsh file
 * @param entMap 
 */
void GmshParser::parseEntities(std::ifstream& f, std::map<int, std::map<int, std::vector<int>>>& entMap) {
    size_t counts[4]; f >> counts[0] >> counts[1] >> counts[2] >> counts[3];
    for (int dim = 0; dim < 4; ++dim) {
        for (size_t i = 0; i < counts[dim]; ++i) {
            int tag; double x, y, z; 
            f >> tag >> x >> y >> z;
            if (dim > 0) {
                double X, Y, Z; 
                f >> X >> Y >> Z; 
            }
            size_t nPhys; f >> nPhys;
            for (size_t p = 0; p < nPhys; ++p) {
                int pTag; 
                f >> pTag;
                entMap[dim][tag].push_back(pTag);
            }
            if (dim > 0) {
                size_t nBound; f >> nBound;
                for (size_t b = 0; b < nBound; ++b) { 
                    int bTag; 
                    f >> bTag; 
                }
            }
        }
    }
    consumeFooter(f, "$EndEntities");
}

/**
 * @brief Parses the $Nodes section and populates the Mesh vertices.
 * @param f Input stream.
 * @param mesh The mesh object where vertices are added.
 * @param nodeMap A translation map: nodeMap[GmshID] = InternalIndex.
 * * This method uses reserve() to prevent reallocations and maps 
 * non-sequential Gmsh IDs to a zero-based internal index system.
 */
void GmshParser::parseNodes(std::ifstream& f, const std::shared_ptr<Mesh>& mesh, std::map<int, int>& nodeMap) {
    size_t nBlocks, totalNodes, minT, maxT;
    f >> nBlocks >> totalNodes >> minT >> maxT;
    int idx = (int)mesh->getNumVertices();
    mesh->getVertices().reserve(idx + totalNodes);
    for (size_t b = 0; b < nBlocks; ++b) {
        int d, t, p; 
        size_t n; 
        f >> d >> t >> p >> n;
        std::vector<int> tags(n);
        for (size_t i = 0; i < n; ++i){
            f >> tags[i];
        }
        for (size_t i = 0; i < n; ++i) {
            double x, y, z; 
            f >> x >> y >> z;
            mesh->addVertex(x, y, z);
            nodeMap[tags[i]] = idx++;
        }
    }
    consumeFooter(f, "$EndNodes");
}

/**
 * @brief Parses the $Elements section into a temporary buffer.
 * @param f Input stream.
 * @param raw Vector of GMSHElemData to store element types and connectivity.
 * @param maxDim Output parameter to capture the highest dimension found (2D vs 3D).
 * * Elements are stored in a raw format first to allow the Mesh to build its 
 * topology (facets) before groups are distributed.
 */
void GmshParser::parseElements(std::ifstream& f, std::vector<GMSHElemData>& raw, int& maxDim) {
    size_t nBlocks, totalElems, minT, maxT;
    f >> nBlocks >> totalElems >> minT >> maxT;
    raw.reserve(raw.size() + totalElems);
    for (size_t b = 0; b < nBlocks; ++b) {
        int d, t, type; 
        size_t n; 
        f >> d >> t >> type >> n;
        if (d > maxDim) maxDim = d;
        for (size_t i = 0; i < n; ++i) {
            int tag; f >> tag;
            int nN = getNodesPerType(type);
            std::vector<int> v(nN);
            for (int j = 0; j < nN; ++j) {
                f >> v[j];
            }
            raw.push_back({type, t, v});
        }
    }
    consumeFooter(f, "$EndElements");
}


std::shared_ptr<Mesh> GmshParser::load(const std::string &filename) 
{
    auto mesh = Mesh::create();
    // Parsing strategy:
    // 1. I/O Initialization
    // 2. Section Parsing ($PhysicalNames, $Entities, $Nodes, $Elements)
    // 3. ID Mapping & Universe Setup
    
    // 4. Volume/Boundary Processing
    //   - High-order permutations (TET10/HEX27)
    //   - Corner Harvesting: boundary group tags are stored per-node (corners only)
    
    // 5. Topological Reconstruction
    //   - calls mesh->buildMesh() to generate internal facets and neighbors.
    
    // 6. Corner-Based Facet Distribution
    //   - For each boundary facet, check if its N primary corner nodes
    //     share a Physical Group tag from Gmsh.
    
    
    // 7. Finalization
    //   - finalize() unique-ifies vertices per entity for fast DoF access.


    std::ifstream file(filename);
    if (!file.is_open()) {
        CFEM_ERROR( ("Failed to open " + filename + ". Make sure the file is in the specified path or is not corrupted.").c_str() );
        return nullptr;
    }

    mesh->clear();
    std::map<std::pair<int, int>, std::string> physicalNames;
    std::map<int, std::map<int, std::vector<int>>> entityToPhysical;
    std::map<int, int> nodeIndexMap;
    std::vector<GMSHElemData> rawElements;
    int maxDim = 0;

    // Parsing sequential sections
    std::string line;
    while (std::getline(file, line)) {
        if (line.find("$PhysicalNames") != std::string::npos) parsePhysicalNames(file, physicalNames);
        else if (line.find("$Entities") != std::string::npos) parseEntities(file, entityToPhysical);
        else if (line.find("$Nodes") != std::string::npos) parseNodes(file, mesh, nodeIndexMap);
        else if (line.find("$Elements") != std::string::npos) parseElements(file, rawElements, maxDim);
    }

    mesh->setDimension(maxDim);

    // Mapping IDs Internes & Validation Universe
    std::map<std::pair<int, int>, int> gmshToInternal;
    bool universeFound = false;

    for (auto const& [key, name] : physicalNames) {
        if (name == "Universe" && key.first == maxDim) {
            gmshToInternal[key] = UNIVERSE_ID;
            universeFound = true;
            mesh->createInternalEntity(UNIVERSE_ID, name, maxDim);
        } else {
            gmshToInternal[key] = mesh->createEntity(name, key.first);
        }
    }

    if (!universeFound) {
        CFEM_ERROR("Fatal Error: Physical entity 'Universe' not found.");
    }

    // Tag Interning & Corner-based Mapping
    std::vector<std::vector<int>> tagGroups;
    std::map<std::vector<int>, int> comboToId;
    std::map<CfemInt, std::vector<int>> cornerToBoundaryGroups;

    auto getTagGroupId = [&](int entityTag, int dim) -> int {
        if (entityToPhysical[dim].count(entityTag) == 0) return -1;
        std::vector<int> pTags = entityToPhysical[dim][entityTag];
        std::sort(pTags.begin(), pTags.end());
        if (comboToId.find(pTags) == comboToId.end()) {
            comboToId[pTags] = (int)tagGroups.size();
            tagGroups.push_back(pTags);
        }
        return comboToId[pTags];
    };

    // Creation of Cells and Harvest of Corner Vertices
    for (const auto &e : rawElements) {
        int dim = getDimFromGmshType(e.type);
        int gId = getTagGroupId(e.entityTag, dim);
        if (gId == -1) continue;

        std::vector<CfemInt> vids;
        for (int gid : e.gmshNodeIds) vids.push_back(nodeIndexMap.at(gid));

        CellType cfem_cell_type = mapGmshType(e.type);

        if (dim == maxDim) {
            if (mesh->getOrder() == MeshOrder::None) mesh->setOrder(mesh->getCellOrderFromType(cfem_cell_type));
            
            // Permutations for high-order elements (TET10/HEX27) to conform with CFEM++/Paraview ordering standard
            if (cfem_cell_type == CellType::TET10) {
                std::swap(vids[8], vids[9]);
            } else if (cfem_cell_type == CellType::HEX27) {
                std::vector<CfemInt> v_orig = vids;
                for (int i = 0; i < 27; ++i) vids[i] = v_orig[GMSH_TO_ODE_HEX27[i]];
            }

            CfemInt c_idx = mesh->addCell(cfem_cell_type, vids);
            mesh->getCell(c_idx).tag_list_id = gId;
        } 
        else if (dim == maxDim - 1) {
            int numCorners = getNumMainCorners(mapGmshType(e.type));
            for (int i = 0; i < numCorners; ++i) {
                cornerToBoundaryGroups[vids[i]].push_back(gId);
            }
        }
    }

    mesh->buildMesh();
    mesh->setTagGroups(tagGroups);

    auto distribute = [&](int gId, CfemInt elementIdx, int dimElem, const std::vector<CfemInt>& v_ids) {
        if (gId < 0) return;
        for (int pTagGmsh : tagGroups[gId]) {
            int internalId = gmshToInternal[{dimElem, pTagGmsh}];
            MeshEntity *ent = mesh->getEntity(internalId);
            if (dimElem == maxDim) ent->addCell(elementIdx);
            else ent->addFacet(elementIdx);
            for (CfemInt vid : v_ids) ent->addVertex(vid);
        }
    };

    for (CfemInt i = 0; i < (CfemInt)mesh->getNumCells(); ++i) {
        const Cell& c = mesh->getCell(i);
        distribute(c.tag_list_id, i, maxDim, c.vertex_ids);
    }

    for (size_t i = 0; i < mesh->getNumFacets(); ++i) {
        Facet &f = const_cast<Facet &>(mesh->getFacet(i));
        if (!f.isBoundary()) continue;

        int numCorners = getNumMainCorners(f.type);

        std::map<int, int> cornerOccurrence;
        for (int j = 0; j < numCorners; ++j) {
            CfemInt v = f.vertex_ids[j];
            if (cornerToBoundaryGroups.count(v)) {
                std::set<int> uniqueGroups(cornerToBoundaryGroups[v].begin(), cornerToBoundaryGroups[v].end());
                for (int gId : uniqueGroups) cornerOccurrence[gId]++;
            }
        }

        for (auto const& [gId, count] : cornerOccurrence) {
            if (count == numCorners) { 
                distribute(gId, (CfemInt)i, maxDim - 1, f.vertex_ids);
            }
        }

        MeshEntity* sysBound = mesh->getEntity(SYS_BOUNDARY_ID);
        sysBound->addFacet((CfemInt)i);
        for (CfemInt vid : f.vertex_ids) sysBound->addVertex(vid);
    }

    // finalize()
    for (auto const& [key, internalId] : gmshToInternal) {
        mesh->getEntity(internalId)->finalize();
    }
    mesh->getEntity(SYS_BOUNDARY_ID)->finalize();

    return mesh;
}

} //cfem