#include "GmshParser.h"

/**
 * @brief 
 * 
 * @param t 
 * @return CellType 
 */
inline static CellType mapGmshType(int t)
{
    switch (t)
    {
    case 1:
        return CellType::LINE2;
    case 8:
        return CellType::LINE3;
    case 2:
        return CellType::TRI3;
    case 9:
        return CellType::TRI6;
    case 3:
        return CellType::QUAD4;
    case 10:
        return CellType::QUAD9;
    case 4:
        return CellType::TET4;
    case 11:
        return CellType::TET10;
    case 5:
        return CellType::HEX8;
    case 12:
        return CellType::HEX27;
    default:
        return CellType::TRI3;
    }
}

/**
 * @brief Get the Nodes Per Type object
 * 
 * @param t 
 * @return int 
 */
inline static int getNodesPerType(int t)
{
    switch (t)
    {
    case 1:
        return 2;
    case 8:
        return 3;
    case 2:
        return 3;
    case 9:
        return 6;
    case 3:
        return 4;
    case 10:
        return 9;
    case 4:
        return 4;
    case 11:
        return 10;
    case 5:
        return 8;
    case 12:
        return 27;
    default:
        return 0;
    }
}



/** 
 * @brief 
 * 
 * @param s 
 */
inline static void clean(std::string &s)
{
    while (!s.empty() && (s.back() == '\r' || s.back() == '\n' || s.back() == ' ' || s.back() == '\t'))
        s.pop_back();
}

/**
 * @brief 
 * 
 * @param f 
 * @param footer 
 */
inline static void consumeFooter(std::ifstream &f, const std::string &footer)
{
    std::string line;
    while (std::getline(f, line))
    {
        clean(line);
        if (line == footer)
            break;
    }
}




/**
 * @brief Parses the $PhysicalNames section to map human-readable names to integer tags.
 * 
 * @param f 
 * @param names 
 */
void GmshParser::parsePhysicalNames(std::ifstream& f, std::map<int, std::string>& names) 
{
    int count; f >> count;
    for (int i = 0; i < count; ++i) {
        int dim, tag; std::string name;
        f >> dim >> tag >> name;
        if (name.front() == '"') name = name.substr(1, name.size() - 2);
        names[tag] = name;
    }
    consumeFooter(f, "$EndPhysicalNames");

}

/**
 * @brief Parses the $Entities section to establish the hierarchy between geometry and physical groups.
 * 
 * @param f 
 * @param entMap 
 */
void GmshParser::parseEntities(std::ifstream& f, std::map<int, std::map<int, std::vector<int>>>& entMap) 
{
    size_t counts[4]; 
    f >> counts[0] >> counts[1] >> counts[2] >> counts[3];
    for (int dim = 0; dim < 4; ++dim) {
        for (size_t i = 0; i < counts[dim]; ++i) {
            int tag; 
            double x, y, z; 
            f >> tag >> x >> y >> z;
            if (dim > 0) { 
                double X, Y, Z; 
                f >> X >> Y >> Z; 
            }
            size_t nPhys; f >> nPhys;
            for (size_t p = 0; p < nPhys; ++p) {
                int pTag; f >> pTag;
                entMap[dim][tag].push_back(pTag);
            }
            if (dim > 0) {
                size_t nBound; f >> nBound;
                for (size_t b = 0; b < nBound; ++b) { 
                    int bTag; f >> bTag; 
                }
            }
        }
    }
    consumeFooter(f, "$EndEntities");
}

/**
 * @brief 
 * 
 * @param f 
 * @param mesh 
 * @param nodeMap 
 */
void GmshParser::parseNodes(std::ifstream& f, Mesh& mesh, std::map<int, int>& nodeMap) 
{
    size_t nBlocks, totalNodes, minT, maxT;
    f >> nBlocks >> totalNodes >> minT >> maxT;
    int idx = (int)mesh.getNumVertices();
    for (size_t b = 0; b < nBlocks; ++b) {
        int d, t, p;
         size_t n; 
         f >> d >> t >> p >> n;
        std::vector<int> tags(n);
        for (size_t i = 0; i < n; ++i) 
            f >> tags[i];
        for (size_t i = 0; i < n; ++i) {
            double x, y, z; 
            f >> x >> y >> z;
            mesh.addVertex(x, y, z);
            nodeMap[tags[i]] = idx++;
        }
    }
    consumeFooter(f, "$EndNodes");
}

/**
 * @brief 
 * 
 * @param f 
 * @param raw 
 * @param maxDim 
 */
void GmshParser::parseElements(std::ifstream& f, std::vector<GMSHElemData>& raw, int& maxDim) 
{
    size_t nBlocks, totalElems, minT, maxT;
    f >> nBlocks >> totalElems >> minT >> maxT;

    raw.reserve(raw.size() + totalElems);

    for (size_t b = 0; b < nBlocks; ++b) {
        int d, t, type; 
        size_t n; 
        f >> d >> t >> type >> n;
        if (d > maxDim) {
            maxDim = d;
        }
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

/**
 * @brief Get the Dim From Gmsh Type object
 * 
 * @param t 
 * @return int 
 */
inline static int getDimFromGmshType(int t)
{
    if (t == 1 || t == 8)
        return 1;
    if (t == 2 || t == 3 || t == 9 || t == 10)
        return 2;
    if (t == 4 || t == 5 || t == 11 || t == 12)
        return 3;
    return 0;
}


/**
 * @brief Main entry point to read a Gmsh file and populate a Mesh object.
 * @param mesh Reference to the Mesh object to be populated.
 * @param filename Path to the .msh file.
 * @return true If the mesh was successfully loaded and built.
 * @return false If the file could not be opened or is malformed.
 * @note This function performs an automatic 'buildMesh()' and can upgrade to quadratic 
 *       if the Mesh object's order is set to 2.
 */
bool GmshParser::load(Mesh &mesh, const std::string &filename)
{
    std::ifstream file(filename);
    if (!file.is_open())
    {
        std::cerr << "Error: Could not open file " << filename << std::endl;
        return false;
    }

    mesh.clear();
    std::map<int, std::string> physicalNames;
    std::map<int, std::map<int, std::vector<int>>> entityToPhysical;
    std::map<int, int> nodeIndexMap;
    std::vector<GMSHElemData> rawElements;
    int maxDim = 0;

    // *** 1. LECTURE DU FICHIER ***
    std::string line;
    while (std::getline(file, line))
    {
        if (line.find("$PhysicalNames") != std::string::npos)
            parsePhysicalNames(file, physicalNames);
        else if (line.find("$Entities") != std::string::npos)
            parseEntities(file, entityToPhysical);
        else if (line.find("$Nodes") != std::string::npos)
            parseNodes(file, mesh, nodeIndexMap);
        else if (line.find("$Elements") != std::string::npos)
            parseElements(file, rawElements, maxDim);
    }

    // *** 2. GESTION DES GROUPES DE TAGS (INTERNING) ***
    std::vector<std::vector<int>> tagGroups;
    std::map<std::vector<int>, int> comboToId;

    auto getTagGroupId = [&](int entityTag, int dim) -> int
    {
        if (entityToPhysical[dim].find(entityTag) == entityToPhysical[dim].end())
            return -1;
        std::vector<int> pTags = entityToPhysical[dim][entityTag];
        std::sort(pTags.begin(), pTags.end());
        if (comboToId.find(pTags) == comboToId.end())
        {
            comboToId[pTags] = (int)tagGroups.size();
            tagGroups.push_back(pTags);
        }
        return comboToId[pTags];
    };

    // *** 3. CREATION DES CELLULES ET CARTOGRAPHIE DES FRONTIERES ***
    using FaceMap = std::unordered_map<std::vector<CfemInt>, int, VectorHash>;
    FaceMap boundaryFaceTags;

    for (const auto &e : rawElements)
    {
        int dim = getDimFromGmshType(e.type);
        std::vector<CfemInt> vids;
        for (int gid : e.gmshNodeIds)
            vids.push_back(nodeIndexMap.at(gid));

        if (dim == maxDim)
        {
            CfemInt c_idx = mesh.addCell(mapGmshType(e.type), vids);
            // On injecte directement l'ID de groupe dans la cellule
            mesh.getCell(c_idx).tag_list_id = getTagGroupId(e.entityTag, dim);
        }
        else if (dim == maxDim - 1)
        {
            // On stocke les tags des éléments de surface Gmsh pour un matching ultérieur
            std::vector<CfemInt> sortedV = vids;
            std::sort(sortedV.begin(), sortedV.end());
            boundaryFaceTags[sortedV] = getTagGroupId(e.entityTag, dim);
        }
    }

    // *** 4. CONSTRUCTION DU MAILLAGE (FACETTES & CONNECTIVITÉ) ***
    mesh.buildMesh();

    // *** 5. MATCHING DES FACETTES (O(N) grâce à l'unordered_map) ***
    for (size_t i = 0; i < mesh.getNumFacets(); ++i)
    {
        Facet &f = const_cast<Facet &>(mesh.getFacet(i));
        std::vector<CfemInt> sortedV = f.vertex_ids;
        std::sort(sortedV.begin(), sortedV.end());

        auto it = boundaryFaceTags.find(sortedV);
        if (it != boundaryFaceTags.end())
        {
            f.tag_list_id = it->second;
        }
        else
        {
            f.tag_list_id = f.isBoundary() ? -2 : -1;
        }
    }

    // *** 6. REMPLISSAGE DES ENTITÉS ***
    std::map<int, CfemInt> internalTagToEntity;
    for (auto const &[tag, name] : physicalNames)
    {
        internalTagToEntity[tag] = mesh.createEntity(name, 0);
        //std::cout << "Entity " << name << " created" << std::endl;
    }

    auto distribute = [&](int gId, CfemInt elementIdx, bool isCell)
    {
        if (gId < 0)
            return;
        for (int pTag : tagGroups[gId])
        {
            MeshEntity *ent = mesh.getEntity(internalTagToEntity[pTag]);
            if (isCell)
                ent->addCell(elementIdx);
            else
                ent->addFacet(elementIdx);
        }
    };

    for (CfemInt i = 0; i < (CfemInt)mesh.getNumCells(); ++i)
        distribute(mesh.getCell(i).tag_list_id, i, true);

    for (CfemInt i = 0; i < (CfemInt)mesh.getNumFacets(); ++i)
        distribute(mesh.getFacet(i).tag_list_id, i, false);

    if (mesh.getOrder() == 2)
        mesh.upgradeToQuadratic();

    return true;
}

