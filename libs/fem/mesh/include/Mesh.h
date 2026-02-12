/**
 * @file Mesh.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.2
 * @date 2026-01-31
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#ifndef CFEMPP_MESH
#define CFEMPP_MESH

#include <iostream>
#include <vector>
#include <array>
#include <string>
#include <unordered_set>
#include <stdexcept>
#include <memory>
#include <utility>
#include <map>
#include <set>

#include "cfemutils.h"
#include "cfem_topology.h"
#include "MeshEntity.h"
#include "Vector3D.h"
#include "DynamicVector.h"


namespace cfem{

class FESpace;

struct TagGroup {
    std::vector<int> gmsh_tags;
    int internal_id;
};


enum class MeshOrder{ None, Linear = 1, Quadratic = 2};

class Mesh {
private:
    struct ConstructorToken {};

public:
    explicit Mesh(ConstructorToken);
    explicit Mesh(ConstructorToken, CfemInt p_meshDim, MeshOrder p_meshOrder);

    Mesh(const Mesh&) = delete;
    Mesh& operator=(const Mesh&) = delete;

    Mesh(Mesh&&) = default;
    Mesh& operator=(Mesh&&) = default;

    /**
     * @brief Factory pour un maillage vide
     */
    static std::shared_ptr<Mesh> create() {
        return std::make_shared<Mesh>(ConstructorToken{});
    }

    /**
     * @brief Factory pour un maillage paramétré
     */
    static std::shared_ptr<Mesh> create(CfemInt dim, MeshOrder order) {
        auto mesh = std::make_shared<Mesh>(ConstructorToken{}, dim, order);
        return mesh;
    }

    CfemInt addVertex(CfemReal x,
                      CfemReal y = 0.0,
                      CfemReal z = 0.0);

    CfemInt addCell(CellType type,
                    const std::vector<CfemInt>& v_ids);

    void buildInverseConnectivities();
    void buildFacets();

    void buildMesh();

    CfemInt getDimension() const { return m_meshDim; }
    MeshOrder getOrder() const { return m_meshOrder; }
    bool isQuadratic() const { return (m_meshOrder == MeshOrder::Quadratic); }
    bool isLinear() const { return (m_meshOrder == MeshOrder::Linear); }
    size_t getNumVertices() const { return m_vertices.size(); }
    size_t getNumCells()    const { return m_cells.size(); }
    size_t getNumFacets()   const { return m_facets.size(); }

    const Vertex& getVertex(CfemInt id) const { return m_vertices.at(id); }
    Vertex& getVertex(CfemInt vid) {  return m_vertices.at(vid); }
    const Cell&   getCell(CfemInt id)   const { return m_cells.at(id); }
    Cell&   getCell(CfemInt id) { return m_cells.at(id); }
    const Facet&  getFacet(CfemInt id)  const { return m_facets.at(id); }
    Facet&  getFacet(CfemInt id) { return m_facets.at(id); }


    void getElementCoordinates(CfemInt elementId, cfem::DynamicVector<cfem::Vector3D>& localCoords) const;

    const std::vector<Vertex>& getVertices() const { return m_vertices; }
    std::vector<Vertex>& getVertices() { return m_vertices; }
    const std::vector<Cell>&   getCells()   const { return m_cells; }
    std::vector<Cell>&   getCells()  { return m_cells; }
    const std::vector<Facet>&  getFacets()  const { return m_facets; }
    std::vector<Facet>&  getFacets() { return m_facets; }
    const std::vector<std::unique_ptr<MeshEntity>>&  getEntities()  const { return m_entities; }
    std::vector<std::unique_ptr<MeshEntity>>&  getEntities() { return m_entities; }

    void setOrder(MeshOrder newOrder);
    void setDimension(int newDim);

    CfemInt createEntity(const std::string& name, CfemInt dim);

    MeshEntity& createInternalEntity(CfemInt id,
                                     const std::string& name,
                                     CfemInt dim);


    MeshEntity* getEntity(const std::string& name);
    MeshEntity* getEntity(CfemInt id);

    std::string info() const;
    std::string getMeshFormat() const;

    void upgradeToQuadratic();
    void downgradeToLinear();
    void updateEntities();

    // HELPERS
    std::string cellTypeToString(CellType type) const;
    CfemReal computeCellMeasure(CellType type, const std::vector<CfemInt>& v_ids) const; 
    CfemReal compute2NodesLineLength(CfemInt p_idx0, CfemInt p_idx1) const;
    CfemReal computeSigned3NodesTriangleArea(CfemInt p_idx0, CfemInt p_idx1, CfemInt p_idx2, cfem::Vector3D &pReferenceNormale) const;
    CfemReal computeSigned4NodesTetraVolume(CfemInt p_idx0, CfemInt p_idx1, CfemInt p_idx2, CfemInt p_idx3) const;
    CfemInt getDimFromType(CellType type) const;
    MeshOrder getCellOrderFromType(CellType type) const;

    CfemReal computeCellQuality(CfemInt cellId) const;

    void clear() {
        m_vertices.clear();
        m_cells.clear();
        m_facets.clear();
        m_entities.clear();
        m_entity_ids.clear();
        // Reset any other internal state if necessary
    }
    
    
    void setTagGroups(const std::vector<std::vector<int>>& groups) { m_tagGroups = groups; }
    const std::vector<int>& getPhysicalTags(int groupId) const { return m_tagGroups.at(groupId); }

    void relaxMesh(int iterations, double lambda = 0.4);

    void addObserver(std::shared_ptr<FESpace> space);




private:
    CfemInt m_meshDim {0};
    MeshOrder m_meshOrder {MeshOrder::None};

    std::vector<Vertex> m_vertices;
    std::vector<Cell>   m_cells;
    std::vector<Facet>  m_facets;

    std::vector<std::unique_ptr<MeshEntity>> m_entities;

    std::unordered_set<CfemInt> m_entity_ids;

    std::vector<std::vector<int>> m_tagGroups; // Table de "interning"

    std::vector<std::weak_ptr<FESpace>> m_observers; // list of FESpace attached to the mesh

    void checkInvariant() const{
        CFEM_INVARIANT(m_meshDim >= 0 && m_meshDim <= 3);
        if (m_meshDim==0 || m_meshOrder == MeshOrder::None){
            CFEM_INVARIANT(m_cells.empty() && m_facets.empty());
        }
    }

    void validateConnectivity(CellType type, const std::vector<CfemInt>& v_ids) const;
    void validateEntityName(const std::string& name) const;

    CfemInt getMidpoint(CfemInt p_edgeFirstVertexId, 
                        CfemInt p_edgeSecondVertexId, 
                        std::map<std::pair<CfemInt, CfemInt>, CfemInt>& edge_to_midpoint);

    void upgradeLine2ToLine3(Cell& pLineCell, std::map<std::pair<CfemInt, CfemInt>, CfemInt>& edge_to_midpoint);
    void upgradeTri3ToTri6(Cell& pTriCell, std::map<std::pair<CfemInt, CfemInt>, CfemInt>& edge_to_midpoint);
    void upgradeQuad4ToQuad9(Cell& pQuadCell, std::map<std::pair<CfemInt, CfemInt>, CfemInt>& edge_to_midpoint);
    void upgradeTetra4ToTetra10(Cell& pTetraCell, std::map<std::pair<CfemInt, CfemInt>, CfemInt>& edge_to_midpoint);
    void upgradePrism6ToPrism18(Cell& pPrismCell, std::map<std::pair<CfemInt, CfemInt>, CfemInt>& edge_to_midpoint);
    void upgradeHex8ToHex27(Cell& pHexCell, std::map<std::pair<CfemInt, CfemInt>, CfemInt>& edge_to_midpoint);

    
    CfemReal computeTriQuality(const std::vector<CfemInt>& v) const;
    CfemReal computeQuadQuality(const std::vector<CfemInt>& v) const;
    CfemReal computeTetQuality(const std::vector<CfemInt>& v) const;

};

class LineMeshBuilder {
public:
    static void build(Mesh& mesh, int nx, const cfem::Vector3D& start, const cfem::Vector3D& end) {
        mesh.clear();
        for (int i = 0; i <= nx; ++i) {
            double t = static_cast<double>(i) / nx;
            cfem::Vector3D p = start * (1.0 - t) + end * t;
            mesh.addVertex(p.x, p.y, p.z);
        }

        for (int i = 0; i < nx; ++i) {
            mesh.addCell(CellType::LINE2, { (CfemInt)i, (CfemInt)i + 1 });
        }

        mesh.buildMesh();
        if (mesh.getOrder() == MeshOrder::Quadratic) mesh.upgradeToQuadratic(); // Becomes LINE3
    }
};

class QuadMeshBuilder {
private:
    struct ConstructorToken{};
    

public:
    enum class ElementType { 
        None, FirstDiagonal, SecondDiagonal, CrissCross 
    };

    struct Corners{
        Vector3D P0 = {0,0,0};
        Vector3D P1 = {1,0,0};
        Vector3D P2 = {1,1,0};
        Vector3D P3 = {0,1,0};
    };

    struct Subdivisions{
        CfemInt Nx {1};
        CfemInt Ny {1};
    };

    struct BoundaryNames {
        std::string bottom = "Bottom";
        std::string top    = "Top";
        std::string left   = "Left";
        std::string right  = "Right";
    };

    QuadMeshBuilder(ConstructorToken){};

    static std::shared_ptr<Mesh> create(MeshOrder order,
                                        const Subdivisions& subDiv,
                                        const Corners& corners,
                                        ElementType elemType = ElementType::None) 
    {
        return QuadMeshBuilder::create(order, subDiv, corners, elemType, BoundaryNames{});
    }

    static std::shared_ptr<Mesh> create(MeshOrder order,
                                        Subdivisions subDiv,
                                        const Corners& corners,
                                        ElementType elemType, 
                                        const BoundaryNames& names ) {
        auto mesh = Mesh::create(); 

        QuadMeshBuilder builder(ConstructorToken{});
        
        builder.build(*mesh, order,subDiv, corners, elemType, names );
        
        return mesh;
    }

private:
    void applyBoundaryTags(Mesh& mesh, int nx, int ny, 
                           const std::vector<CfemInt>& grid, 
                           const BoundaryNames& names);
                           
    void build(Mesh& mesh, 
               MeshOrder order,
               const Subdivisions& subDiv,
               const Corners& corners,
               ElementType elemType,
               const BoundaryNames& names);
};

class CubicMeshBuilder {
private:
    struct ConstructorToken{};
    
public:
    enum class ElementType {
        None,       // Keep as HEX8
        SixTet      // Split each Hex into 6 Tetrahedra (standard Kuhn triangulation)
    };

    struct Subdivisions {
        CfemInt Nx {1};
        CfemInt Ny {1};
        CfemInt Nz {1};
    };

    struct Corners {
        Vector3D P0 {0,0,0};
        Vector3D P1 {1,0,0};
        Vector3D P2 {1,1,0};
        Vector3D P3 {0,1,0};
        Vector3D P4 {0,0,1};
        Vector3D P5 {1,0,1};
        Vector3D P6 {1,1,1};
        Vector3D P7 {0,1,1};
    };

    CubicMeshBuilder(ConstructorToken){}

    static std::shared_ptr<Mesh> create(MeshOrder order,
                                        const Subdivisions& subDiv, 
                                        const Corners& corners,
                                        ElementType elemType = ElementType::None) 
    {
        auto mesh = Mesh::create(); 

        CubicMeshBuilder builder(ConstructorToken{});
        
        builder.build(*mesh, order, subDiv, corners, elemType );
        
        return mesh;
    }

private:

    void build(Mesh& mesh, 
               MeshOrder order, 
               const Subdivisions& subDiv, 
               const Corners& corners,
               ElementType elemType = ElementType::None);

};

class DiskMeshBuilder {
private:
    struct ConstructorToken{};

public:
    enum class ElementType { 
        Quadrilateral, 
        Triangle 
    };

    DiskMeshBuilder(ConstructorToken){}

    static std::shared_ptr<Mesh> create(MeshOrder order, 
                                        CfemReal target_mesh_size,
                                        const cfem::Vector3D& center,
                                        double radius,
                                        ElementType type = ElementType::Triangle) 
    {
        auto mesh = Mesh::create(); 

        DiskMeshBuilder builder(ConstructorToken{});
        
        builder.build(*mesh, order, target_mesh_size, center,  radius, type);
        
        return mesh;
    }

private:

    void build(Mesh& mesh,
               MeshOrder order,
               CfemReal target_mesh_size,
               const cfem::Vector3D& center,
               double radius,
               ElementType type = ElementType::Triangle);

    static void projectBoundaryMidNodes(Mesh& mesh, const cfem::Vector3D& center, double radius);

};

class SphereMeshBuilder {
private:
    struct ConstructorToken{};

public:
    enum class ElementType { Hexahedron, Tetrahedron };

    SphereMeshBuilder(ConstructorToken){}

    static std::shared_ptr<Mesh> create(MeshOrder order, 
                                        CfemReal target_mesh_size,
                                        const cfem::Vector3D& center,
                                        double radius,
                                        ElementType type = ElementType::Tetrahedron) 
    {
        auto mesh = Mesh::create(); 

        SphereMeshBuilder builder(ConstructorToken{});
        
        builder.build(*mesh, order, target_mesh_size, center,  radius, type);
        
        return mesh;
    }

private:
    static void build(Mesh& mesh,
                      MeshOrder order,
                      CfemReal target_mesh_size,
                      const cfem::Vector3D& center,
                      double radius,
                      ElementType type = ElementType::Tetrahedron);

    static void projectBoundaryNodes(Mesh& mesh,
                                     const cfem::Vector3D& center,
                                     double radius);

};


struct QualityReport {
    double minVolume = std::numeric_limits<double>::max();
    double maxVolume = -std::numeric_limits<double>::max();
    int invertedCount = 0;
    int degenerateCount = 0; // Near-zero volume
};

class MeshQualityTool {
public:
    static QualityReport analyze(const Mesh& mesh) {
        QualityReport report;

        for (size_t i = 0; i < mesh.getNumCells(); ++i) {
            const auto& cell = mesh.getCell(i);
            
            // Signed Volume (Jacobian)
            double vol = mesh.computeCellMeasure(cell.type, cell.vertex_ids);
            
            if (vol < -cfem::ZERO_TOLERANCE) {
                report.invertedCount++;
            } else if (std::abs(vol) <= cfem::ZERO_TOLERANCE) {
                report.degenerateCount++;
            }

            report.minVolume = std::min(report.minVolume, vol);
            report.maxVolume = std::max(report.maxVolume, vol);
        }

        return report;
    }

    static void printSummary(const QualityReport& report) {
        std::cout << "\n--- CFEM MESH QUALITY REPORT ---" << std::endl;
        std::cout << "Min Volume: " << report.minVolume << std::endl;
        std::cout << "Max Volume: " << report.maxVolume << std::endl;
        
        if (report.invertedCount > 0) {
            CFEM_WARNING("Mesh contains " + std::to_string(report.invertedCount) + " inverted elements!");
        }
        if (report.degenerateCount > 0) {
            CFEM_WARNING("Mesh contains " + std::to_string(report.degenerateCount) + " degenerate (zero-volume) elements!");
        }
        
        if (report.invertedCount == 0 && report.degenerateCount == 0) {
            std::cout << "Result: Mesh is valid for simulation." << std::endl;
        }
        std::cout << "--------------------------------\n" << std::endl;
    }
};

}  // namespace

#include "Mesh.tpp"

#endif


