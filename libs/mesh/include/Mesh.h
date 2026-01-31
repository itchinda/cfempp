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
#include "MeshEntity.h"
#include "Vector3D.h"


enum class CellType {
    LINE2, TRI3, QUAD4, TET4, PRISM6, HEX8,
    LINE3, TRI6, QUAD9, TET10, PRISM18, HEX27
};

struct Vertex {
    Vector3D coords;
    std::vector<CfemInt> connected_cells;
    int tag_list_id = -1;
};

struct Cell {
    CellType type;
    std::vector<CfemInt> vertex_ids;
    int tag_list_id = -1;
};

struct Facet {
    CfemInt id {-1};
    std::vector<CfemInt> vertex_ids {};
    CfemInt number_of_vertices {0};
    CfemInt left_cell  {-1};
    CfemInt right_cell {-1};
    int tag_list_id = -1;

    bool isBoundary() const { return right_cell < 0; }
};

class Mesh {
public:
    explicit Mesh(CfemInt p_meshDim, CfemInt p_meshOrder);

    CfemInt addVertex(CfemReal x,
                      CfemReal y = 0.0,
                      CfemReal z = 0.0);

    CfemInt addCell(CellType type,
                    const std::vector<CfemInt>& v_ids);

    void buildInverseConnectivities();
    void buildFacets();

    void buildMesh();

    CfemInt getDimension() const { return m_meshDim; }
    CfemInt getOrder() const { return m_meshOrder; }
    size_t getNumVertices() const { return m_vertices.size(); }
    size_t getNumCells()    const { return m_cells.size(); }
    size_t getNumFacets()   const { return m_facets.size(); }

    const Vertex& getVertex(CfemInt id) const { return m_vertices.at(id); }
    Vertex& getVertex(CfemInt vid) {  return m_vertices.at(vid); }
    const Cell&   getCell(CfemInt id)   const { return m_cells.at(id); }
    const Facet&  getFacet(CfemInt id)  const { return m_facets.at(id); }

    void setOrder(int newOrder);

    CfemInt createEntity(const std::string& name, CfemInt dim);

    MeshEntity& createInternalEntity(CfemInt id,
                                     const std::string& name,
                                     CfemInt dim);


    MeshEntity* getEntity(const std::string& name);
    MeshEntity* getEntity(CfemInt id);

    std::string getMeshFormat() const;

    void upgradeToQuadratic();
    void downgradeToLinear();
    void updateEntities();

    // HELPERS
    std::string cellTypeToString(CellType type) const;
    CfemReal computeCellMeasure(CellType type, const std::vector<CfemInt>& v_ids) const; 
    CfemReal compute2NodesLineLength(CfemInt p_idx0, CfemInt p_idx1) const;
    CfemReal computeSigned3NodesTriangleArea(CfemInt p_idx0, CfemInt p_idx1, CfemInt p_idx2, Vector3D &pReferenceNormale) const;
    CfemReal computeSigned4NodesTetraVolume(CfemInt p_idx0, CfemInt p_idx1, CfemInt p_idx2, CfemInt p_idx3) const;
    CfemInt getDimFromType(CellType type) const;
    CfemInt getOrderFromType(CellType type) const;

    CfemReal computeCellQuality(CfemInt cellId) const;

    void clear() {
        m_vertices.clear();
        m_cells.clear();
        m_facets.clear();
        m_entities.clear();
        m_entity_ids.clear();
        // Reset any other internal state if necessary
    }
    


private:
    CfemInt m_meshDim;
    CfemInt m_meshOrder;

    std::vector<Vertex> m_vertices;
    std::vector<Cell>   m_cells;
    std::vector<Facet>  m_facets;

    std::vector<std::unique_ptr<MeshEntity>> m_entities;

    std::unordered_set<CfemInt> m_entity_ids;

    void checkInvariant() const{
        CFEM_INVARIANT(m_meshDim >= 1 && m_meshDim <= 3);
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
    static void build(Mesh& mesh, int nx, const Vector3D& start, const Vector3D& end) {
        mesh.clear();
        for (int i = 0; i <= nx; ++i) {
            double t = static_cast<double>(i) / nx;
            Vector3D p = start * (1.0 - t) + end * t;
            mesh.addVertex(p.x, p.y, p.z);
        }

        for (int i = 0; i < nx; ++i) {
            mesh.addCell(CellType::LINE2, { (CfemInt)i, (CfemInt)i + 1 });
        }

        mesh.buildMesh();
        if (mesh.getOrder() == 2) mesh.upgradeToQuadratic(); // Becomes LINE3
    }
};

class QuadMeshBuilder {
public:
    enum class SubdivisionType { 
        None, FirstDiagonal, SecondDiagonal, CrissCross 
    };

    struct BoundaryNames {
        std::string bottom = "Bottom";
        std::string top    = "Top";
        std::string left   = "Left";
        std::string right  = "Right";
    };

    void build(Mesh& mesh, int nx, int ny, 
               const std::array<Vector3D, 4>& corners, 
               SubdivisionType subDiv = SubdivisionType::None) {
        build(mesh, nx, ny, corners, subDiv, BoundaryNames{});
    }

    void build(Mesh& mesh, int nx, int ny, 
               const std::array<Vector3D, 4>& corners, 
               SubdivisionType subDiv,
               const BoundaryNames& names);

private:
    void applyBoundaryTags(Mesh& mesh, int nx, int ny, 
                           const std::vector<CfemInt>& grid, 
                           const BoundaryNames& names);
};

class CubicMeshBuilder {
public:
    enum class SubdivisionType {
        None,       // Keep as HEX8
        SixTet      // Split each Hex into 6 Tetrahedra (standard Kuhn triangulation)
    };

    void build(Mesh& mesh, int nx, int ny, int nz, 
               const std::array<Vector3D, 8>& corners,
               SubdivisionType subDiv = SubdivisionType::None);

};

class DiskMeshBuilder {
public:
    enum class ElementType { 
        Quadrilateral, 
        Triangle 
    };

    void build(Mesh& mesh,
                      CfemInt order,
                      CfemReal target_mesh_size,
                      const Vector3D& center,
                      double radius,
                      ElementType type = ElementType::Triangle);

    static void projectBoundaryMidNodes(Mesh& mesh, const Vector3D& center, double radius);

    void relaxMesh(Mesh& mesh, int iterations = 5);
};

class SphereMeshBuilder {
public:
    enum class ElementType { Hexahedron, Tetrahedron };

    static void build(Mesh& mesh,
                      CfemInt order,
                      CfemReal h,
                      const Vector3D& center,
                      double radius,
                      ElementType type = ElementType::Tetrahedron);

private:
    static void projectBoundaryNodes(Mesh& mesh,
                                     const Vector3D& center,
                                     double radius);

    static void relaxMesh(Mesh& mesh, int iterations = 5);
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

#include "Mesh.tpp"

#endif


