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

struct TagGroup {
    std::vector<int> gmsh_tags;
    int internal_id;
};

class Mesh {
public:
    explicit Mesh();
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

    void setOrder(int newOrder);
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
    
    
    void setTagGroups(const std::vector<std::vector<int>>& groups) { m_tagGroups = groups; }
    const std::vector<int>& getPhysicalTags(int groupId) const { return m_tagGroups.at(groupId); }

    void relaxMesh(int iterations, double lambda = 0.4);


private:
    CfemInt m_meshDim {0};
    CfemInt m_meshOrder {0};

    std::vector<Vertex> m_vertices;
    std::vector<Cell>   m_cells;
    std::vector<Facet>  m_facets;

    std::vector<std::unique_ptr<MeshEntity>> m_entities;

    std::unordered_set<CfemInt> m_entity_ids;

    std::vector<std::vector<int>> m_tagGroups; // Table de "interning"

    void checkInvariant() const{
        CFEM_INVARIANT(m_meshDim >= 0 && m_meshDim <= 3);
        CFEM_INVARIANT(m_meshOrder >= 0 && m_meshOrder <= 2);
        if (m_meshDim==0 || m_meshOrder == 0){
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

    void build(Mesh& mesh, int order, int nx, int ny, 
               const std::array<cfem::Vector3D, 4>& corners, 
               SubdivisionType subDiv = SubdivisionType::None) {
        build(mesh, order, nx, ny, corners, subDiv, BoundaryNames{});
    }

    void build(Mesh& mesh, int order, int nx, int ny, 
               const std::array<cfem::Vector3D, 4>& corners, 
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

    void build(Mesh& mesh, int order, int nx, int ny, int nz, 
               const std::array<cfem::Vector3D, 8>& corners,
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
                      const cfem::Vector3D& center,
                      double radius,
                      ElementType type = ElementType::Triangle);

    static void projectBoundaryMidNodes(Mesh& mesh, const cfem::Vector3D& center, double radius);

    void relaxMesh(Mesh& mesh, int iterations = 5);
};

class SphereMeshBuilder {
public:
    enum class ElementType { Hexahedron, Tetrahedron };

    static void build(Mesh& mesh,
                      CfemInt order,
                      CfemReal h,
                      const cfem::Vector3D& center,
                      double radius,
                      ElementType type = ElementType::Tetrahedron);

private:
    static void projectBoundaryNodes(Mesh& mesh,
                                     const cfem::Vector3D& center,
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

}  // namespace

#include "Mesh.tpp"

#endif


