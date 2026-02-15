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
#include <string>
#include <unordered_set>
#include <stdexcept>
#include <memory>
#include <utility>
#include <map>
#include <set>
#include <span>
#include <optional>

#include "cfemutils.h"
#include "cfem_topology.h"
#include "MeshEntity.h"
#include "GmshParser.h"
#include "Vector3D.h"
#include "DynamicVector.h"
#include "GmshParser.h"


constexpr int MAX_CELL_TYPES {static_cast<int>(cfem::CellType::COUNT)};
namespace cfem{
    class QuadMeshBuilder;
    class LineMeshBuilder;
    class DiskMeshBuilder;
    class BlockMeshBuilder;
    class SphereMeshBuilder;
    class AnnulusMeshBuilder;


class FESpace;


struct TagGroup {
    std::vector<int> gmsh_tags;
    int internal_id;
};

/**
 * @brief Stack-allocated facet representation to avoid heap allocation.
 */
struct LocalFace {
    CfemInt v[9]; // Capacity for QUAD9 nodes
    int size = 0;
    CellType type = CellType::NODE;

    bool contains(CfemInt v_id) const {
        for (int i = 0; i < size; ++i) if (v[i] == v_id) return true;
        return false;
    }
    
    // Check if another LocalFace matches this one (topology-independent)
    bool matches(const LocalFace& other) const {
        if (size != other.size) return false;
        // Since nodes in facets are usually few, this nested loop is faster 
        // than sorting for 8 million elements.
        for (int i = 0; i < size; ++i) {
            if (!other.contains(v[i])) return false;
        }
        return true;
    }
};

template <typename T, typename MeshPtr>
class MeshIterator {
public:
    using iterator_category = std::forward_iterator_tag;
    using value_type = T;
    using pointer = T*;
    using reference = T&;

    MeshIterator(MeshPtr mesh, CfemInt id) : m_mesh(mesh), m_id(id) {}

    // Dereference returns the actual object via the Big Block accessors
    reference operator*() const { return m_mesh->template getSpecific<T>(m_id); }
    
    // Increment simply moves the index
    MeshIterator& operator++() { ++m_id; return *this; }
    bool operator!=(const MeshIterator& other) const { return m_id != other.m_id; }
    
    CfemInt id() const { return m_id; } // Access the raw index if needed for FE Assembler

private:
    MeshPtr m_mesh;
    CfemInt m_id;
};

class Mesh {
    friend class GmshParser;
    friend class QuadMeshBuilder;
    friend class BlockMeshBuilder;
    friend class DiskMeshBuilder;
    friend class AnnulusMeshBuilder;
    friend class SphereMeshBuilder;

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
    void buildAdjacency();
    void buildFacets();

    void buildMesh();
    void compactVertices();

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

    // Returns view of Cell nodes
    std::span<const CfemInt> getCellConnectivity(CfemInt c_id) const {
        const auto& c = m_cells[c_id];
        return { &m_allCellVertexIds[c.offset], c.node_count };
    }

    // Returns view of Facet nodes
    std::span<const CfemInt> getFacetConnectivity(CfemInt f_id) const {
        const auto& f = m_facets[f_id];
        return { &m_allFacetVertexIds[f.offset], f.node_count };
    }

    void getElementCoordinates(CfemInt elementId, cfem::DynamicVector<cfem::Vector3D>& localCoords) const;

    std::span<const CfemInt> getCellNeighbors(CfemInt c_id) const {
        size_t start = m_neighborOffsets[c_id];
        size_t count = m_neighborOffsets[c_id + 1] - start;
        
        return { &m_allCellNeighbors[start], count };
    }

    // Returns view of a Cell's facet IDs
    std::span<const CfemInt> getCellFacetIds(CfemInt c_id) const {
        size_t start = m_neighborOffsets[c_id];
        size_t count = m_neighborOffsets[c_id + 1] - start;
        return { &m_allCellFacetIds[start], count };
    }

    const std::vector<Vertex>& getVertices() const { return m_vertices; }
    const std::vector<Cell>&   getCells()   const { return m_cells; }
    const std::vector<Facet>&  getFacets()  const { return m_facets; }
    const std::vector<std::unique_ptr<MeshEntity>>&  getEntities()  const { return m_entities; }
    std::vector<std::unique_ptr<MeshEntity>>&  getEntities() { return m_entities; }

    void setOrder(MeshOrder newOrder);
    void setDimension(int newDim);
    void setReferenceNormal(const cfem::Vector3D& norm);
    void clearReferenceNormal();

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

    template<typename T>
    CfemReal computeCellMeasure(CellType type, const T& v_ids) const;
    //CfemReal computeCellMeasure(CellType type, const auto& v_ids) const; 
    CfemReal compute2NodesLineLength(CfemInt p_idx0, CfemInt p_idx1) const;
    CfemReal computeSigned3NodesTriangleArea(CfemInt p_idx0, CfemInt p_idx1, CfemInt p_idx2) const;
    CfemReal computeSigned4NodesTetraVolume(CfemInt p_idx0, CfemInt p_idx1, CfemInt p_idx2, CfemInt p_idx3) const;
    CfemInt getDimFromType(CellType type) const;
    MeshOrder getCellOrderFromType(CellType type) const;
    void getCellVerticesCoords(CfemInt cell_id, DynamicVector<Vector3D>& v_coords) const;
    const std::vector<CellType>& getCellTypesPresent() const { return m_uniqueTypes; }
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

    void relaxMesh(CfemInt iterations, CfemReal lambda = 0.4);

    void addObserver(std::shared_ptr<FESpace> space);

    void registerFieldName(const std::string& name){
        NamingPolicy::validate(name);
        m_globalFieldRegistry.insert(name);
    }

    bool nameAlreadyUsed(const std::string& name){
        return m_globalFieldRegistry.find(name) != m_globalFieldRegistry.end();
    }

    static std::shared_ptr<Mesh> loadFromGmsh(const std::string& filename, std::optional<Vector3D> pReferenceNormal = std::nullopt){
        GmshParser parser;
        return parser.load(filename, pReferenceNormal);
    }



    struct CellRange {
        const Mesh* mesh;
        auto begin() { return MeshIterator<const Cell, const Mesh*>(mesh, 0); }
        auto end()   { return MeshIterator<const Cell, const Mesh*>(mesh, mesh->getNumCells()); }
    };

    CellRange cells() const { return { this }; }

    // Range Proxy for Facets
    struct FacetRange {
        const Mesh* mesh;
        auto begin() { return MeshIterator<const Facet, const Mesh*>(mesh, 0); }
        auto end()   { return MeshIterator<const Facet, const Mesh*>(mesh, mesh->getNumFacets()); }
    };

    FacetRange facets() const { return { this }; }

    // Mesh builders
    struct QuadBounds;
    struct BlockBounds;
    struct QuadCorners;
    struct BlockCorners;

    enum class QuadElementType { 
        Quadrilateral, TriangleFirstDiag, TriangleSecondDiag, TriangleCrissCross
    };

    enum class BlockElementType { 
        Hexahedron, Tetrahedron, Edge
    };

    enum class DiskElementType
    {
        Quadrilateral,
        Triangle
    };

    enum class SphereElementType
    {
        Hexahedron,
        Tetrahedron
    };

    enum class CylinderElementType
    {
        Wedge,
        Tetrahedron
    };
    
    using TJunctionElementType = CylinderElementType;

    struct QuadSubdivisions{
        CfemInt nx {1};
        CfemInt ny {1};
        QuadSubdivisions(CfemInt p_nx=1, CfemInt p_ny = 1): nx(p_nx), ny(p_ny){
            if (nx <= 0 || ny <= 0){
                CFEM_ERROR("The number of subdivisions in both directions must be strictly positive.");
            }
        }
    };

    struct BlockSubdivisions{
        CfemInt nx {1};
        CfemInt ny {1};
        CfemInt nz {1};
        BlockSubdivisions(CfemInt p_nx=1, CfemInt p_ny = 1, CfemInt p_nz = 1): nx(p_nx), ny(p_ny), nz(p_nz){
            if (nx <= 0 || ny <= 0 || nz <= 0){
                CFEM_ERROR("The number of subdivisions in the three directions must be strictly positive.");
            }
        }
    };
    
    struct QuadBoundaryNames 
    {
        std::string xmin = "xmin";
        std::string xmax = "xmax";
        std::string ymin = "ymin";
        std::string ymax = "ymax";

    };

    struct BlockBoundaryNames 
    {
        std::string xmin = "xmin";
        std::string xmax = "xmax";
        std::string ymin = "ymin";
        std::string ymax = "ymax";
        std::string zmin = "zmin";
        std::string zmax = "zmax";

    };

    struct AnnulusBoundaryNames {
        std::string inner = "inner_circle";
        std::string outer = "outer_circle";
    };

    struct CylinderBoundaryNames {
        std::string lateral = "Lateral";
        std::string bottom  = "Bottom";
        std::string top     = "Top";
    };

    struct HollowCylinderBoundaryNames {
        std::string inner  = "Inner";
        std::string outer  = "Outer";
        std::string bottom = "Bottom";
        std::string top    = "Top";
    };

struct TJunctionBoundaryNames
    {
        std::string main_inlet  = "main_inlet";
        std::string main_outlet = "main_outlet";
        std::string branch_inlet = "branch_inlet";
        std::string wall = "wall";
    };

    static std::shared_ptr<Mesh> quad(MeshOrder order,
                                      const QuadSubdivisions& subDiv,
                                      const QuadBounds& bounds,
                                      QuadElementType elemType,
                                      const QuadBoundaryNames& names);
    static std::shared_ptr<Mesh> quad(MeshOrder order,
                                      const QuadSubdivisions& subDiv,
                                      const QuadBounds& bounds,
                                      QuadElementType elemType);

    static std::shared_ptr<Mesh> quad(MeshOrder order,
                                      const QuadSubdivisions& subDiv,
                                      const QuadCorners& corners,
                                      QuadElementType elemType,
                                      const QuadBoundaryNames& names);
    static std::shared_ptr<Mesh> quad(MeshOrder order,
                                      const QuadSubdivisions& subDiv,
                                      const QuadCorners& corners,
                                      QuadElementType elemType);

    static std::shared_ptr<Mesh> block(
                                        MeshOrder order,
                                        const BlockSubdivisions& subDiv, 
                                        const BlockBounds& bounds,
                                        BlockElementType elemType,
                                        const BlockBoundaryNames& names
                                    );

    static std::shared_ptr<Mesh> block(
                                        MeshOrder order,
                                        const BlockSubdivisions& subDiv, 
                                        const BlockBounds& bounds,
                                        BlockElementType elemType
                                    );

    static std::shared_ptr<Mesh> block(
                                        MeshOrder order,
                                        const BlockSubdivisions& subDiv, 
                                        const BlockCorners& corners,
                                        BlockElementType elemType,
                                        const BlockBoundaryNames& names
                                    );

    static std::shared_ptr<Mesh> block(
                                        MeshOrder order,
                                        const BlockSubdivisions& subDiv, 
                                        const BlockCorners& corners,
                                        BlockElementType elemType
                                    );

    static std::shared_ptr<Mesh> disk(
                                        MeshOrder order, 
                                        CfemReal mesh_size, 
                                        const Vector3D& center, 
                                        CfemReal radius,
                                        const cfem::Vector3D &normal = {0,0,1},
                                        DiskElementType elemtype = DiskElementType::Triangle
                                    );

    static std::shared_ptr<Mesh> disk(
                                        MeshOrder order, 
                                        CfemReal mesh_size, 
                                        const Vector3D& center, 
                                        CfemReal radius,
                                        DiskElementType elemtype
                                    );

    static std::shared_ptr<Mesh> annulus(
                                            MeshOrder order,
                                            CfemReal h,
                                            const cfem::Vector3D &center,
                                            CfemReal radius_in,
                                            CfemReal radius_out,
                                            DiskElementType type,
                                            const AnnulusBoundaryNames &names
                                        );

    static std::shared_ptr<Mesh> annulus(
                                            MeshOrder order,
                                            CfemReal h,
                                            const cfem::Vector3D &center,
                                            CfemReal radius_in,
                                            CfemReal radius_out,
                                            DiskElementType type = DiskElementType::Triangle
                                        );

    static std::shared_ptr<Mesh> sphere(
                                            MeshOrder order, 
                                            CfemReal mesh_size, 
                                            const Vector3D& center, 
                                            CfemReal radius,
                                            SphereElementType elemtype = SphereElementType::Hexahedron
                                        );

    static std::shared_ptr<Mesh> cylinder( 
                                            MeshOrder order,
                                            CfemReal target_mesh_size,
                                            const cfem::Vector3D &center,
                                            CfemReal radius,
                                            CfemReal height,
                                            const cfem::Vector3D &axis,
                                            CylinderElementType type,
                                            const CylinderBoundaryNames &names
                                        );

    static std::shared_ptr<Mesh> cylinder( 
                                            MeshOrder order,
                                            CfemReal target_mesh_size,
                                            const cfem::Vector3D &center,
                                            CfemReal radius,
                                            CfemReal height,
                                            const cfem::Vector3D &axis = {0,0,1},
                                            CylinderElementType type = CylinderElementType::Wedge
                                        );

    static std::shared_ptr<Mesh> hollow_cylinder(
                                                    MeshOrder order,
                                                    CfemReal target_mesh_size,
                                                    const cfem::Vector3D &center,
                                                    CfemReal r_in,
                                                    CfemReal r_out,
                                                    CfemReal height,
                                                    const cfem::Vector3D &axis,
                                                    CylinderElementType type,
                                                    const HollowCylinderBoundaryNames &names
                                                );
    
    static std::shared_ptr<Mesh> hollow_cylinder(
                                                    MeshOrder order,
                                                    CfemReal target_mesh_size,
                                                    const cfem::Vector3D &center,
                                                    CfemReal r_in,
                                                    CfemReal r_out,
                                                    CfemReal height,
                                                    const cfem::Vector3D &axis,
                                                    CylinderElementType type = CylinderElementType::Wedge
                                                );

    static void projectBoundaryMidNodes(Mesh &mesh, const cfem::Vector3D &center, CfemReal radius);

    int getLocalFacesCount(CellType type) const;
    LocalFace getLocalFaceData(CfemInt cell_id, int face_idx) const;
private:
    CfemInt m_meshDim {0};
    MeshOrder m_meshOrder {MeshOrder::None};

    std::vector<Vertex> m_vertices;

    std::vector<CfemInt> m_allCellVertexIds; // The "Big Block" for cells
    std::vector<size_t>  m_cellOffsets;      // Offsets for mixed types (Tet/Hex)
    std::vector<Cell>   m_cells;

    std::vector<CfemInt> m_allFacetVertexIds;
    std::vector<size_t>  m_facetOffsets;
    std::vector<Facet>  m_facets;

    std::vector<CfemInt> m_allCellNeighbors;  // Stores IDs of adjacent Cells
    std::vector<CfemInt> m_allCellFacetIds;   // Stores IDs of the Facets (edges/faces)
    std::vector<size_t>  m_neighborOffsets;   // Offsets into the neighbor blocks

    std::vector<std::unique_ptr<MeshEntity>> m_entities;
    std::unordered_set<CfemInt> m_entity_ids;

    std::vector<std::vector<int>> m_tagGroups; // Table de "interning"

    std::vector<std::weak_ptr<FESpace>> m_observers; ///> list of FESpace attached to the mesh

    std::bitset<MAX_CELL_TYPES> m_typeFlags;
    std::vector<CellType> m_uniqueTypes;

    std::vector<CfemInt> m_globalConnectivity;

    std::optional<cfem::Vector3D> m_global_reference_normal;
    std::map<CellType, size_t> m_invertedCellCounts;

    void checkInvariant() const{
        CFEM_INVARIANT(m_meshDim >= 0 && m_meshDim <= 3);
        if (m_meshDim==0 || m_meshOrder == MeshOrder::None){
            CFEM_INVARIANT(m_cells.empty() && m_facets.empty());
        }
    }

    void checkConsistency(CellType type, const std::vector<CfemInt>& v_ids) const;
    void fixElementOrientation(CellType type, std::vector<CfemInt> &vids);
    void validateEntityName(const std::string& name) const;

    CfemInt getMidpoint(CfemInt v1, CfemInt v2, auto& edge_map);
    CfemInt getQuadFaceCenter(CfemInt v1, CfemInt v2, CfemInt v3, CfemInt v4, auto &quadface_map);

    void fillLine3Nodes   (std::span<const CfemInt> v, auto& edge_map, std::vector<CfemInt>& quad_buffer_out);
    void fillTri6Nodes    (std::span<const CfemInt> v, auto& edge_map, std::vector<CfemInt>& quad_buffer_out);
    void fillQuad9Nodes   (std::span<const CfemInt> v, auto& edge_map, std::vector<CfemInt>& quad_buffer_out);
    void fillTet10Nodes   (std::span<const CfemInt> v, auto& edge_map, std::vector<CfemInt>& quad_buffer_out);
    void fillPrism18Nodes (std::span<const CfemInt> v, auto& edge_map, auto &quadface_to_midpoint_map, std::vector<CfemInt>& quad_buffer_out);
    void fillHex27Nodes   (std::span<const CfemInt> v, auto& edge_map, auto &quadface_to_midpoint_map, std::vector<CfemInt>& quad_buffer_out);

    CfemReal computeTriQuality(const std::span<const CfemInt>& v) const;
    CfemReal computeQuadQuality(const std::span<const CfemInt>& v) const;
    CfemReal computeTetQuality(const std::span<const CfemInt>& v) const;
    CfemReal computePrismQuality(const std::span<const CfemInt>& v) const;

    // int getLocalFacesCount(CellType type) const;
    // LocalFace getLocalFaceData(CfemInt cell_id, int face_idx) const;
    int findLocalFaceIndex(CfemInt cell_id, const LocalFace& target) const;

    std::set<std::string> m_globalFieldRegistry;


};



struct QualityReport {
    CfemReal minVolume = std::numeric_limits<CfemReal>::max();
    CfemReal maxVolume = -std::numeric_limits<CfemReal>::max();
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
            CfemReal vol = mesh.computeCellMeasure(cell.type, mesh.getCellConnectivity(cell.id) );
            
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


