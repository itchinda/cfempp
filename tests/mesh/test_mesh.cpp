#include <gtest/gtest.h>
#include <vector>
#include <map>
#include <set>
#include <cmath>
#include <algorithm>
#include <array>
#include "cfemutils.h"
#include "Mesh.h"

using namespace cfem;

// -----------------------------------------------------------------------------
// EXTENDED TEST SUITE FOR MESH GENERATION
// -----------------------------------------------------------------------------

class MeshSystemTest : public ::testing::Test {
protected:
    std::shared_ptr<Mesh> mesh;

    // Helper: Calculate the distance from a point to a line segment
    double distToSegment(Vector3D p, Vector3D a, Vector3D b) {
        Vector3D ab = b - a;
        Vector3D ap = p - a;
        double t = (ap.x * ab.x + ap.y * ab.y + ap.z * ab.z) / (ab.norm() * ab.norm());
        t = std::max(0.0, std::min(1.0, t));
        Vector3D proj = a + ab * t;
        return (p - proj).norm();
    }
};

/**
 * @test Verify that quadratic mid-nodes are placed EXACTLY on the linear edges.
 */
TEST_F(MeshSystemTest, QuadraticNodePlacementAccuracy) {
    QuadMeshBuilder::Corners corners{ Vector3D{0,0,0}, Vector3D{10,0,0}, Vector3D{10,10,0}, Vector3D{0,10,0} };
    mesh = QuadMeshBuilder::create( MeshOrder::Quadratic, 
                                    QuadMeshBuilder::Subdivisions{1,1}, corners, 
                                    QuadMeshBuilder::ElementType::None, {});

    // For a QUAD9, we expect nodes 4, 5, 6, 7 to be mid-edges.
    const auto& cell = mesh->getCell(0);
    ASSERT_EQ(cell.vertex_ids.size(), 9);

    // Mid-node between v0(0,0) and v1(10,0) should be (5,0)
    Vector3D mid01 = mesh->getVertex(cell.vertex_ids[4]).coords;
    EXPECT_NEAR(mid01.x, 5.0, 1e-9);
    EXPECT_NEAR(mid01.y, 0.0, 1e-9);

    // Mid-node between v1(10,0) and v2(10,10) should be (10,5)
    Vector3D mid12 = mesh->getVertex(cell.vertex_ids[5]).coords;
    EXPECT_NEAR(mid12.x, 10.0, 1e-9);
    EXPECT_NEAR(mid12.y, 5.0, 1e-9);
}

/**
 * @test Ensures that the 'getMidpoint' map prevents duplicate vertices on shared edges.
 */
TEST_F(MeshSystemTest, QuadraticNodeSharingBetweenElements) {
    QuadMeshBuilder::Corners corners{ Vector3D{0,0,0}, Vector3D{2,0,0}, Vector3D{2,1,0}, Vector3D{0,1,0} };
    mesh = QuadMeshBuilder::create( MeshOrder::Quadratic, 
                                    QuadMeshBuilder::Subdivisions{2,1}, corners, 
                                    QuadMeshBuilder::ElementType::None, {});

    // Linear nodes: (2+1)*(1+1) = 6
    // Edges: 7 (3 horizontal, 4 vertical). 
    // Quad Centers: 2
    // Total expected: 6 + 7 + 2 = 15
    EXPECT_EQ(mesh->getNumVertices(), 15);
}

/**
 * @test Validates the 'CrissCross' subdivision logic (1 Quad -> 4 Tris).
 */
TEST_F(MeshSystemTest, CrissCrossTopologyTest) {
    QuadMeshBuilder::Corners corners = {
        Vector3D{0,0,0}, Vector3D{1,0,0}, Vector3D{1,1,0}, Vector3D{0,1,0}
    };
    mesh = QuadMeshBuilder::create(MeshOrder::Linear, QuadMeshBuilder::Subdivisions{1,1}, corners, QuadMeshBuilder::ElementType::CrissCross);

    EXPECT_EQ(mesh->getNumCells(), 4);
    EXPECT_EQ(mesh->getNumVertices(), 5); // 4 corners + 1 center
    for(int i=0; i<4; ++i) {
        EXPECT_EQ(mesh->getCell(i).type, CellType::TRI3);
    }
}

/**
 * @test 3D Mapping: Checks if Trilinear Interpolation inside CubicMeshBuilder is correct.
 */
TEST_F(MeshSystemTest, CubicTrilinearInterpolation) {
    CubicMeshBuilder::Corners corners = {
        Vector3D{0,0,0}, Vector3D{2,0,0}, Vector3D{2,2,0}, Vector3D{0,2,0},
        Vector3D{0,0,2}, Vector3D{2,0,2}, Vector3D{2,2,2}, Vector3D{0,2,2}
    };
    mesh = CubicMeshBuilder::create(MeshOrder::Linear, CubicMeshBuilder::Subdivisions{2,2,2}, corners, CubicMeshBuilder::ElementType::None);

    // The center vertex of a 2x2x2 grid should be exactly at (1,1,1)
    // In a (nx+1) grid, index for (1,1,1) is k*(3*3) + j*3 + i => 1*9 + 1*3 + 1 = 13
    Vector3D centerV = mesh->getVertex(13).coords;
    EXPECT_NEAR(centerV.x, 1.0, 1e-9);
    EXPECT_NEAR(centerV.y, 1.0, 1e-9);
    EXPECT_NEAR(centerV.z, 1.0, 1e-9);
}

/**
 * @test Verification of the Boundary Tagging logic.
 */
TEST_F(MeshSystemTest, BoundaryTaggingCompleteness) {
    QuadMeshBuilder::BoundaryNames names = {"Bottom", "Right", "Top", "Left"};
    QuadMeshBuilder::Corners corners = {
        Vector3D{0,0,0}, Vector3D{1,0,0}, Vector3D{1,1,0}, Vector3D{0,1,0}
    };
    mesh = QuadMeshBuilder::create(MeshOrder::Linear, QuadMeshBuilder::Subdivisions{10,10}, corners, QuadMeshBuilder::ElementType::None, names);

    auto* top = mesh->getEntity("Top");
    auto* bottom = mesh->getEntity("Bottom");
    
    ASSERT_NE(top, nullptr);
    ASSERT_NE(bottom, nullptr);
    
    // 10x10 grid has 10 segments on each boundary
    EXPECT_EQ(top->getFacetIds().size(), 10);
    EXPECT_EQ(bottom->getFacetIds().size(), 10);
}

/**
 * @test Stress test for Relaxation (Laplacian).
 * Ensures that a highly distorted interior node moves toward the average of its neighbors.
 */
TEST_F(MeshSystemTest, LaplacianSmoothingConvergence) {
    // Manually create a single QUAD4 with a distorted center node
    mesh = Mesh::create();
    mesh->setDimension(2);
    mesh->setOrder(MeshOrder::Linear);
    CfemInt v0 = mesh->addVertex(0,0,0);
    CfemInt v1 = mesh->addVertex(2,0,0);
    CfemInt v2 = mesh->addVertex(2,2,0);
    CfemInt v3 = mesh->addVertex(0,2,0);
    CfemInt vCenter = mesh->addVertex(0.1, 0.1, 0); // Distorted: should be at (1,1)

    mesh->addCell(CellType::TRI3, {v0, v1, vCenter});
    mesh->addCell(CellType::TRI3, {v1, v2, vCenter});
    mesh->addCell(CellType::TRI3, {v2, v3, vCenter});
    mesh->addCell(CellType::TRI3, {v3, v0, vCenter});

    mesh->buildMesh();

    std::set<CfemInt> boundaries = {v0, v1, v2, v3};
    
    mesh->relaxMesh(50);

    Vector3D finalPos = mesh->getVertex(vCenter).coords;
    EXPECT_NEAR(finalPos.x, 1.0, 1e-4);
    EXPECT_NEAR(finalPos.y, 1.0, 1e-4);
}

/**
 * @test High-Order 3D: Upgrade Hex8 to Hex27
 */
TEST_F(MeshSystemTest, Hex8ToHex27Upgrade) {
    CubicMeshBuilder::Corners corners = {
        Vector3D{0,0,0}, Vector3D{1,0,0}, Vector3D{1,1,0}, Vector3D{0,1,0},
        Vector3D{0,0,1}, Vector3D{1,0,1}, Vector3D{1,1,1}, Vector3D{0,1,1}
    };
    mesh = CubicMeshBuilder::create(MeshOrder::Linear, CubicMeshBuilder::Subdivisions{1,1,1}, corners, CubicMeshBuilder::ElementType::None);
    
    mesh->upgradeToQuadratic();

    ASSERT_EQ(mesh->getNumCells(), 1);
    EXPECT_EQ(mesh->getCell(0).type, CellType::HEX27);
    EXPECT_EQ(mesh->getNumVertices(), 27);
    
    // Check if the body center node (v26) is actually at (0.5, 0.5, 0.5)
    Vector3D bodyCenter = mesh->getVertex(mesh->getCell(0).vertex_ids[26]).coords;
    EXPECT_NEAR(bodyCenter.x, 0.5, 1e-9);
    EXPECT_NEAR(bodyCenter.y, 0.5, 1e-9);
    EXPECT_NEAR(bodyCenter.z, 0.5, 1e-9);
}

// -----------------------------------------------------------------------------