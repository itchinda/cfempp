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

#include <algorithm>
#include <array>
#include <cmath>
#include <gtest/gtest.h>
#include <map>
#include <set>
#include <vector>

#include "cfemutils.h"
#include "Mesh.h"

using namespace cfem;

// -----------------------------------------------------------------------------
// EXTENDED TEST SUITE FOR MESH GENERATION
// -----------------------------------------------------------------------------

class MeshSystemTest : public ::testing::Test
{
protected:
    std::shared_ptr<Mesh> mesh;

    // Helper: Calculate the distance from a point to a line segment
    double distToSegment(Vector3D p, Vector3D a, Vector3D b)
    {
        Vector3D ab = b - a;
        Vector3D ap = p - a;
        double t = (ap.x * ab.x + ap.y * ab.y + ap.z * ab.z) / (ab.norm() * ab.norm());
        t = std::max(0.0, std::min(1.0, t));
        Vector3D proj = a + ab * t;
        return (p - proj).norm();
    }
};

/**
 * @brief Validates self-healing for an inverted TRI3 (Linear Triangle)
 */
TEST_F(MeshSystemTest, Tri3ConnectivityCheck)
{
    mesh = Mesh::create();
    mesh->setTopologicalDimension(2);
    mesh->setOrder(MeshOrder::Linear);
    mesh->setReferenceNormal({0, 0, 1});

    mesh->addVertex(0.0, 0.0, 0.0); // 0
    mesh->addVertex(1.0, 0.0, 0.0); // 1
    mesh->addVertex(0.0, 1.0, 0.0); // 2

    // Inverted: Swap 0 and 1
    std::vector<CfemInt> inverted_tri3 = {1, 0, 2};
    mesh->addCell(CellType::TRI3, inverted_tri3);

    ASSERT_NO_THROW(mesh->buildMesh());

#ifndef NDEBUG
    // Should be healed to {0, 1, 2}
    ASSERT_EQ(mesh->getCellConnectivity(0)[0], 0);
    ASSERT_EQ(mesh->getCellConnectivity(0)[1], 1);
    ASSERT_EQ(mesh->getCellConnectivity(0)[2], 2);
    ASSERT_GT(mesh->computeCellMeasure(CellType::TRI3, mesh->getCellConnectivity(0)), 0.0);
#endif
}

/**
 * @brief Validates self-healing for an inverted TRI6 (Quadratic Triangle)
 */
TEST_F(MeshSystemTest, Tri6ConnectivityCheck)
{
    mesh = Mesh::create();
    mesh->setTopologicalDimension(2);
    mesh->setOrder(MeshOrder::Quadratic);
    mesh->setReferenceNormal({0, 0, 1});

    // Corners
    mesh->addVertex(0.0, 0.0, 0.0); // 0
    mesh->addVertex(1.0, 0.0, 0.0); // 1
    mesh->addVertex(0.0, 1.0, 0.0); // 2
    // Mid-nodes (VTK Order: 3 on 0-1, 4 on 1-2, 5 on 2-0)
    mesh->addVertex(0.5, 0.0, 0.0); // 3
    mesh->addVertex(0.5, 0.5, 0.0); // 4
    mesh->addVertex(0.0, 0.5, 0.0); // 5

    // Inverted: Swap corners 0,1 AND mid-nodes 4,5
    std::vector<CfemInt> inverted_tri6 = {1, 0, 2, 3, 5, 4};
    mesh->addCell(CellType::TRI6, inverted_tri6);

    ASSERT_NO_THROW(mesh->buildMesh());

#ifndef NDEBUG
    for (size_t i = 0; i < 6; ++i)
        ASSERT_EQ(mesh->getCellConnectivity(0)[i], i);
    ASSERT_GT(mesh->computeCellMeasure(CellType::TRI6, mesh->getCellConnectivity(0)), 0.0);
#endif
}

/**
 * @brief Validates self-healing for an inverted QUAD4 (Linear Quad)
 */
TEST_F(MeshSystemTest, Quad4ConnectivityCheck)
{
    mesh = Mesh::create();
    mesh->setTopologicalDimension(2);
    mesh->setOrder(MeshOrder::Linear);
    mesh->setReferenceNormal({0, 0, 1});

    mesh->addVertex(0.0, 0.0, 0.0); // 0
    mesh->addVertex(1.0, 0.0, 0.0); // 1
    mesh->addVertex(1.0, 1.0, 0.0); // 2
    mesh->addVertex(0.0, 1.0, 0.0); // 3

    // Inverted: Swap 1 and 3
    std::vector<CfemInt> inverted_quad4 = {0, 3, 2, 1};
    mesh->addCell(CellType::QUAD4, inverted_quad4);

    ASSERT_NO_THROW(mesh->buildMesh());

#ifndef NDEBUG
    for (size_t i = 0; i < 4; ++i)
        ASSERT_EQ(mesh->getCellConnectivity(0)[i], i);
    ASSERT_GT(mesh->computeCellMeasure(CellType::QUAD4, mesh->getCellConnectivity(0)), 0.0);
#endif
}

/**
 * @brief Validates self-healing for an inverted QUAD9 (Bi-quadratic Quad)
 */
TEST_F(MeshSystemTest, Quad9ConnectivityCheck)
{
    mesh = Mesh::create();
    mesh->setTopologicalDimension(2);
    mesh->setOrder(MeshOrder::Quadratic);
    mesh->setReferenceNormal({0, 0, 1});

    mesh->addVertex(0.0, 0.0, 0.0); // 0
    mesh->addVertex(1.0, 0.0, 0.0); // 1
    mesh->addVertex(1.0, 1.0, 0.0); // 2
    mesh->addVertex(0.0, 1.0, 0.0); // 3
    // Edges
    mesh->addVertex(0.5, 0.0, 0.0); // 4 (0-1)
    mesh->addVertex(1.0, 0.5, 0.0); // 5 (1-2)
    mesh->addVertex(0.5, 1.0, 0.0); // 6 (2-3)
    mesh->addVertex(0.0, 0.5, 0.0); // 7 (3-0)
    // Center
    mesh->addVertex(0.5, 0.5, 0.0); // 8

    // Inverted: Swap corners (1,3) and edge-nodes (4,7) and (5,6)
    std::vector<CfemInt> inverted_quad9 = {0, 3, 2, 1, 7, 6, 5, 4, 8};
    mesh->addCell(CellType::QUAD9, inverted_quad9);

    ASSERT_NO_THROW(mesh->buildMesh());

#ifndef NDEBUG
    for (size_t i = 0; i < 9; ++i)
        ASSERT_EQ(mesh->getCellConnectivity(0)[i], i);
#endif
}

/**
 * @brief
 *
 */
TEST_F(MeshSystemTest, Hex27ConnectivityCheck)
{
    mesh = Mesh::create();

    mesh->setTopologicalDimension(3);
    mesh->setOrder(MeshOrder::Quadratic);
    // --- 1. Create the physical vertices in perfect VTK HEX27 order ---
    // Corners (0-7)
    mesh->addVertex(0.0, 0.0, 0.0); // 0 (Bottom)
    mesh->addVertex(1.0, 0.0, 0.0); // 1
    mesh->addVertex(1.0, 1.0, 0.0); // 2
    mesh->addVertex(0.0, 1.0, 0.0); // 3
    mesh->addVertex(0.0, 0.0, 1.0); // 4 (Top)
    mesh->addVertex(1.0, 0.0, 1.0); // 5
    mesh->addVertex(1.0, 1.0, 1.0); // 6
    mesh->addVertex(0.0, 1.0, 1.0); // 7
    // Horizontal Edges (8-15)
    mesh->addVertex(0.5, 0.0, 0.0); // 8 (Bottom 0-1)
    mesh->addVertex(1.0, 0.5, 0.0); // 9 (Bottom 1-2)
    mesh->addVertex(0.5, 1.0, 0.0); // 10 (Bottom 2-3)
    mesh->addVertex(0.0, 0.5, 0.0); // 11 (Bottom 3-0)
    mesh->addVertex(0.5, 0.0, 1.0); // 12 (Top 4-5)
    mesh->addVertex(1.0, 0.5, 1.0); // 13 (Top 5-6)
    mesh->addVertex(0.5, 1.0, 1.0); // 14 (Top 6-7)
    mesh->addVertex(0.0, 0.5, 1.0); // 15 (Top 7-4)
    // Vertical Edges (16-19)
    mesh->addVertex(0.0, 0.0, 0.5); // 16 (0-4)
    mesh->addVertex(1.0, 0.0, 0.5); // 17 (1-5)
    mesh->addVertex(1.0, 1.0, 0.5); // 18 (2-6)
    mesh->addVertex(0.0, 1.0, 0.5); // 19 (3-7)
    // Face Centers (20-25) matching your map: Left, Right, Front, Back, Bottom, Top
    mesh->addVertex(0.0, 0.5, 0.5); // 20 Left (x=0)
    mesh->addVertex(1.0, 0.5, 0.5); // 21 Right (x=1)
    mesh->addVertex(0.5, 0.0, 0.5); // 22 Front (y=0)
    mesh->addVertex(0.5, 1.0, 0.5); // 23 Back (y=1)
    mesh->addVertex(0.5, 0.5, 0.0); // 24 Bottom (z=0)
    mesh->addVertex(0.5, 0.5, 1.0); // 25 Top (z=1)
    // Body Center (26)
    mesh->addVertex(0.5, 0.5, 0.5); // 26 Center
    // --- 2. Add the INVERTED Cell ---
    // We feed the node IDs intentionally "upside down" relative to the Z-axis.
    std::vector<CfemInt> inverted_hex27 = {
        4, 5, 6, 7,     // Corners (Feeding Top nodes as Bottom)
        0, 1, 2, 3,     // Corners (Feeding Bottom nodes as Top)
        12, 13, 14, 15, // Horiz Edges (Feeding Top edges as Bottom)
        8, 9, 10, 11,   // Horiz Edges (Feeding Bottom edges as Top)
        16, 17, 18, 19, // Vertical Edges (Mid-height, stay the same)
        20, 21, 22, 23, // Side Face Centers (Mid-height, stay the same)
        25,             // Bottom Face Center (Intentionally giving it the Top center!)
        24,             // Top Face Center (Intentionally giving it the Bottom center!)
        26              // Body Center
    };
    mesh->addCell(CellType::HEX27, inverted_hex27);
    ASSERT_NO_THROW(mesh->buildMesh());

#ifndef NDEBUG
    for (size_t i = 0; i < 27; ++i) // Validate the internal connectivity fix in debug mode
    {
        ASSERT_EQ(mesh->getCellConnectivity(0)[i], i);
    }

    ASSERT_GT(mesh->computeCellMeasure(CellType::HEX27, mesh->getCellConnectivity(0)), 0.0);
#endif
}

/**
 * @test Verify that quadratic mid-nodes are placed EXACTLY on the linear edges.
 */
TEST_F(MeshSystemTest, QuadraticNodePlacementAccuracy)
{
    Mesh::QuadCorners corners({Vector3D{0, 0, 0}, Vector3D{10, 0, 0}, Vector3D{10, 10, 0}, Vector3D{0, 10, 0}});
    mesh = Mesh::quad(MeshOrder::Quadratic,
                      Mesh::QuadSubdivisions(1, 1), corners,
                      Mesh::QuadElementType::Quadrilateral, {});

    // For a QUAD9, we expect nodes 4, 5, 6, 7 to be mid-edges.
    const auto &cell_0_vids = mesh->getCellConnectivity(0);
    ASSERT_EQ(cell_0_vids.size(), 9);

    // Mid-node between v0(0,0) and v1(10,0) should be (5,0)
    Vector3D mid01 = mesh->getVertex(cell_0_vids[4]).coords;
    EXPECT_NEAR(mid01.x, 5.0, 1e-9);
    EXPECT_NEAR(mid01.y, 0.0, 1e-9);

    // Mid-node between v1(10,0) and v2(10,10) should be (10,5)
    Vector3D mid12 = mesh->getVertex(cell_0_vids[5]).coords;
    EXPECT_NEAR(mid12.x, 10.0, 1e-9);
    EXPECT_NEAR(mid12.y, 5.0, 1e-9);
}

/**
 * @test Ensures that the 'getMidpoint' map prevents duplicate vertices on shared edges.
 */
TEST_F(MeshSystemTest, QuadraticNodeSharingBetweenElements)
{
    Mesh::QuadCorners corners({Vector3D{0, 0, 0}, Vector3D{2, 0, 0}, Vector3D{2, 1, 0}, Vector3D{0, 1, 0}});
    mesh = Mesh::quad(MeshOrder::Quadratic,
                      Mesh::QuadSubdivisions(2, 1), corners,
                      Mesh::QuadElementType::Quadrilateral, {});

    // Linear nodes: (2+1)*(1+1) = 6
    // Edges: 7 (3 horizontal, 4 vertical).
    // Quad Centers: 2
    // Total expected: 6 + 7 + 2 = 15
    EXPECT_EQ(mesh->getNumVertices(), 15);
}

/**
 * @test Validates the 'CrissCross' subdivision logic (1 Quad -> 4 Tris).
 */
TEST_F(MeshSystemTest, CrissCrossTopologyTest)
{
    Mesh::QuadCorners corners({Vector3D{0, 0, 0}, Vector3D{1, 0, 0}, Vector3D{1, 1, 0}, Vector3D{0, 1, 0}});
    mesh = Mesh::quad(MeshOrder::Linear, Mesh::QuadSubdivisions(1, 1), corners, Mesh::QuadElementType::TriangleCrissCross);

    EXPECT_EQ(mesh->getNumCells(), 4);
    EXPECT_EQ(mesh->getNumVertices(), 5); // 4 corners + 1 center
    for (CfemInt i = 0; i < 4; ++i)
    {
        EXPECT_EQ(mesh->getCell(i).type, CellType::TRI3);
    }
}

/**
 * @test 3D Mapping: Checks if Trilinear Interpolation inside Mesh is correct.
 */
TEST_F(MeshSystemTest, CubicTrilinearInterpolation)
{
    Mesh::BlockCorners corners({Vector3D{0, 0, 0}, Vector3D{2, 0, 0}, Vector3D{2, 2, 0}, Vector3D{0, 2, 0},
                                Vector3D{0, 0, 2}, Vector3D{2, 0, 2}, Vector3D{2, 2, 2}, Vector3D{0, 2, 2}});
    mesh = Mesh::block(MeshOrder::Linear, Mesh::BlockSubdivisions(2, 2, 2), corners, Mesh::BlockElementType::Hexahedron);

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
TEST_F(MeshSystemTest, BoundaryTaggingCompleteness)
{
    Mesh::QuadBoundaryNames names = {"xmin", "xmax", "ymin", "ymax"};
    Mesh::QuadCorners corners({Vector3D{0, 0, 0}, Vector3D{1, 0, 0}, Vector3D{1, 1, 0}, Vector3D{0, 1, 0}});
    mesh = Mesh::quad(MeshOrder::Linear, Mesh::QuadSubdivisions(10, 10), corners, Mesh::QuadElementType::Quadrilateral, names);

    auto *ymax = mesh->getEntity("ymax");
    auto *ymin = mesh->getEntity("ymin");

    ASSERT_NE(ymax, nullptr);
    ASSERT_NE(ymin, nullptr);

    // 10x10 grid has 10 segments on each boundary
    EXPECT_EQ(ymax->getFacetIds().size(), 10);
    EXPECT_EQ(ymin->getFacetIds().size(), 10);
}

/**
 * @test Stress test for Relaxation (Laplacian).
 * Ensures that a highly distorted interior node moves toward the average of its neighbors.
 */
TEST_F(MeshSystemTest, LaplacianSmoothingConvergence)
{
    // Manually create a single QUAD4 with a distorted center node
    mesh = Mesh::create();
    mesh->setTopologicalDimension(2);
    mesh->setOrder(MeshOrder::Linear);
    CfemInt v0 = mesh->addVertex(0, 0, 0);
    CfemInt v1 = mesh->addVertex(2, 0, 0);
    CfemInt v2 = mesh->addVertex(2, 2, 0);
    CfemInt v3 = mesh->addVertex(0, 2, 0);
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
TEST_F(MeshSystemTest, Hex8ToHex27Upgrade)
{
    Mesh::BlockCorners corners({Vector3D{0, 0, 0}, Vector3D{1, 0, 0}, Vector3D{1, 1, 0}, Vector3D{0, 1, 0},
                                Vector3D{0, 0, 1}, Vector3D{1, 0, 1}, Vector3D{1, 1, 1}, Vector3D{0, 1, 1}});
    mesh = Mesh::block(MeshOrder::Linear, Mesh::BlockSubdivisions(1, 1, 1), corners, Mesh::BlockElementType::Hexahedron);

    mesh->upgradeToQuadratic();

    ASSERT_EQ(mesh->getNumCells(), 1);
    EXPECT_EQ(mesh->getCell(0).type, CellType::HEX27);
    EXPECT_EQ(mesh->getNumVertices(), 27);

    // Check if the body center node (v26) is actually at (0.5, 0.5, 0.5)
    Vector3D bodyCenter = mesh->getVertex(mesh->getCellConnectivity(0)[26]).coords;
    EXPECT_NEAR(bodyCenter.x, 0.5, 1e-9);
    EXPECT_NEAR(bodyCenter.y, 0.5, 1e-9);
    EXPECT_NEAR(bodyCenter.z, 0.5, 1e-9);
}

// -----------------------------------------------------------------------------

#include <memory>
TEST(MeshTest, InverseConnectivityCSRValidation)
{
    auto mesh = Mesh::create();

    mesh->setTopologicalDimension(2);
    mesh->setOrder(MeshOrder::Linear);

    mesh->addVertex(0.0, 0.0, 0.0); // ID 0
    mesh->addVertex(1.0, 0.0, 0.0); // ID 1
    mesh->addVertex(1.0, 1.0, 0.0); // ID 2
    mesh->addVertex(0.0, 1.0, 0.0); // ID 3
    mesh->addVertex(0.5, 0.5, 0.0); // ID 4

    mesh->addCell(CellType::TRI3, {0, 1, 4}); // Cell 0
    mesh->addCell(CellType::TRI3, {1, 2, 4}); // Cell 1
    mesh->addCell(CellType::TRI3, {2, 3, 4}); // Cell 2
    mesh->addCell(CellType::TRI3, {3, 0, 4}); // Cell 3

    mesh->buildMesh();

    // 4. Validation: Check Vertex 1 (should be connected to Cell 0 and Cell 1)
    const auto &v1 = mesh->getVertices()[1];

    ASSERT_EQ(v1.num_connected_cells, 2)
        << "Vertex 1 should be connected to exactly 2 cells.";
    ASSERT_NE(v1.connected_cells_ptr, nullptr);

    // Collect IDs from the CSR pointer
    std::vector<CfemInt> connected_ids;
    for (CfemInt i = 0; i < v1.num_connected_cells; ++i)
    {
        connected_ids.push_back(v1.connected_cells_ptr[i]);
    }

    // Sort to ensure order doesn't fail the test
    std::sort(connected_ids.begin(), connected_ids.end());

    EXPECT_EQ(connected_ids[0], 0);
    EXPECT_EQ(connected_ids[1], 1);

    const auto &v4 = mesh->getVertices()[4];
    EXPECT_EQ(v4.num_connected_cells, 4);
}

TEST(MeshTopologyTest, BuildFacetsConnectivity)
{
    auto mesh = Mesh::create();

    mesh->setTopologicalDimension(2);
    mesh->setOrder(MeshOrder::Linear);

    mesh->addVertex(0.0, 0.0, 0.0); // ID 0
    mesh->addVertex(1.0, 0.0, 0.0); // ID 1
    mesh->addVertex(1.0, 1.0, 0.0); // ID 2
    mesh->addVertex(0.0, 1.0, 0.0); // ID 3
    mesh->addVertex(0.5, 0.5, 0.0); // ID 4

    mesh->addCell(CellType::TRI3, {0, 1, 4}); // Cell 0
    mesh->addCell(CellType::TRI3, {1, 2, 4}); // Cell 1
    mesh->addCell(CellType::TRI3, {2, 3, 4}); // Cell 2
    mesh->addCell(CellType::TRI3, {3, 0, 4}); // Cell 3

    mesh->buildMesh();

    const auto &facets = mesh->getFacets();

    // 3. Verification
    EXPECT_EQ(facets.size(), 8);

    CfemInt internalCount = 0;
    CfemInt boundaryCount = 0;
    bool sharedEdgeFound = false;

    for (const auto &f : facets)
    {
        auto f_vids = mesh->getFacetConnectivity(f.id);
        std::set<CfemInt> v_set(f_vids.begin(), f_vids.end());
        if (f.right_cell != -1)
        {
            internalCount++;

            // Check the shared edge {0, 2}
            // Use std::set or sort to check vertex equality regardless of order
            std::set<CfemInt> expected_set = {1, 4};

            if (v_set == expected_set)
            {
                sharedEdgeFound = true;
                // Verify left/right connectivity
                EXPECT_TRUE((f.left_cell == 0 && f.right_cell == 1) ||
                            (f.left_cell == 1 && f.right_cell == 0));
            }
        }
        else
        {
            boundaryCount++;

            std::set<CfemInt> expected_set = {0, 3};
            if (v_set == expected_set)
            {
                // Verify left/right connectivity
                EXPECT_TRUE((f.left_cell == 3 && f.right_cell == -1));
            }
        }
    }

    EXPECT_EQ(internalCount, 4) << "There should be exactly 4 shared edge.";
    EXPECT_EQ(boundaryCount, 4) << "There should be 4 boundary edges.";
    EXPECT_TRUE(sharedEdgeFound) << "The facet connecting vertex 1 and 4 was not found.";
}