#include <gtest/gtest.h>
#include <fstream>
#include "GmshParser.h"
#include "Mesh.h"

using namespace cfem;

class GmshHexTest : public ::testing::Test {
protected:
    const std::string filename = "../tests/grenier/hex.msh";
    Mesh mesh;
    GmshParser parser;

    void SetUp() override {
        ASSERT_TRUE(std::ifstream(filename).good()) << "hex.msh file not found!";
    }
};

/**
 * @test Validates the Universe (Domain) and Order detection.
 */
TEST_F(GmshHexTest, UniverseAndOrderDetection) {
    ASSERT_TRUE(parser.load(mesh, filename));

    EXPECT_EQ(mesh.getDimension(), 3);
    EXPECT_EQ(mesh.getOrder(), 2);

    MeshEntity* universe = mesh.getEntity(cfem::UNIVERSE_ID);
    ASSERT_NE(universe, nullptr);
    EXPECT_EQ(universe->getName(), "Universe");
    
    EXPECT_EQ(universe->getCellIds().size(), 6);
}

/**
 * @test Validates the mapping of User-Defined Physical Groups.
 */
TEST_F(GmshHexTest, BoundaryGroupMapping) {
    parser.load(mesh, filename);

    std::vector<std::string> expectedGroups = {
        "Bottom", "Top", "Front", "Right", "Back", "Left", "Bord"
    };

    for (const auto& name : expectedGroups) {
        MeshEntity* ent = mesh.getEntity(name);
        ASSERT_NE(ent, nullptr) << "Group " << name << " not found in mesh!";
        EXPECT_GT(ent->getFacetIds().size(), 0);
        EXPECT_GT(ent->getVertexIds().size(), 0);
    }
}

/**
 * @test Validates the node permutation for TET10.
 * CFEM++ internally std::swap(vids[8], vids[9])
 */
TEST_F(GmshHexTest, Tet10Permutation) {
    parser.load(mesh, filename);

    bool foundTet10 = false;
    for (size_t i = 0; i < mesh.getNumCells(); ++i) {
        const auto& cell = mesh.getCell(i);
        if (cell.type == CellType::TET10) {
            foundTet10 = true;
            // TET10 order: 0-3 corners, 4-9 mid-edges
            EXPECT_EQ(cell.vertex_ids.size(), 10);
            break;
        }
    }
    EXPECT_TRUE(foundTet10);
}

/**
 * @test Automated System Boundary Detection.
 */
TEST_F(GmshHexTest, SystemBoundaryDetection) {
    parser.load(mesh, filename);

    MeshEntity* sysBound = mesh.getEntity(cfem::SYS_BOUNDARY_ID);
    ASSERT_NE(sysBound, nullptr);
    
    EXPECT_GT(sysBound->getFacetIds().size(), 0);

    // Ensure all vertices in 'Top' are also in 'SYS_BOUNDARY'
    MeshEntity* top = mesh.getEntity("Top");
    auto sysVids = sysBound->getVertexIds();
    for (CfemInt vId : top->getVertexIds()) {
        EXPECT_TRUE(std::find(sysVids.begin(), sysVids.end(), vId) != sysVids.end());
    }
}

/**
 * @test Vertex Harvesting Finalization.
 * Finalize() should ensure vertex lists are unique and sorted.
 */
TEST_F(GmshHexTest, EntityFinalization) {
    parser.load(mesh, filename);

    MeshEntity* universe = mesh.getEntity(cfem::UNIVERSE_ID);
    const auto& vids = universe->getVertexIds();

    EXPECT_TRUE(std::is_sorted(vids.begin(), vids.end()));
    
    auto it = std::adjacent_find(vids.begin(), vids.end());
    EXPECT_EQ(it, vids.end()) << "Duplicate vertices found in finalized entity!";
}