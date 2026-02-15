#ifndef CFEMPP_MESH_CONFORMITY_CHECKER
#define CFEMPP_MESH_CONFORMITY_CHECKER

#include "Mesh.h"
#include <string>
#include <vector>
#include <sstream>

namespace cfem {

struct ConformityReport {
    bool is_valid = true;
    CfemInt isolated_vertices = 0;
    CfemInt non_manifold_faces = 0;
    CfemInt inverted_elements = 0;
    CfemInt degenerate_elements = 0;
    CfemInt unreachable_cells = 0;

    std::string details;

    /**
     * @brief Generates a formatted string of the report.
     */
    std::string toString() const {
        std::ostringstream oss;
        oss << "\n========== MESH CONFORMITY REPORT ==========\n";
        oss << "Overall Status:      " << (is_valid ? "PASS" : "FAIL") << "\n";
        oss << "Isolated Vertices:   " << isolated_vertices << "\n";
        oss << "Non-Manifold Faces:  " << non_manifold_faces << "\n";
        oss << "Inverted Elements:   " << inverted_elements << "\n";
        oss << "Degenerate Elements: " << degenerate_elements << "\n";
        oss << "Unreachable Cells:   " << unreachable_cells << "\n";
        
        if (!details.empty()) {
            oss << "--------------------------------------------\n";
            oss << "Detailed Errors:\n" << details;
        }
        oss << "============================================\n";
        return oss.str();
    }
};

class MeshConformityChecker {
public:
    /**
     * @brief Performs a full topological and geometric audit of the mesh.
     * @param mesh The mesh to check.
     * @param check_geometry If true, computes Jacobian/measures for all cells.
     */
    static ConformityReport check(const Mesh& mesh, bool check_geometry = true);

private:
    // Generates a sorted node-ID key to identify unique faces
    static std::vector<CfemInt> getCanonicalFaceKey(const LocalFace& local_face);
};

} // namespace cfem

#endif