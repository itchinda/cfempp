#ifndef ANNULUS_MESH_BUILDER_H
#define ANNULUS_MESH_BUILDER_H

#include "Mesh.h"
#include <vector>
#include <cmath>

namespace cfem {



using DiskElementType = Mesh::DiskElementType;
using AnnulusBoundaryNames = Mesh::AnnulusBoundaryNames;


/**
 * @ingroup MeshBuilders
 * @class AnnulusMeshBuilder
 * @brief Generates a quasi-uniform 2D annular mesh (disk with hole) in the XY-plane.
 *
 * Constructs the mesh by creating concentric radial rings and connecting
 * points in counter-clockwise order to form triangular elements.
 *
 * Features:
 * - Supports linear (P1) and quadratic (P2) elements.
 * - Automatic boundary tagging for inner and outer circular boundaries.
 * - Quadratic upgrade projects mid-side nodes onto exact circular boundaries.
 *
 * Usage example:
 * @code
 * auto mesh = cfem::AnnulusMeshBuilder::create(
 *     cfem::MeshOrder::Quadratic,  // Desired element order
 *     0.1,                         // Target element size
 *     {0,0,0},                     // Center of annulus
 *     1.0,                         // Inner radius
 *     2.0                          // Outer radius
 * );
 * @endcode
 */ 
class AnnulusMeshBuilder {

private:
    struct ConstructorToken {};
    AnnulusMeshBuilder (ConstructorToken) {};

public:
    /**
     * @brief Static factory method for creating an annulus mesh.
     *
     * Internally calls the private `build()` method.
     *
     * @param order Mesh order (Linear or Quadratic).
     * @param h Target element size (approximate).
     * @param center Center of the annulus in 3D space (Z is constant).
     * @param radius_in Inner radius.
     * @param radius_out Outer radius.
     * @param type Element type (Triangle by default).
     * @param names Optional names for boundary entities (inner, outer).
     * @return std::shared_ptr<Mesh> The generated 2D mesh.
     */
    static std::shared_ptr<Mesh> create(
        MeshOrder order,
        CfemReal h,
        const cfem::Vector3D &center,
        CfemReal radius_in,
        CfemReal radius_out,
        DiskElementType type = DiskElementType::Triangle,
        const AnnulusBoundaryNames &names = {}
    )
    {
        auto mesh = Mesh::create();
        AnnulusMeshBuilder builder(ConstructorToken{});
        builder.build(*mesh, order, h, center, radius_in, radius_out, type, names);
        return mesh;
    }

private:

    /**
     * @brief Internal mesh generation routine.
     *
     * Builds a concentric ring triangular mesh for an annulus at the given
     * center. Mesh is generated linearly first, then optionally upgraded
     * to quadratic order.
     *
     * Steps:
     * 1. Compute number of radial subdivisions based on target element size.
     * 2. Generate vertices for each concentric ring.
     * 3. Connect vertices in counter-clockwise order to form triangular cells.
     * 4. Build the linear mesh (Mesh::buildMesh()).
     * 5. Apply boundary tags for inner and outer circles.
     * 6. Upgrade to quadratic elements and project mid-side nodes to exact
     *    inner/outer radius if requested.
     *
     * @param mesh Reference to the Mesh object to populate.
     * @param order Requested mesh order (Linear or Quadratic).
     * @param h Target element size (used for radial and circumferential subdivisions).
     * @param center Center of the annulus in 3D space.
     * @param r_in Inner radius.
     * @param r_out Outer radius.
     * @param type Element type (Triangle).
     * @param names Names for boundary entities (inner, outer).
     *
     * @note The mesh is always built linearly first. Quadratic upgrade,
     *       if requested, is applied after linear mesh construction.
     * @note Counter-clockwise orientation is enforced for all triangles.
     */
    void build(
                Mesh &mesh,
                MeshOrder order,
                CfemReal h,
                const Vector3D &center,
                CfemReal r_in,
                CfemReal r_out,
                DiskElementType type = DiskElementType::Triangle,
                const AnnulusBoundaryNames &names = {}
            );
    
    /**
     * @brief Assigns boundary facets to inner and outer entities.
     *
     * Each boundary vertex is identified by its ring index.
     * Facets connecting only inner/outer ring vertices are tagged accordingly.
     */
    static void applyRadialBoundaryTags(
                                            Mesh &mesh, 
                                            const std::vector<std::vector<CfemInt>> &ring_nodes,
                                            const AnnulusBoundaryNames &names
                                        );

    /**
     * @brief Projects mid-side nodes of quadratic annulus elements
     *        onto exact circular boundaries.
     *
     * @param mesh The mesh to modify.
     * @param center Center of the annulus.
     * @param r_in Inner radius.
     * @param r_out Outer radius.
     * @param names Boundary entity names.
     */
    static void projectQuadraticAnnulus(
                                            Mesh &mesh, 
                                            const Vector3D &center, 
                                            CfemReal r_in, 
                                            CfemReal r_out,
                                            const AnnulusBoundaryNames &names
                                        );


};

} // namespace cfem
#endif