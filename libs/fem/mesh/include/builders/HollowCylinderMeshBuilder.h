#ifndef CFEM_HOLLOW_CYLINDER_MESH_BUILDER_H
#define CFEM_HOLLOW_CYLINDER_MESH_BUILDER_H

#include "cfemutils.h"
#include "cfem_topology.h"
#include "Vector3D.h"
#include "Mesh.h"
#include <map>
#include <vector>
#include <cmath>
#include <string>

namespace cfem
{
    using CylinderElementType = Mesh::CylinderElementType; // Assumes Wedge or Tetrahedron
    
    // Assumes your Mesh class has a struct for hollow cylinder boundaries
    using HollowCylinderBoundaryNames = Mesh::HollowCylinderBoundaryNames;

    /**
     * @ingroup MeshBuilders
     * @class HollowCylinderMeshBuilder
     * @brief Generates a quasi-uniform 3D mesh of a hollow cylinder (pipe).
     *
     * This builder constructs a structured annular discretization in the
     * cross-section (XY-plane) using concentric rings, then extrudes it
     * along the local Z-axis to generate a full 3D volume mesh.
     *
     * The mesh is initially generated as a **linear mesh** (TET4 or PRISM6),
     * and optionally upgraded to quadratic order via Mesh::upgradeToQuadratic().
     *
     * ------------------------------------------------------------------------
     * Geometry Definition
     * ------------------------------------------------------------------------
     *
     * The hollow cylinder is defined by:
     *
     * - Inner radius  : r_in
     * - Outer radius  : r_out
     * - Height        : height
     * - Axis direction: axis (normalized internally)
     * - Center        : center
     *
     * The mesh is first generated at the origin, aligned with the global Z-axis,
     * then rotated to match the provided axis direction using Rodrigues'
     * rotation formula, and finally translated to the specified center.
     *
     * ------------------------------------------------------------------------
     * Element Types
     * ------------------------------------------------------------------------
     *
     * Supported 3D element types:
     *
     * - CylinderElementType::Wedge
     *     -> PRISM6 (linear)
     *
     * - CylinderElementType::Tetrahedron
     *     -> 3 x TET4 per extruded triangle
     *
     * If MeshOrder::Quadratic is requested:
     *
     * - The mesh is first built as linear
     * - Mesh::upgradeToQuadratic() is called
     * - Curved boundary nodes are projected exactly onto cylindrical surfaces
     *
     * The resulting elements become:
     *
     * - PRISM18 (quadratic wedge)
     * - TET10 (quadratic tetrahedron)
     *
     * No mixed linear/quadratic elements are permitted.
     *
     * ------------------------------------------------------------------------
     * Mesh Resolution
     * ------------------------------------------------------------------------
     *
     * The target_mesh_size parameter (h) controls:
     *
     * - Radial subdivision:
     *       num_rings \simeq (r_out - r_in) / h
     *
     * - Axial subdivision:
     *       nz \simeq height / h
     *
     * - Angular subdivision:
     *       n_\theta \simeq 2\pi r / h (minimum 8 points per ring)
     *
     * This produces a quasi-uniform mesh in radial, circumferential,
     * and axial directions.
     *
     * ------------------------------------------------------------------------
     * Boundary Tagging
     * ------------------------------------------------------------------------
     *
     * After mesh.buildMesh(), boundary facets are classified into
     * four surface entities:
     *
     * - Inner  : r = r_in
     * - Outer  : r = r_out
     * - Bottom : z = -height/2
     * - Top    : z = +height/2
     *
     * ------------------------------------------------------------------------
     * Quadratic Geometry Handling
     * ------------------------------------------------------------------------
     *
     * If quadratic order is requested:
     *
     * 1. Mesh::upgradeToQuadratic() inserts mid-edge nodes.
     * 2. Nodes belonging to the inner and outer cylindrical surfaces
     *    are projected radially to enforce exact geometry.
     *
     * This ensures optimal geometric accuracy and restores
     * high-order convergence rates.
     *
     * ------------------------------------------------------------------------
     * Usage Example
     * ------------------------------------------------------------------------
     *
     * auto mesh = HollowCylinderMeshBuilder::create(
     *     MeshOrder::Quadratic,
     *     0.05,
     *     {0.0, 0.0, 0.0},
     *     0.5,
     *     1.0,
     *     2.0,
     *     {0.0, 0.0, 1.0},
     *     CylinderElementType::Wedge
     * );
     *
     */
    class HollowCylinderMeshBuilder
    {
    private:
        struct ConstructorToken {};

    public:
        HollowCylinderMeshBuilder(ConstructorToken) {}

        /**
         * @brief Creates a new hollow cylinder mesh.
         * @param order Mesh order (Linear or Quadratic).
         * @param target_mesh_size The target element size (h).
         * @param center The 3D center of the hollow cylinder.
         * @param r_in The inner radius.
         * @param r_out The outer radius.
         * @param height The total height in the Z direction.
         * @param type The 3D element type to generate (Wedge or Tetrahedron).
         * @param names Custom names for the generated boundary entities.
         * @return std::shared_ptr<Mesh> The generated mesh.
         */
        static std::shared_ptr<Mesh> create(
            MeshOrder order,
            CfemReal target_mesh_size,
            const cfem::Vector3D &center,
            CfemReal r_in,
            CfemReal r_out,
            CfemReal height,
            const cfem::Vector3D &axis,
            CylinderElementType type = CylinderElementType::Wedge,
            const HollowCylinderBoundaryNames &names = HollowCylinderBoundaryNames{}
        );

    private:
        /**
         * @brief Internal mesh generation routine.
         *
         * Generates a hollow cylinder mesh at the origin (Z-axis aligned),
         * using concentric radial rings and axial extrusion.
         * The mesh is built as **linear** initially, then optionally upgraded
         * to quadratic if requested.
         *
         * Steps performed:
         * 1. Normalize axis and compute radial, circumferential, and axial subdivisions.
         * 2. Generate vertices for each ring and axial layer.
         * 3. Extrude triangles to form 3D elements (TET4 or PRISM6).
         * 4. Identify and tag boundary facets (inner, outer, top, bottom).
         * 5. Upgrade to quadratic order if requested and project boundary nodes.
         * 6. Apply rigid-body rotation and translation to match user-defined axis and center.
         *
         * @param mesh Reference to the Mesh object to populate.
         * @param order Requested mesh order (Linear or Quadratic).
         * @param h Target element size (used for radial, angular, and axial subdivisions).
         * @param center 3D center of the cylinder.
         * @param r_in Inner radius.
         * @param r_out Outer radius.
         * @param height Cylinder height along the axis.
         * @param axis Cylinder axis vector (arbitrary orientation).
         * @param type CylinderElementType (Wedge or Tetrahedron).
         * @param names Names for boundary entities (inner, outer, bottom, top).
         *
         * @note The mesh is always built linearly first. Quadratic upgrade,
         *       if requested, is applied after linear mesh construction.
         * @note Degenerate or reversed triangles are automatically corrected
         *       or skipped using NumericPolicy::isZero() with scale-aware tolerance.
         */
        void build(
            Mesh &mesh,
            MeshOrder order,
            CfemReal h,
            const cfem::Vector3D &center,
            CfemReal r_in,
            CfemReal r_out,
            CfemReal height,
            const cfem::Vector3D &axis,
            CylinderElementType type,
            const HollowCylinderBoundaryNames &names
        );

        static void projectQuadraticSurfaces(
            Mesh &mesh, 
            CfemReal r_in, 
            CfemReal r_out,
            const HollowCylinderBoundaryNames &names
        );
    };
}

#endif // CFEM_HOLLOW_CYLINDER_MESH_BUILDER_H