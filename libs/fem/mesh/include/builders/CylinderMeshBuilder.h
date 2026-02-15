#ifndef CFEM_UNIFORM_CYLINDER_MESH_BUILDER_H
#define CFEM_UNIFORM_CYLINDER_MESH_BUILDER_H

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
    using CylinderBoundaryNames = Mesh::CylinderBoundaryNames;
    

    /**
     * @ingroup MeshBuilders
     * @class CylinderMeshBuilder
     * @brief Generates a high-quality quasi-uniform 3D mesh for a cylinder along an arbitrary axis.
     *
     * Constructs the mesh by generating a quasi-uniform 2D disk pattern 
     * using concentric rings and extruding it along a given axis.  
     * The mesh can be linear or quadratic and supports PRISM6 or TET4 elements.
     *
     * ### Features:
     * - Mesh along arbitrary axis (`axis` parameter)
     * - Automatic rigid-body rotation + translation to place the cylinder
     * - Quadratic upgrade with exact radial projection of mid-edge nodes
     * - High-quality element connectivity with orientation checks
     * - Boundary tagging:
     *     - Bottom, Top, Lateral surfaces
     *
     * ### Supported Element Types:
     * - `PRISM6` (extruded triangles)
     * - `TET4` (3 TET per extruded triangle)
     *
     * ### Mesh Generation Steps:
     * 1. Generate concentric rings of vertices in XY-plane around origin.
     * 2. Extrude rings along Z-axis to form prism/tetrahedron elements.
     * 3. Correct element orientation to ensure positive volume.
     * 4. Tag boundaries (bottom, top, lateral) for use in FEM conditions.
     * 5. Upgrade to quadratic elements if requested, projecting mid-nodes radially.
     * 6. Apply rotation to align cylinder with user-defined axis.
     * 7. Apply translation to move cylinder to user-defined center.
     */
    class CylinderMeshBuilder
    {
    private:
        struct ConstructorToken {};

    public:
        CylinderMeshBuilder(ConstructorToken) {}

        /**
         * @brief Creates a shared pointer to a cylinder mesh.
         * @param order Mesh order (linear/quadratic)
         * @param target_mesh_size Approximate edge size
         * @param center Cylinder center (translation after generation)
         * @param radius Cylinder radius
         * @param height Cylinder height
         * @param axis Axis vector along which the cylinder extends (normalized internally)
         * @param type Element type: Wedge or Tetrahedron
         * @param names Names for lateral, top, bottom boundaries
         * @return std::shared_ptr<Mesh> Constructed cylinder mesh
         */
        static std::shared_ptr<Mesh> create(
            MeshOrder order,
            CfemReal target_mesh_size,
            const cfem::Vector3D &center,
            CfemReal radius,
            CfemReal height,
            const cfem::Vector3D &axis = {0,0,1},
            CylinderElementType type = CylinderElementType::Wedge,
            const CylinderBoundaryNames &names = CylinderBoundaryNames{}
        );

    private:
        /**
         * @brief Builds the cylinder mesh at the origin, then applies axis rotation and translation
         * @param mesh Mesh object to populate
         * @param order Mesh order
         * @param h Target mesh size
         * @param center Desired center of the cylinder after transformation
         * @param radius Cylinder radius
         * @param height Cylinder height
         * @param axis Direction vector of cylinder axis (will be normalized)
         * @param type Element type
         * @param names Boundary names
         */
        void build(
            Mesh &mesh,
            MeshOrder order,
            CfemReal h,
            const cfem::Vector3D &center,
            CfemReal radius,
            CfemReal height,
            const cfem::Vector3D &axis,
            CylinderElementType type,
            const CylinderBoundaryNames &names
        );
    };
}

#endif // CFEM_UNIFORM_CYLINDER_MESH_BUILDER_H