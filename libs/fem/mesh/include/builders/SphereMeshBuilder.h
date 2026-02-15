#ifndef CFEM_SPHERE_MESH_BUILDER_H
#define CFEM_SPHERE_MESH_BUILDER_H

#include "cfemutils.h"
#include "Vector3D.h"
#include "Mesh.h"

namespace cfem
{
    using SphereElementType = Mesh::SphereElementType;

    /**
     * @ingroup MeshBuilders
     * @brief Structured mesh generator for a solid sphere.
     *
     * This builder generates a volumetric mesh of a sphere of given center and radius.
     * The mesh is constructed from a structured Cartesian grid mapped onto the sphere using a smooth radial
     * projection.
     *
     * The resulting mesh can consist of:
     * - Tetrahedra
     * - Hexahedra
     *
     * The mesh is suitable for:
     * - Elliptic problems
     * - Stokes/Navier–Stokes simulations
     * - CutFEM background grids
     *
     * Quadratic meshes are supported through automatic upgrade.
     */
    class SphereMeshBuilder
    {
    private:
        struct ConstructorToken { };

    public:
    
        SphereMeshBuilder(ConstructorToken) {}

        /**
         * @brief Creates a sphere mesh.
         *
         * @param order            Polynomial order of the mesh (Linear or Quadratic).
         * @param target_mesh_size Target approximate mesh size.
         * @param center           Center of the sphere.
         * @param radius           Sphere radius.
         * @param type             Element type (Tetrahedron or Hexahedron).
         *
         * @return Shared pointer to the generated mesh.
         *
         * The mesh is generated from a structured Cartesian grid
         * enclosing the sphere and projected onto the spherical geometry.
         *
         * If quadratic order is requested, mid-edge nodes are projected
         * onto the exact spherical boundary.
         */
        static std::shared_ptr<Mesh> create(MeshOrder order,
                                            CfemReal target_mesh_size,
                                            const cfem::Vector3D &center,
                                            CfemReal radius,
                                            SphereElementType type = SphereElementType::Tetrahedron);

    private:
        /**
         * @brief Internal build routine.
         *
         * @param mesh   Mesh to populate.
         * @param order  Mesh order.
         * @param target_mesh_size Approximate mesh size.
         * @param center Sphere center.
         * @param radius Sphere radius.
         * @param type   Element type.
         *
         * The algorithm:
         * 1. Builds a structured background grid.
         * 2. Keeps cells inside the sphere.
         * 3. Projects boundary vertices onto the exact sphere.
         * 4. Upgrades to quadratic if requested.
         */
        static void build(Mesh &mesh,
                          MeshOrder order,
                          CfemReal target_mesh_size,
                          const cfem::Vector3D &center,
                          CfemReal radius,
                          SphereElementType type = SphereElementType::Tetrahedron);
    };

}

#endif