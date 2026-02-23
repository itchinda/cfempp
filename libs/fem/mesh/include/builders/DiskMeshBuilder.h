#ifndef CFEM_DISK_MESH_BUILDER_H
#define CFEM_DISK_MESH_BUILDER_H

#include "cfemutils.h"
#include "cfem_topology.h"
#include "Vector3D.h"
#include "Mesh.h"
#include <map>
#include <vector>
#include <cmath>

namespace cfem
{
    using DiskElementType = Mesh::DiskElementType;

    /**
     * @ingroup MeshBuilders
     * @class DiskMeshBuilder
     * @brief Generates a high-quality quasi-uniform triangular mesh for a disk in 3D.
     *
     * This builder constructs the mesh using concentric rings:
     * - Radial subdivision:     \f$ \displaystyle R \simeq N_r  h \f$
     * - Angular subdivision:    \f$ \displaystyle 2\pi r \simeq N_{\theta} h \f$
     *
     * Triangles are generated with guaranteed CCW orientation to avoid inverted elements.
     * 
     * New feature: You can now specify a normal vector for the disk, allowing it to lie
     * in any plane.
     *
     * After connectivity construction:
     *  - Mesh::buildMesh() is called
     *  - Optional Laplacian smoothing improves element quality
     *  - Quadratic upgrade is supported
     */
    class DiskMeshBuilder
    {
        friend class CylinderMeshBuilder;

    private:
        struct ConstructorToken {};

    public:
        DiskMeshBuilder(ConstructorToken) {}

        static std::shared_ptr<Mesh> create(
            MeshOrder order,
            CfemReal target_mesh_size,
            const cfem::Vector3D &center,
            CfemReal radius,
            const cfem::Vector3D &normal = {0,0,1}, // Default XY plane
            DiskElementType type = DiskElementType::Triangle
        );

    private:
        void build(
            Mesh &mesh,
            MeshOrder order,
            CfemReal target_mesh_size,
            const cfem::Vector3D &center,
            CfemReal radius,
            const cfem::Vector3D &normal,
            DiskElementType type
        );

        // Helper to compute orthonormal basis (u,v) given normal
        static void computeOrthonormalBasis(const cfem::Vector3D &n, cfem::Vector3D &u, cfem::Vector3D &v)
        {
            cfem::Vector3D n_normalized = n.normalized();

            // Pick an arbitrary vector not parallel to n
            cfem::Vector3D tmp = (std::abs(n_normalized.x) < 0.9) ? cfem::Vector3D{1,0,0} : cfem::Vector3D{0,1,0};

            u = n_normalized.cross(tmp).normalized();
            v = n_normalized.cross(u);
        }
    };
}

#endif // CFEM_DISK_MESH_BUILDER_H