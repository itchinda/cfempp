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

#ifndef CFEM_DISK_MESH_BUILDER_H
#define CFEM_DISK_MESH_BUILDER_H

#include <cmath>
#include <map>
#include <vector>

#include "cfem_topology.h"
#include "cfemutils.h"
#include "Mesh.h"
#include "Vector3D.h"

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