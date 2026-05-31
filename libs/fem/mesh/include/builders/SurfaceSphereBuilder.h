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

#ifndef CFEM_SURFACE_SPHERE_BUILDER_H
#define CFEM_SURFACE_SPHERE_BUILDER_H

#include "cfemutils.h"
#include "Mesh.h"
#include "Vector3D.h"

namespace cfem
{

    /**
    * @ingroup MeshBuilders
    * @brief Surface mesh generator for a sphere (Icosphere).
    *
    * This builder generates a triangular surface mesh of a sphere using recursive 
    * subdivision of a base polyhedron (typically an icosahedron or octahedron). 
    * After each subdivision, new vertices are radially projected onto the exact 
    * spherical surface of radius R.
    *
    * Unlike UV-spheres (lat/long grids), this method produces a quasi-uniform 
    * distribution of triangles, avoiding the "pinched" elements at the poles.
    *
    * The generated mesh is suitable for:
    * - Surface transport problems.
    * - Shell and membrane finite element simulations.
    * - Spherical boundary conditions in acoustics or electromagnetics.
    */
    class SurfaceSphereBuilder
    {
    private:
        struct ConstructorToken { };

    public:
    
        SurfaceSphereBuilder(ConstructorToken) {}

        /**
         * @brief Creates a spherical surface mesh.
         * @param order           Polynomial order (Linear TRI3 or Quadratic TRI6).
         * @param target_mesh_size The desired approximate size of the elements.
         * @param center          Center coordinates of the sphere.
         * @param radius          Radius of the sphere.
         *
         * @return Shared pointer to the generated Mesh object.
         * * @note The number of subdivision levels is automatically calculated to 
         * match the requested target_mesh_size as closely as possible.
         */
        static std::shared_ptr<Mesh> create(MeshOrder order,
                                            CfemReal target_mesh_size,
                                            const cfem::Vector3D &center,
                                            CfemReal radius);

    private:
        /**
         * @brief Internal recursive construction algorithm.
         *
         * The algorithm follows these steps:
         * 1. Initialize a base unit polyhedron (e.g., Icosahedron).
         * 2. Calculate the required subdivision depth based on target_mesh_size.
         * 3. Subdivide each triangular face into 4 smaller triangles.
         * 4. Project new mid-edge nodes onto the sphere: P_new = C + (P_mid - C).norm() * R.
         * 5. If quadratic order is requested, mid-edge nodes are preserved and mapped.
         */
        static void build(Mesh &mesh,
                          MeshOrder order,
                          CfemReal target_mesh_size,
                          const cfem::Vector3D &center,
                          CfemReal radius);
    };

}

#endif