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

#ifndef CFEM_QUAD_MESH_BUILDER_H
#define CFEM_QUAD_MESH_BUILDER_H

#include "cfemutils.h"
#include "Vector3D.h"
#include "cfem_topology.h"
#include "Mesh.h"

namespace cfem{
    using QuadElementType   = Mesh::QuadElementType;   /// Alias to Mesh element type for 2D quad meshes.
    using QuadCorners       = Mesh::QuadCorners;       /// 4 ordered corner points defining the quadrilateral domain. The ordering must be consistent (counter-clockwise recommended).
    using QuadSubdivisions  = Mesh::QuadSubdivisions;  /// Number of structured subdivisions in x and y directions.
    using QuadBoundaryNames = Mesh::QuadBoundaryNames; /// Optional names assigned to the four boundaries: xmin, xmax, ymin, ymax.


/**
 * @ingroup MeshBuilders
 * @brief Structured 2D quadrilateral mesh generator.
 *
 * QuadMeshBuilder generates structured meshes over a general
 * convex quadrilateral domain defined by four corner points.
 *
 * The geometry is defined through bilinear interpolation between
 * the provided corners. The resulting mesh is logically structured
 * (nx × ny grid) but geometrically mapped to an arbitrary quad.
 *
 * Supported element types:
 *  - Quadrilateral elements
 *  - Triangular elements (via structured splitting)
 *
 * The builder:
 *  - Generates vertices
 *  - Creates cells according to the chosen MeshOrder
 *  - Assigns boundary tags
 *
 * This class cannot be instantiated directly.
 * Use the static create() factory methods.
 */
class QuadMeshBuilder {
private:
    struct ConstructorToken{}; /// Tag type used to restrict direct construction. 

public:

    QuadMeshBuilder(ConstructorToken){}; /// Internal constructor (factory pattern enforcement).

    /**
     * @brief Create a structured quadrilateral mesh.
     *
     * @param order Polynomial order of the mesh (P1, P2, ...).
     * @param subDiv Number of subdivisions in x and y.
     * @param corners Geometric corners of the domain.
     * @param elemType Type of elements (Quad or Tri split).
     *
     * @return Shared pointer to the generated Mesh.
     *
     * Boundaries receive default names.
     */
    static std::shared_ptr<Mesh> create(MeshOrder order,
                                        const QuadSubdivisions& subDiv,
                                        const QuadCorners& corners,
                                        QuadElementType elemType = QuadElementType::Quadrilateral) 
    {
        return QuadMeshBuilder::create(order, subDiv, corners, elemType, QuadBoundaryNames{});
    }

    static std::shared_ptr<Mesh> create(MeshOrder order,
                                        const QuadSubdivisions subDiv,
                                        const QuadCorners& corners,
                                        QuadElementType elemType, 
                                        const QuadBoundaryNames& names );
private:
    void applyBoundaryTags(Mesh& mesh, 
                           CfemInt nx, 
                           CfemInt ny, 
                           const std::vector<CfemInt>& grid, 
                           const QuadBoundaryNames& names);
                           
    void build(Mesh&                          mesh, 
               MeshOrder                      order,
               const QuadSubdivisions  &subDiv,
               const QuadCorners       &corners,
               QuadElementType          elemType,
               const QuadBoundaryNames &names);
};


}

#endif