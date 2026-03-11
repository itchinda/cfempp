#ifndef CFEM_BLOCK_MESH_BUILDER_H
#define CFEM_BLOCK_MESH_BUILDER_H

#include "cfemutils.h"
#include "Vector3D.h"
#include "cfem_topology.h"
#include "Mesh.h"

namespace cfem
{
    using BlockElementType   = Mesh::BlockElementType;
    using BlockCorners       = Mesh::BlockCorners;
    using BlockSubdivisions  = Mesh::BlockSubdivisions;
    using BlockBoundaryNames = Mesh::BlockBoundaryNames;


    /**
     * @class BlockMeshBuilder
     * @brief Structured 3D block mesh generator.
     *
     * Builds a structured hexahedral or tetrahedral mesh over a
     * general trilinear block defined by eight corner vertices.
     *
     * The geometry is constructed via trilinear interpolation of the
     * provided corner coordinates.
     *
     * Supported element types:
     *  - Hexahedron (HEX8)
     *  - Tetrahedron (6 TET4 per voxel)
     *
     * The mesh can optionally be upgraded to quadratic order.
     *
     * Boundary faces are automatically tagged into six entities:
     *  - xmin, xmax
     *  - ymin, ymax
     *  - zmin, zmax
     *
     * This builder is intended for structured meshes and benchmarking.
     * 
     * @ingroup MeshBuilders
     */
    class BlockMeshBuilder
    {
    private:
        /// Internal token to prevent direct public construction.
        struct ConstructorToken {};

    public:

        /**
         * @brief Internal constructor (not intended for direct use).
         */
        BlockMeshBuilder(ConstructorToken) {}

        /**
         * @brief Create a structured 3D block mesh.
         *
         * @param order     Desired polynomial order of the mesh (Linear or Quadratic).
         * @param subDiv    Number of subdivisions in each coordinate direction.
         * @param corners   Eight corner vertices defining the block geometry.
         * @param elemType  Element type (Hexahedron or Tetrahedron).
         *
         * @return Shared pointer to the constructed Mesh.
         *
         * @note Boundary entities will be created with default names.
         */
        static std::shared_ptr<Mesh> create(
            MeshOrder order,
            const BlockSubdivisions &subDiv,
            const BlockCorners &corners,
            BlockElementType elemType = BlockElementType::Hexahedron
        )
        {
            return create(order, subDiv, corners, elemType, BlockBoundaryNames{});
        }

        /**
         * @brief Create a structured 3D block mesh with custom boundary names.
         *
         * @param order     Desired polynomial order.
         * @param subDiv    Number of subdivisions (nx, ny, nz).
         * @param corners   Eight corner vertices defining the block.
         * @param elemType  Element type (HEX8 or TET4).
         * @param names     Custom boundary entity names.
         *
         * @return Shared pointer to the constructed Mesh.
         */
        static std::shared_ptr<Mesh> create(
            MeshOrder order,
            const BlockSubdivisions &subDiv,
            const BlockCorners &corners,
            BlockElementType elemType,
            const BlockBoundaryNames &names
        );

    private:

        /**
         * @brief Build the mesh topology and geometry.
         *
         * Generates vertices via trilinear interpolation and constructs
         * structured cells according to the chosen element type.
         *
         * @param mesh      Target mesh object.
         * @param order     Polynomial order.
         * @param subDiv    Subdivision counts.
         * @param corners   Block corner geometry.
         * @param elemType  Element type.
         * @param names     Boundary naming structure.
         */
        void build(
            Mesh &mesh,
            MeshOrder order,
            const BlockSubdivisions &subDiv,
            const BlockCorners &corners,
            BlockElementType elemType,
            const BlockBoundaryNames &names
        );

        void build_omp(
            Mesh &mesh,
            MeshOrder order,
            const BlockSubdivisions &subDiv,
            const BlockCorners &corners,
            BlockElementType elemType,
            const BlockBoundaryNames &names
        );

        /**
         * @brief Assign boundary facets to named entities.
         *
         * Detects boundary facets using structured indexing and
         * assigns them to their corresponding boundary entity.
         *
         * @param mesh   Target mesh.
         * @param nx     Subdivisions in x-direction.
         * @param ny     Subdivisions in y-direction.
         * @param nz     Subdivisions in z-direction.
         * @param grid   Structured grid vertex ID mapping.
         * @param names  Boundary entity names.
         */
        void applyBoundaryTags(
            Mesh &mesh,
            CfemInt nx,
            CfemInt ny,
            CfemInt nz,
            const std::vector<CfemInt> &grid,
            const BlockBoundaryNames &names
        );

        void applyBoundaryTags_omp(
            Mesh &mesh,
            CfemInt nx,
            CfemInt ny,
            CfemInt nz,
            const std::vector<CfemInt> &grid,
            const BlockBoundaryNames &names
        );
    };

} // namespace cfem

#endif
