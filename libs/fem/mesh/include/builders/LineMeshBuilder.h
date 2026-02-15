#ifndef CFEM_LINE_MESH_BUILDER_H
#define CFEM_LINE_MESH_BUILDER_H

#include "cfemutils.h"
#include "cfem_topology.h"
#include "Vector3D.h"
#include "Mesh.h"


namespace cfem
{

    /**
     * @class LineMeshBuilder
     * @brief 
     * 
     */
    class LineMeshBuilder
    {
    private:
        struct ConstructorToken{};

    public:
        LineMeshBuilder(ConstructorToken) {};
        static std::shared_ptr<Mesh> create(MeshOrder order,
                                            CfemInt num_divisions,
                                            const cfem::Vector3D &starting_point,
                                            const cfem::Vector3D &end_point);

    private:
        void build(Mesh &mesh, 
                   MeshOrder order, 
                   CfemInt num_divisions, 
                   const cfem::Vector3D &starting_position, 
                   const cfem::Vector3D &ending_position);
    };
}

#endif