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