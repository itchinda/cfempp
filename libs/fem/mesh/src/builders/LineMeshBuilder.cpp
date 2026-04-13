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

#include "LineMeshBuilder.h"
#include "Mesh.h"

namespace cfem{

    /**
     * @brief 
     * 
     * @param order 
     * @param num_divisions 
     * @param starting_point 
     * @param end_point 
     * @return std::shared_ptr<Mesh> 
     */
    std::shared_ptr<Mesh> LineMeshBuilder::create(MeshOrder order,
                                                  CfemInt num_divisions,
                                                  const cfem::Vector3D &starting_point,
                                                  const cfem::Vector3D &end_point)
    {
        auto mesh = Mesh::create();

        LineMeshBuilder builder(ConstructorToken{});

        builder.build(*mesh, order, num_divisions, starting_point, end_point);

        return mesh;
    }


    /**
     * @brief 
     * 
     * @param mesh 
     * @param order 
     * @param num_divisions 
     * @param starting_point 
     * @param ending_point 
     */
    void LineMeshBuilder::build(Mesh &mesh, 
                                MeshOrder order, 
                                CfemInt num_divisions, 
                                const cfem::Vector3D &starting_point, 
                                const cfem::Vector3D &ending_point)
    {
        mesh.setOrder(MeshOrder::Linear);
        mesh.setDimension(1);

        for (CfemInt i = 0; i <= num_divisions; ++i)
        {
            CfemReal t = static_cast<CfemReal>(i) / num_divisions;
            cfem::Vector3D p = starting_point * (1.0 - t) + ending_point * t;
            mesh.addVertex(p.x, p.y, p.z);
        }

        for (int i = 0; i < num_divisions; ++i)
        {
            mesh.addCell(CellType::LINE2, {(CfemInt)i, (CfemInt)i + 1});
        }

        mesh.buildMesh();
        if (order == MeshOrder::Quadratic)
            mesh.upgradeToQuadratic(); // Becomes LINE3
    }

}