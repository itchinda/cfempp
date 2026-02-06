#include <petsc.h>
#include "PetscMatrix.h"
#include "PetscVector.h"
#include "PetscKSPSolver.h"

#include "cfemutils.h"
#include "Mesh.h"
#include "Vector3D.h"
#include "VTKExporter.h"
#include "GmshParser.h"

#include <iostream>

int main(int argc, char **argv) {

    //-----------------------------------------------------------------------------
    //------------------------------ Mesh class -----------------------------------
    Mesh lMail;

    // LineMeshBuilder lLine;
    // lLine.build(lMail, 10, Vector3D{1.0, 0.0,0.0}, Vector3D{0.5, 0, 1});

    // QuadMeshBuilder lRectangle;
    // lRectangle.build(lMail, 1, 1, {Vector3D{0,0,0}, Vector3D{1,0,0}, Vector3D{1,1,0}, Vector3D{0,0.5,0} }, QuadMeshBuilder::SubdivisionType::CrissCross);
    // lMail.upgradeToQuadratic();

    // CubicMeshBuilder lCube;
    // CfemInt N {5};
    // CfemInt order {2};
    // lCube.build(lMail, order, N, N, N, {Vector3D{0,0,0}, Vector3D{1,0,0}, Vector3D{1,1,0}, Vector3D{0,1,0},
    //                                 Vector3D{0,0,1}, Vector3D{1,0,1}, Vector3D{1,1,1}, Vector3D{0,1,1} },
    //                                 CubicMeshBuilder::SubdivisionType::None);

    // DiskMeshBuilder lDisk;
    // double radius {1.0};
    // double target_mesh_size {radius / 1000};
    // CfemInt order {1};
    // lDisk.build(lMail, order, target_mesh_size, Vector3D{0,0,0}, radius, DiskMeshBuilder::ElementType::Triangle);
    
    // SphereMeshBuilder lSphere;
    // double radius {1.0};
    // double target_mesh_size {radius / 50};
    // CfemInt order {1};
    // lSphere.build(lMail, order, target_mesh_size, Vector3D{0,0,0}, radius, SphereMeshBuilder::ElementType::Tetrahedron);

    GmshParser lMshParser;
    lMshParser.load(lMail, "../tests/grenier/hex.msh");
    //lMail.downgradeToLinear();
    VTKExporter vtk_exporter (lMail);
    vtk_exporter.setExportedEntities({
        "Bottom", "Top", "Front", "Right", "Back", "Left", "Bord", "SYS_BOUNDARY"
    } );

    vtk_exporter.setExportQuality(true);
    vtk_exporter.write("bunny.vtk");

    // auto report = MeshQualityTool::analyze(lMail);
    // MeshQualityTool::printSummary(report);

    std::cout << lMail.info() << std::endl;

    // for (const auto& ent : lMail.getEntities()){
    //     std::cout << ent->getEntityFormat() << std::endl;
    // }

    //-------------------------------------------------------------------------------
    //----------------------------- Vector3D struct --------------------------------- cx
    // Vector3D v1 {1.0, 2.0, 3.0};
    // Vector3D v2 {3.0, 2.0, 1.0};
    // Vector3D v3 = v1+v2;
    // Vector3D v4 = 4.0*v2;
    // Vector3D v5 = v2*4.0;

    // std::cout << "v1 = " << v1 << std::endl;
    // std::cout << "v2 = " << v2 << std::endl;
    // std::cout << "v1+v2 = " << v3 << std::endl;
    // std::cout << "4*v2 = " << v4 << std::endl;
    // std::cout << "v2*4 = " << v5 << std::endl;

    // std::cout << "norm(v1) = " << norm(v1) << std::endl;

    // std::cout << "dot(v1,v2) = " << dot(v1,v2) << std::endl;
    // std::cout << "v1.dot(v2) = " << v1.dot(v2) << std::endl;

    // std::cout << "cross(v1,v2) = " << cross(v1,v2) << std::endl;
    // std::cout << "v2.cross(v1) = " << v2.cross(v1) << std::endl;
    //---------------------------------------------------------------------------------





    return 0;


}

