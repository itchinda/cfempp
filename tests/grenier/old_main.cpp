/**
 * @file main.cpp
 * @brief CFEM demonstration: Quadratic Mesh Generation and P2/P1 Field Projection.
 * * This example demonstrates:
 * 1. Initializing the RAII PETSc Environment.
 * 2. Building a structured Quadratic Quadrilateral mesh.
 * 3. Defining Lagrange Function Spaces (Taylor-Hood P2-P1 style).
 * 4. Projecting analytical lambdas onto Finite Element Fields.
 * 5. Exporting data to VTK for Paraview visualization.
 */

#include "cfemutils.h"
#include "Vector3D.h"
#include "Mesh.h"
#include "MeshConformityChecker.h"
#include "LagrangeSpace.h"
#include "FEField.h"
#include "VTKExporter.h"

#include <iostream>
#include <cmath>
#include <memory>

using namespace cfem;

static char help[] = "CFEM Example: Field projection on a quadratic mesh.\n\n";

int main(int argc, char **argv)
{
    cfem::Environment::initialize(argc, argv, help);
    {
        for (int i = 0; i < 1; ++i)
        {

            std::shared_ptr<Mesh> mesh;

            {
                PROFILE_SCOPE("MeshBuilding");

                mesh = Mesh::block(
                    MeshOrder::Quadratic,
                    Mesh::BlockSubdivisions{10, 10, 10}, // Mesh::BlockSubdivisions{5*(2<<i), 5*(2<<i), 5*(2<<i)},
                    Mesh::BlockBounds(
                        Vector3D{-1., -1., -1.},
                        Vector3D{1., 1., 1.}),
                    Mesh::BlockElementType::Tetrahedron);

                //mesh = Mesh::disk(MeshOrder::Quadratic, 0.25/(2<<i),  {0,0,0},  1.0, {0,1,1}, Mesh::DiskElementType::Triangle);
                // mesh = Mesh::annulus( MeshOrder::Quadratic, 0.05, Vector3D{0,0,0}, 0.75, 1.0, Mesh::DiskElementType::Triangle);
                // mesh = Mesh::cylinder(MeshOrder::Quadratic, 0.25/(2<<i), Vector3D{0.0, 0.0, 0.0}, 0.5, 2.0, Vector3D{0,1,0}, Mesh::CylinderElementType::Wedge );
                // mesh = Mesh::hollow_cylinder(MeshOrder::Quadratic, 0.025, Vector3D{0,0,0}, 0.5, 0.75, 2.0, Vector3D{0,0,1}, Mesh::CylinderElementType::Wedge);
                // mesh = Mesh::create(2, MeshOrder::Linear);
                // mesh->addVertex(0,0,0); //0 
                // mesh->addVertex(1,0,0); //1
                // mesh->addVertex(1,1,0); 
                // mesh->addVertex(0,1,0);
                // mesh->addVertex(0.5,0.5,0);

                // mesh->addCell(CellType::TRI3, {0,1,3});
                // mesh->addCell(CellType::TRI3, {4,2,3});
                // mesh->buildMesh();

                // Mesh::QuadBoundaryNames names = {"xmin", "xmax", "ymin", "ymax"};
                // Mesh::QuadCorners corners ({
                //     Vector3D{0,0,0}, Vector3D{1,0,0}, Vector3D{1,1,0}, Vector3D{0,1,0}
                // });
                // mesh =  Mesh::quad(MeshOrder::Linear, Mesh::QuadSubdivisions(2,2), corners, Mesh::QuadElementType::Quadrilateral, names);

                // Mesh::BlockCorners corners ({
                //                         Vector3D{0,0,0}, Vector3D{2,0,0}, Vector3D{2,2,0}, Vector3D{0,2,0},
                //                         Vector3D{0,0,2}, Vector3D{2,0,2}, Vector3D{2,2,2}, Vector3D{0,2,2}
                //                     });
                // mesh = Mesh::block(MeshOrder::Linear, Mesh::BlockSubdivisions(3,3,3), corners, Mesh::BlockElementType::Hexahedron);

                // std::cout << mesh->getMeshFormat() << std::endl;
                //std::cout << MeshConformityChecker::check(*mesh, true).toString() << std::endl;

                // mesh = Mesh::loadFromGmsh("../tests/grenier/mail.msh", Vector3D{0,0,1});


            }

            //std::cout << MeshConformityChecker::check(*mesh, true).toString() << std::endl;

            std::shared_ptr<LagrangeSpace> p1Space, p2Space;
            {
                PROFILE_SCOPE("FESpaceCreation");
                p1Space = LagrangeSpace::create(mesh, Order::P2, FieldType::Scalar);   // Pressure (Linear)
                p2Space = LagrangeSpace::create(mesh, Order::P1, FieldType::Vector3D); // Velocity (Quadratic)
            }

            std::shared_ptr<FEField> pressure, velocity, un;

            {
                PROFILE_SCOPE("FEFieldCreation");
                pressure = FEField::create("P", p1Space);
                un = FEField::create("Un", p1Space);
                velocity = FEField::create("U", p2Space);
            }

            const CfemReal PI = cfem::pi;

            auto analytical_p = [&](const Vector3D &X) -> CfemReal
            {
                return X.x * X.x + X.y * X.y + X.z * X.z + X.x * X.y + X.x * X.z + X.z * X.y + X.x + X.y + X.z + 1;
            };

            auto analytical_u = [&](const Vector3D &X) -> Vector3D
            {
                return {
                    std::sin(2 * PI * X.x) * std::sin(2 * PI * X.y) * std::sin(2 * PI * X.z),
                    std::cos(2 * PI * X.x) * std::cos(2 * PI * X.y) * std::cos(2 * PI * X.z),
                    std::cos(2 * PI * X.x) * std::sin(2 * PI * X.y) * std::cos(2 * PI * X.z)};
            };

            auto analytical_un = [&](const Vector3D &X) -> CfemReal
            {
                return 1.0;
            };

            {
                PROFILE_SCOPE("FEFieldProjection");
                pressure->project(analytical_p);
                velocity->project(analytical_u);
                un->project(analytical_un);
            }

            IntegrationSchemeParams schemes_params{};
            schemes_params.NumIntgPtTriElem = 16;
            {
                PROFILE_SCOPE("L2NormDistance");
                CfemReal lNorm_l2_p{FEField::L2NormDistance(*pressure, analytical_p, schemes_params)};
                CfemReal lNorm_l2_u{FEField::L2NormDistance(*velocity, analytical_u, schemes_params)};
                std::cout << "L2 error on the pressure: " << lNorm_l2_p  <<", L2 error on the velocity: " << lNorm_l2_u  << std::endl;
            //     CfemReal lVolume_sqrt{un->normL2(schemes_params)};
            //     CfemReal lVolume {lVolume_sqrt*lVolume_sqrt};
            //     std::cout << "Mesh size: " << 0.5/(2<<i) << ", Volume: " << lVolume  << ", L2 Error: " << std::abs(lVolume - cfem::pi*1.0*1.0) << std::endl;
            }

            // VTKExporter vtk_exporter(*mesh);
            // vtk_exporter.setExportedEntities( {"Lateral", "Top", "Bottom" });
            // vtk_exporter.setExportQuality(true);

            // vtk_exporter.addField(*pressure);
            // vtk_exporter.addField(*velocity);
            // vtk_exporter.write("cube.vtk");

            //std::cout << mesh->getMeshFormat() << std::endl;
        }
    }

    Profiler::instance().saveProfiler();

    return 0;
}