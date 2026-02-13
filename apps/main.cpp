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
    // Initializes PETSc/MPI via RAII. Finalization is handled automatically at exit.
    cfem::Environment::initialize(argc, argv, help);
    {
        PROFILE_SCOPE("Main");
        
        // Creates a [-1, 1] x [-1, 1] square mesh with 10,000 quadratic elements.
        // auto mesh = QuadMeshBuilder::create(
        //     MeshOrder::Linear,
        //     QuadMeshBuilder::Subdivisions{2000, 2000},
        //     QuadMeshBuilder::Corners{
        //         Vector3D{-1., -1., 0.}, // Bottom-Left
        //         Vector3D{ 1., -1., 0.}, // Bottom-Right
        //         Vector3D{ 1.,  1., 0.}, // Top-Right
        //         Vector3D{-1.,  1., 0.}  // Top-Left
        //     },
        //     QuadMeshBuilder::ElementType::CrissCross
        // );

        auto mesh = DiskMeshBuilder::create(
            MeshOrder::Linear,
            1.0/1000,
            Vector3D{0., 0., 0.},
            1.0,
            DiskMeshBuilder::ElementType::Triangle
        );

        // Define the polynomial basis on the mesh.
        auto p1Space = LagrangeSpace::create(mesh, Order::P1, FieldType::Scalar);   // Pressure (Linear)
        auto p2Space = LagrangeSpace::create(mesh, Order::P2, FieldType::Vector3D); // Velocity (Quadratic)

        // Associate physical fields with their respective function spaces.
        auto pressure = FEField::create("P", p1Space);
        auto velocity = FEField::create("U", p2Space);

        // Define exact solutions and project them onto the FE basis.
        const CfemReal PI = cfem::pi;

        auto analytical_p = [&](const Vector3D &X) -> CfemReal {
            // return std::sin(2 * PI * X.x) * std::cos(2 * PI * X.y);
            return X.x * X.x + X.y * X.y + X.x*X.y + X.x + X.y + 1;
        };

        auto analytical_u = [&](const Vector3D &X) -> Vector3D {
            return {
                std::sin(2 * PI * X.x) * std::sin(2 * PI * X.x), 
                std::cos(2 * PI * X.y) * std::cos(2 * PI * X.y), 
                0.0
            };
        };

        pressure->project(analytical_p);
        velocity->project(analytical_u);

        //std::cout << mesh->info() << std::endl;

        #include "Quadrature.h"

        IntegrationSchemeParams schemes_params{};
        schemes_params.NumIntgPtTriElem = 12;

        CfemReal lNorm_l2_p {FEField::L2NormDistance(*pressure, analytical_p, schemes_params)};
        CfemReal lNorm_l2_u {FEField::L2NormDistance(*velocity, analytical_u, schemes_params)};
        std::cout << "L2 Error norm for P: "  << lNorm_l2_p << std::endl;
        std::cout << "L2 Error norm for U: "  << lNorm_l2_u << std::endl;

        // Prepare the VTK exporter and add fields for visualization.
        VTKExporter vtk_exporter(*mesh);

        // // Define which boundary entities to include in the output
        vtk_exporter.setExportedEntities({"Bottom", "Top", "Right", "Left", "SYS_BOUNDARY"});

        vtk_exporter.addField(*pressure);
        vtk_exporter.addField(*velocity);
        vtk_exporter.write("square.vtk");

        // // Print mesh statistics and entity formatting for verification.
        std::cout << mesh->info() << std::endl;

        // for (const auto &ent : mesh->getEntities()) {
        //     std::cout << "Entity: " << ent->getEntityFormat() << std::endl;
        // }
    }

    Profiler::instance().saveAndPlot();

    return 0;
}