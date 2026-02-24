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
            }

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
            }

            VTKExporter vtk_exporter(*mesh);
            vtk_exporter.setExportedEntities( {"Lateral", "Top", "Bottom" });
            vtk_exporter.setExportQuality(true);

            vtk_exporter.addField(*pressure);
            vtk_exporter.addField(*velocity);
            vtk_exporter.write("cube.vtk");

    }

    Profiler::instance().saveProfiler();

    return 0;
}