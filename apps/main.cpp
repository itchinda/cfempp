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
#include "DomainIntegrator.h"
#include "VTKExporter.h"
#include "VTUExporter.h"
#include "MathOperators.h"

#include <iostream>
#include <cmath>
#include <memory>
#include <tuple>

using namespace cfem;
using namespace cfem::sym;
static char help[] = "CFEM Example: Field projection on a quadratic mesh.\n\n";

int main(int argc, char **argv)
{
    //omp_set_num_threads(1);
    Environment::initialize(argc, argv, help);
    {

        std::shared_ptr<Mesh> mesh;

        {
            PROFILE_SCOPE("MeshBuilding");

            mesh = Mesh::block(
                MeshOrder::Quadratic,
                Mesh::BlockSubdivisions{10, 10, 10}, // Mesh::BlockSubdivisions{5*(2<<i), 5*(2<<i), 5*(2<<i)},
                Mesh::BlockBounds(
                    Vector3D{0., 0., 0.},
                    Vector3D{1., 1., 1.}),
                Mesh::BlockElementType::Tetrahedron);

            // mesh = Mesh::annulus(
            //     MeshOrder::Linear,
            //     0.1,
            //     {0., 0., .0},
            //     0.75,
            //     1.,
            //     Mesh::DiskElementType::Quadrilateral
            // );

            // mesh = Mesh::hollow_cylinder(
            //     MeshOrder::Quadratic,
            //     0.02,
            //     {.0, .0, .0},
            //     0.75,
            //     1.0,
            //     4.0,
            //     {0., 0., 1.},
            //     Mesh::CylinderElementType::Wedge
            //  );
        }


        std::shared_ptr<LagrangeSpace> pressureSpace, velocitySpace;
        {
            PROFILE_SCOPE("FESpaceCreation");
            pressureSpace = LagrangeSpace::create(mesh, Order::P1, FieldType::Scalar);
            velocitySpace = LagrangeSpace::create(mesh, Order::P2, FieldType::Vector3D);
            //std::tie(velocitySpace, pressureSpace) = LagrangeSpace::TaylorHood(mesh);
        }

        std::shared_ptr<FEField> p, u; 

        {
            PROFILE_SCOPE("FEFieldCreation");
            p = FEField::create("P", pressureSpace);
            u = FEField::create("U", velocitySpace);
        }

        CFEM_INFO("projection completed");
        
        auto analytical_p = scalar(x * x + y * y + z * z + x * y + x * z + z * y + x + y + z + 1); // SymEngine expressions
        auto analytical_u = vector(sin(2 * pi * x) * sin(2 * pi * y) * sin(2 * pi * z),
                                   cos(2 * pi * x) * cos(2 * pi * y) * cos(2 * pi * z),
                                   cos(2 * pi * x) * sin(2 * pi * y) * cos(2 * pi * z));

        // auto grad_p_analyt = cfem::sym::vector([](const Vector3D& X) {
        //     return Vector3D{2.0 * X.x, 2.0 * X.y, 2.0 * X.z};
        // });

        // auto analytical_p = scalar([](const Vector3D& X) {
        //     return X.x * X.x + X.y * X.y + X.z * X.z;
        // }, grad_p_analyt);

        // auto analytical_grad_u = tensor([](const Vector3D& X) -> TensorO2 {
        //     return TensorO2({2*X.x, 1., 0.,   X.z, 2*X.y, X.x, 2*X.x, 0., 1.});
        // } );
        // auto analytical_u = vector([](const Vector3D& X) {
        //     return Vector3D{X.x * X.x + X.y, X.y*X.y + X.x*X.z, X.z + X.x*X.x};
        // }, analytical_grad_u);

        {
            PROFILE_SCOPE("FEFieldProjection");
            p->project(*analytical_p);
            u->project(*analytical_u);
        }

        // auto fe_gradp = grad(p);
        // auto sin_u = 2*pow(pi*u, 2.0);
        // auto p_sq = 2*p*p;

        auto sin2 = sym::sin(2.0);

        
        {
            PROFILE_SCOPE("L2Norm");
            IntegrationSchemeParams schemes_params{};
            CfemReal l2Norm_P, l2Norm_U, h1SemiNorm_P, h1SemiNorm_U;
            DomainIntegrator::L2Norm(pressureSpace, analytical_p - p, &l2Norm_P, &h1SemiNorm_P, schemes_params);
            DomainIntegrator::L2Norm(velocitySpace, analytical_u - u, &l2Norm_U, &h1SemiNorm_U, schemes_params);
            std::cout << "L2-norm error on P: " << l2Norm_P << " H1-Seminorm error on P: " << h1SemiNorm_P << std::endl;
            std::cout << "L2-norm error on U: " << l2Norm_U << " H1-Seminorm error on U: " << h1SemiNorm_U << std::endl;

        }


        // {
        //     PROFILE_SCOPE("VTKExport");
        //     VTUExporter vtk_exporter(*mesh);
        //     vtk_exporter.setExportedEntities( {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"});
        //     vtk_exporter.exportCellQuality(true);

        //     vtk_exporter.addExpression("P", p);
        //     vtk_exporter.addExpression("P_exact", analytical_p);
        //     vtk_exporter.addExpression("U", u);
        //     vtk_exporter.addExpression("U_exact", analytical_u);
        //     // vtk_exporter.addExpression("GradP", fe_gradp);
        //     // vtk_exporter.addExpression("sin_u", sin_u);
        //     vtk_exporter.setFormat(VTUExporter::Format::BINARY);
        //     vtk_exporter.write("cube ");
        // }
    }

    Profiler::instance().saveProfiler();

    return 0;
}


