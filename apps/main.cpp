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
#include <omp.h>
#include <span>

using namespace cfem;
using namespace cfem::sym;
using namespace cfem::sym::symengine;
static char help[] = "CFEM Example: Field projection on a quadratic mesh.\n\n";

int main(int argc, char **argv)
{
    // omp_set_num_threads(1);
    Environment::initialize(argc, argv, help);
    {

        std::shared_ptr<Mesh> mesh;

        {
            PROFILE_SCOPE("MeshBuilding");

            // mesh = Mesh::quad(
            //     MeshOrder::Linear,
            //     Mesh::QuadSubdivisions{1, 1}, // Mesh::BlockSubdivisions{5*(2<<i), 5*(2<<i), 5*(2<<i)},
            //     Mesh::QuadBounds(
            //         Vector3D{0., 0., 0.},
            //         Vector3D{1., 1., 0.}
            //     ),
            //     Mesh::QuadElementType::TriangleCrissCross
            // );

            mesh = Mesh::block(
                MeshOrder::Linear,
                Mesh::BlockSubdivisions{10, 10, 10}, // Mesh::BlockSubdivisions{5*(2<<i), 5*(2<<i), 5*(2<<i)},
                Mesh::BlockBounds(
                    Vector3D{-1., -1., -1.},
                    Vector3D{1., 1., 1.}
                ),
                Mesh::BlockElementType::Tetrahedron );

        }

        auto Q = LagrangeSpace::create(mesh, Order::P2, FieldType::Scalar);
        auto V = LagrangeSpace::create(mesh, Order::P2, FieldType::Vector3D);

        auto p = Trial(Q, "p");
        auto u = Trial(V, "u");
        auto q = Test (Q, "q");
        auto v = Test (V, "v");

        auto analytical_p = scalar(x * x + y * y +  x + y); 
        auto analytical_u = vector( x*x+y*y+z*z,
                                   -2*x*y,
                                    x*z);

        // auto analytical_p = scalar(x * x + y * y + z * z + y + x + y + z + 1); // SymEngine expressions
        // auto analytical_u = vector(sin(2 * pi * x ) * sin(2 * pi * z ),
        //                            cos(2 * pi * y ) * cos(2 * pi * x ),
        //                            cos(2 * pi * y ) * cos(2 * pi * z ));

        {
            PROFILE_SCOPE("FieldProjection");
            u->project(analytical_u);
            p->project(analytical_p);
        }

        {
            PROFILE_SCOPE("L2Norm3");
            IntegrationSchemeParams schemes_params{};
            DomainIntegrator integrator(mesh, schemes_params);
            // schemes_params.NumIntgPtHexElem = 1;
            CfemReal l2Norm_P, l2Norm_U, h1SemiNorm_P, h1SemiNorm_U;
            integrator.L2Norm(p - analytical_p, &l2Norm_P, &h1SemiNorm_P);
            integrator.L2Norm(u - analytical_u, &l2Norm_U, &h1SemiNorm_U);
            CFEM_INFO << "L2-norm error on P: " << l2Norm_P << " H1-Seminorm error on P: " << h1SemiNorm_P;
            CFEM_INFO << "L2-norm error on U: " << l2Norm_U << " H1-Seminorm error on U: " << h1SemiNorm_U;



            CFEM_INFO << "sqrt (int_{Omega} |grad(u-analytical_u)|^2 ds ) = " << sqrt(integrator.computeIntegral(mag_sqr(grad(u-analytical_u)), "ymin"));
            // CFEM_INFO << "surface error = " << sqrt( integrator.computeIntegral(mag_sqr(scalar(y)), "xmin") );
            // CFEM_INFO << "surface error = " << sqrt( integrator.computeIntegral(mag_sqr(grad(u - analytical_u)), "xmin") );

        }

    }

    Profiler::instance().saveProfiler();

    return 0;
}

// #include "cfemutils.h"
// #include "Mesh.h"
// #include "FEField.h"
// #include "MathOperators.h"
// #include "DomainIntegrator.h"
// #include "VTUExporter.h"

// using namespace cfem;
// using namespace cfem::sym;
// using namespace cfem::sym::symengine;

// int main(int argc, char **argv) {
//     // 1. One call to rule them all: Handles MPI, PETSc, and CFEM internal state
//     cfem::Environment::initialize(argc, argv);

//     // 2. Generate Mesh
//     auto mesh = Mesh::block(MeshOrder::Linear, 
//                             {10, 10, 10}, 
//                             Mesh::BlockBounds({0,0,0}, {1,1,1}), 
//                             Mesh::BlockElementType::Tetrahedron);

//     // 3. Setup P2 Vector Space and Field
//     auto V = LagrangeSpace::create(mesh, Order::P2, FieldType::Vector3D);
//     auto u = FEField::create("U", V);
        
//     // 4. JIT-compiled Symbolic Definition
//     auto analytical_p = scalar(x * x + y * y + z * z + x * y + x * z + z * y + x + y + z + 1); // SymEngine expressions
//     auto analytical_u = vector(sin(2 * pi * x), cos(2 * pi * y), sin(pi * x * y * z));

//     // 5. High-level Projection
//     u->project(analytical_u);
//     DomainIntegrator integrator(mesh, IntegrationSchemeParams{});
//     CfemReal l2Norm_U, h1SemiNorm_U;
//     integrator.L2Norm(analytical_u-u, &l2Norm_U, &h1SemiNorm_U);
//     auto expr = abs( inner(grad(u), grad(u)) - inner(grad(analytical_u), grad(analytical_u)) );
//     CFEM_INFO << "Domain volume: " << integrator.computeIntegral( expr );
        
//     // CFEM_INFO << "L2-norm error on U: " << l2Norm_U;
//     // CFEM_INFO << "H1-semi norm error on U: " << h1SemiNorm_U;
//     // 6. Error Analysis (Clean syntax, no explicit scheme params needed)
//     // 7. High-performance Binary Export
//     VTUExporter vtu_exporter(*mesh);
//     vtu_exporter.addExpression("U", u);
//     vtu_exporter.setFormat(VTUExporter::Format::BINARY);
//     vtu_exporter.write("simulation_results");

//     return 0; 
//     // cfem::Environment handles finalization automatically on scope exit or via atexit
// }

// #include "cfemutils.h"
// #include "Vector3D.h"
// #include "Mesh.h"
// #include "MeshConformityChecker.h"
// #include "LagrangeSpace.h"
// #include "FEField.h"
// #include "DomainIntegrator.h"
// #include "VTKExporter.h"
// #include "VTUExporter.h"
// #include "MathOperators.h"

// #include <iostream>
// #include <cmath>
// #include <memory>
// #include <tuple>

// using namespace cfem;
// using namespace cfem::sym;
// static char help[] = "CFEM Example: Field projection on a quadratic mesh.\n\n";

// int main(int argc, char **argv)
// {
//     Environment::initialize(argc, argv, help);
//     {

//         std::shared_ptr<Mesh> mesh;

//         {
//             PROFILE_SCOPE("MeshBuilding");

//             mesh = Mesh::block(
//                 MeshOrder::Quadratic,
//                 Mesh::BlockSubdivisions{100, 100, 100}, // Mesh::BlockSubdivisions{5*(2<<i), 5*(2<<i), 5*(2<<i)},
//                 Mesh::BlockBounds(
//                     Vector3D{0., 0., 0.},
//                     Vector3D{1., 1., 1.}),
//                 Mesh::BlockElementType::Tetrahedron);
//         }

//         auto u_raw = scalar( sin(pi * x) * exp(-x*x) ); 
//         auto u2 = cache(u_raw);
//         auto f_u = u2 + sin(u2);       
//         auto g_u = u_raw + sin(u_raw);

//         {
//             PROFILE_SCOPE("UnCachedExport");
//             VTUExporter vtk_exporter(*mesh);
//             vtk_exporter.addExpression("g_u", g_u);
//             vtk_exporter.setFormat(VTUExporter::Format::BINARY);
//             vtk_exporter.write("cube ");
//         }

//         {
//             PROFILE_SCOPE("CachedExport");
//             VTUExporter vtk_exporter(*mesh);
//             vtk_exporter.addExpression("f_u", f_u);
//             vtk_exporter.setFormat(VTUExporter::Format::BINARY);
//             vtk_exporter.write("cube ");
//         }
//     }

//     Profiler::instance().saveProfiler();

//     return 0;
// }


// /**
//  * @file main.cpp
//  * @brief CFEM demonstration: Benchmark for evaluating expressions in Raw, Cached.
//  */

// #include "cfemutils.h"
// #include "Vector3D.h"
// #include "Mesh.h"
// #include "VTUExporter.h"
// #include "MathOperators.h"
// #include <iostream>
// #include <cmath>
// #include <memory>

// using namespace cfem;
// using namespace cfem::sym;
// using namespace cfem::sym::symengine;

// static char help[] = "CFEM Example: Benchmark expression evaluation modes (Raw, Cached, Compiled).\n\n";

// int main(int argc, char **argv)
// {
//     Environment::initialize(argc, argv, help);

//     // Create mesh for testing
//     std::shared_ptr<Mesh> mesh;
//     {
//         PROFILE_SCOPE("MeshBuilding");
//         mesh = Mesh::block(
//             MeshOrder::Quadratic,
//             Mesh::BlockSubdivisions{100, 100, 100},
//             Mesh::BlockBounds(Vector3D{0., 0., 0.},
//                               Vector3D{1., 1., 1.}),
//             Mesh::BlockElementType::Tetrahedron);
//     }

//     // auto Y = sin(x) - sin(x);

//     // std::cout << std::endl;
//     // std::cout << "Y = " << Y->toString() << std::endl;
//     // std::cout << std::endl;
    
//     // Define expression for benchmarking
//     auto raw = sym::pow(sin(pi * x), 4.0)* exp(-x * x) * log(1.0 + x * x) * sqrt(1.0 + y * y);

//     CFEM_INFO << "raw = "<< scalar(raw)->toString();

//     // Create repeated subtree for testing
//     auto subtree =
//         raw + sin(raw) + cos(raw) + exp(raw) + sqrt(abs(raw));

//     std::cout << "subtree = "<< scalar(subtree)<< std::endl;

//     // BIG expression (repeated subtree)
//     auto big = scalar(subtree + subtree + subtree + subtree + subtree + subtree + subtree + subtree);
//     CFEM_INFO << "big = "<< big;

//     // Benchmark Raw Evaluation
//     {
//         PROFILE_SCOPE("RawExport");
//         VTUExporter vtk_exporter(*mesh);

//         // Benchmark: Raw expression
//         vtk_exporter.addExpression("raw", big);
//         vtk_exporter.setFormat(VTUExporter::Format::BINARY);
//         vtk_exporter.write("cube_raw");
//     }

//     // Benchmark Cached Evaluation
//     {
//         PROFILE_SCOPE("CachedExport");

//         // Cache the repeated subtree
//         auto subtree_cached = scalar(subtree);

//         // Benchmark: Cached expression
//         VTUExporter vtk_exporter(*mesh);
//         vtk_exporter.addExpression("cached", 
//             subtree_cached + subtree_cached + subtree_cached + subtree_cached +
//             subtree_cached + subtree_cached + subtree_cached + subtree_cached);
//         vtk_exporter.setFormat(VTUExporter::Format::BINARY);
//         vtk_exporter.write("cube_cached");
//     }

//     auto u = sin(x)/tan(x);
//     auto v = vector(cos(z), y, log(tan(z))*z);
//     CFEM_INFO << v->diff(2);

//     // Output Profiler results
//     Profiler::instance().saveProfiler();

//     return 0;
// }

