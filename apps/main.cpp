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



// #include "cfemutils.h"
// #include "Mesh.h"
// #include "FEField.h"
// #include "MathOperators.h"
// #include "ScalarFunctionalEvaluator.h"
// #include "LagrangeSpace.h"
// #include "VTUExporter.h"

// using namespace cfem;
// using namespace cfem::sym;
// using namespace cfem::sym::symengine;

// int main(int argc, char **argv) {
//     // One call to rule them all: Handles MPI, PETSc, and CFEM internal state
//     cfem::Environment::initialize(argc, argv);

//     // Generate Mesh
//     auto mesh = Mesh::block(MeshOrder::Linear, 
//                             {10, 10, 10}, 
//                             Mesh::BlockBounds({0,0,0}, {1,1,1}), 
//                             Mesh::BlockElementType::Tetrahedron);

//     // Setup P2 Vector Space and Field
//     auto V = LagrangeSpace::create(mesh, Order::P2, FieldType::Vector3D);
//     auto u = FEField::create("u", V);
        
//     auto analytical_u = vector(sin(2 * CFEM_PI * x), cos(2 * CFEM_PI * y), sin(CFEM_PI * x * y * z));
//     u->interpolate(analytical_u);

//     // Error Analysis (Clean syntax, no explicit scheme params needed)
//     CfemReal l2Norm_U, h1SemiNorm_U;
//     IntegrationSchemeParams schemes_params{}; // Use default number of integration points
//     ScalarFunctionalEvaluator integrator(mesh, schemes_params);
//     integrator.L2Norm(u - analytical_u, &l2Norm_U, &h1SemiNorm_U);    
//     CFEM_INFO << "L2-norm error on u: " << l2Norm_U;
//     CFEM_INFO << "H1-semi norm error on u: " << h1SemiNorm_U;
        
//     // High-performance Binary Export
//     VTUExporter vtu_exporter(*mesh);
//     vtu_exporter.addExpression("u", u);
//     vtu_exporter.setFormat(VTKExportFormat::BINARY);
//     vtu_exporter.write("simulation_results");

//     return 0; 
// }

#include "cfemutils.h"
#include "Vector3D.h"
#include "Mesh.h"
#include "MeshConformityChecker.h"
#include "LagrangeSpace.h"
#include "DGSpace.h"
#include "FEField.h"
#include "ScalarFunctionalEvaluator.h"
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
static char help[] = "CFEM Example: Field interpolation on a quadratic mesh.\n\n";

int main(int argc, char **argv)
{
    // omp_set_num_threads(1);
    Environment::initialize(argc, argv, help);
    {
        #define CFEM_USE_DYNAMIC_SPACE_CACHE

        std::shared_ptr<Mesh> mesh;

        CfemReal x0 = 0;
        CfemReal y0 = 0;
        CfemReal z0 = 0;
        CfemReal r0 = 0.5;
        CfemReal h = 0.0625;

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

            // mesh = Mesh::disk(MeshOrder::Quadratic, 0.1, Vector3D{0,0,0}, 0.5);

            mesh = Mesh::block(
                MeshOrder::Linear,
                Mesh::BlockSubdivisions{150, 150, 150}, // Mesh::BlockSubdivisions{5*(2<<i), 5*(2<<i), 5*(2<<i)},
                Mesh::BlockBounds(
                    Vector3D{0., 0., 0.},
                    Vector3D{1., 1., 1.}
                ),
                Mesh::BlockElementType::Tetrahedron
            );

            // mesh = Mesh::loadFromGmsh("../tests/grenier/stanford_bunny/bunny.msh");
            // mesh = Mesh::loadFromGmsh("../tests/grenier/wheel/wheel.msh");

            // mesh = Mesh::sphere( MeshOrder::Linear,
            //                      h,
            //                      Vector3D{0,0,0},
            //                      r0
            //                 );

            
            // mesh = Mesh::surface_sphere( 
            //     MeshOrder::Quadratic,
            //     h,
            //     Vector3D{x0,y0,z0},
            //     r0
            // );
            // mesh = Mesh::sphere( MeshOrder::Quadratic,
            //                      0.25,
            //                      Vector3D{0,0,0},
            //                      0.5,
            //                     Mesh::SolidSphereElementType::Tetrahedron);

            // mesh = Mesh::hollow_cylinder(MeshOrder::Quadratic, 0.01, Vector3D{0, 0, 0}, 0.5, 0.75, 2, Vector3D{0, 1, 0});

        }

        // auto N = DGSpace::create(mesh, Order::P0, FieldType::Vector3D);
        auto Q = LagrangeSpace::create(mesh, Order::P1, FieldType::Scalar);
        auto V = LagrangeSpace::create(mesh, Order::P2, FieldType::Vector3D);

        // auto n = Trial(N, "n");
        auto p = Trial(Q, "p");
        auto u = Trial(V, "u");
        auto q = Test (Q, "q");
        auto v = Test (V, "v");

        // auto analytical_p = scalar(x*x + y + 1.0); 
        // auto analytical_u = vector( x*y + y*y +  x + y,
        //                             y*y  + x*x +  x*y + 1,              
        //                             0.0);

        // auto analytical_p = scalar(x + 10); 
        // auto analytical_u = vector( 0.0, 1, 0.0);
        
        // auto analytical_p = scalar(exp(x * y * z)); 
        // auto analytical_u = vector( x * x + y * y +  x + y,
        //                             y * y + z * z +  y + z,                    
        //                             x * x + z * z +  x + z);
        // auto normaleX = scalar ( (x - x0) / r0 );
        // auto normaleY = scalar ( (y - y0) / r0 );
        // auto normaleZ = scalar ( (z - z0) / r0 );
        // auto normale  = vector ( normaleX, normaleY, normaleZ );
        auto analytical_p = scalar(x * x + y * y + z * z + y + x + y + z + 1); // SymEngine expressions
        auto analytical_u = vector(sin(2 * CFEM_PI * x ) * sin(2 * CFEM_PI * z ),
                                   cos(2 * CFEM_PI * y ) * cos(2 * CFEM_PI * x ),
                                   cos(2 * CFEM_PI * y ) * cos(2 * CFEM_PI * z ));

        {
            PROFILE_SCOPE("FieldProjection");
            u->interpolate(analytical_u);
            p->interpolate(analytical_p);
            // n->interpolate(geometric_normal());
        }

        {
            PROFILE_SCOPE("L2Norm");
            IntegrationSchemeParams schemes_params{};
            // schemes_params.NumIntgPtTetElem = 33;
            ScalarFunctionalEvaluator integrator(mesh, schemes_params);
            CfemReal l2Norm_P, l2Norm_U, h1SemiNorm_P, h1SemiNorm_U;
            integrator.L2Norm(p - analytical_p, &l2Norm_P, &h1SemiNorm_P);
            integrator.L2Norm(u - analytical_u, &l2Norm_U, &h1SemiNorm_U);
            CFEM_INFO << "L2-norm error on P: " << l2Norm_P << " H1-Seminorm error on P: " << h1SemiNorm_P;
            CFEM_INFO << "L2-norm error on U: " << l2Norm_U << " H1-Seminorm error on U: " << h1SemiNorm_U;

            // auto Proj = tensor(1 -normaleX*normaleX,   -normaleX*normaleY,   -normaleX*normaleZ,
            //                   -normaleY*normaleX, 1 -normaleY*normaleY,   -normaleY*normaleZ,
            //                   -normaleZ*normaleX,   -normaleZ*normaleY, 1 -normaleZ*normaleZ);
            // CfemReal s = integrator.computeIntegral(scalar(1), SYS_BOUNDARY_ID);
            // CFEM_INFO << " Boundary surface = " << s << " Exact value: " << 4 * CFEM_PI * r0 * r0  << " Error " << std::abs( s - 4 * CFEM_PI * r0 * r0) ;

            // auto gradPanal = Proj*grad(analytical_p); //grad(analytical_p) - dot(grad(analytical_p), normale) * normale;
            // auto gradP = Proj*grad(p);  //grad(p) - dot(grad(p), normale) * normale;
            // CFEM_INFO << "surface error = " << sqrt( integrator.computeIntegral( mag_sqr( gradP - gradPanal ) ) );

            // CfemReal lseminorm_p = sqrt( integrator.computeIntegral( mag_sqr( grad_s(p) - grad_s(analytical_p) ) ) );
            // CfemReal lseminorm_u = sqrt( integrator.computeIntegral( mag_sqr( grad_s(u) - grad_s(analytical_u) ) ) );
            // CFEM_INFO << "surface gradient error for p  = " << lseminorm_p;
            // CFEM_INFO << "surface gradient error for u  = " << lseminorm_u;

            // auto gradUanal = dot(grad(analytical_u), Proj); //grad(analytical_p) - dot(grad(analytical_p), normale) * normale;
            // auto gradU = dot(grad(u), Proj);  //grad(p) - dot(grad(p), normale) * normale;
            // CFEM_INFO << "surface error = " << sqrt( integrator.computeIntegral( mag_sqr( gradU - gradUanal ) ) );
        }


        // h = 0.5  lseminorm_p = 0.0132748  lseminorm_u = 0.0408421
        // h = 0.25 lseminorm_p = 0.00349812 lseminorm_u = 0.00980175
        // h = 0.125 lseminorm_p = 0.000886251 lseminorm_u = 0.00241861
        // h = 0.0625 lseminorm_p = 0.000222298 lseminorm_u = 0.000602551
        // h = 0.03125 lseminorm_p = 5.56204e-05 lseminorm_u = 0.000150505

        {
            // PROFILE_SCOPE("Exportation");
            // VTKExporter vtu_exporter(*mesh);
            // vtu_exporter.addExpression("U", u);
            // vtu_exporter.addExpression("P", p);
            // // vtu_exporter.setExportedEntities({"Mag_external_surface", "Bande_de_roulement", "Caoutchouc"});
            // vtu_exporter.addExpression("N", geometric_normal());
            
            // vtu_exporter.setFormat(VTKExportFormat::BINARY);
            // vtu_exporter.exportCellQuality(true);
            // vtu_exporter.write("result_global");
            // vtu_exporter.write("result_boundary", SYS_BOUNDARY_ID);
            // // vtu_exporter.write("Gomme", "Gomme");
            // // vtu_exporter.write("Mag_external_surface", "Mag_external_surface");
            // // vtu_exporter.write("Bande_de_roulement", "Bande_de_roulement");
            // // vtu_exporter.write("Caoutchouc", "Caoutchouc");

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
//     auto analytical_u = vector(sin(2 * CFEM_PI * x), cos(2 * CFEM_PI * y), sin(CFEM_PI * x * y * z));

//     // 5. High-level Projection
//     u->interpolate(analytical_u);
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
// static char help[] = "CFEM Example: Field interpolation on a quadratic mesh.\n\n";

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

//         auto u_raw = scalar( sin(CFEM_PI * x) * exp(-x*x) ); 
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
//     auto raw = sym::pow(sin(CFEM_PI * x), 4.0)* exp(-x * x) * log(1.0 + x * x) * sqrt(1.0 + y * y);

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

