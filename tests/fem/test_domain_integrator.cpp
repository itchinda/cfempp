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

#include <gtest/gtest.h>

#include "cfemutils.h"
#include "FEField.h"
#include "LagrangeSpace.h"
#include "MathOperators.h"
#include "Mesh.h"
#include "ScalarFunctionalEvaluator.h"

using namespace cfem;
using namespace cfem::sym;
using namespace cfem::sym::symengine;

// Fixture pour initialiser l'environnement et créer un maillage commun
class ScalarFunctionalEvaluatorTest : public ::testing::Test {
protected:
    void SetUp() override {
        // On crée un cube [-1, 1]^3 (Volume = 8.0, Surface d'une face = 4.0)
        mesh = Mesh::block(
            MeshOrder::Linear,
            Mesh::BlockSubdivisions{4, 4, 4}, 
            Mesh::BlockBounds(Vector3D{-1., -1., -1.}, Vector3D{1., 1., 1.}),
            Mesh::BlockElementType::Tetrahedron
        );
        
        IntegrationSchemeParams params{};
        integrator = std::make_unique<ScalarFunctionalEvaluator>(mesh, params);
    }

    std::shared_ptr<Mesh> mesh;
    std::unique_ptr<ScalarFunctionalEvaluator> integrator;
    const CfemReal TOL = 1e-12;
};

// 1. Validation de l'intégration de volume (Fonction constante)
TEST_F(ScalarFunctionalEvaluatorTest, VolumeIntegralConstant) {
    auto one = scalar(1.0);
    CfemReal volume = integrator->computeIntegral(one);
    // Cube de côté 2 -> Volume = 2^3 = 8
    EXPECT_NEAR(volume, 8.0, TOL);
}

// 2. Validation de l'intégration de surface (Fonction constante sur frontière)
TEST_F(ScalarFunctionalEvaluatorTest, SurfaceIntegralConstant) {
    auto one = scalar(1.0);
    // Face "xmin" est à x = -1, de y=-1..1 et z=-1..1 -> Aire = 4.0
    CfemReal area = integrator->computeIntegral(one, "xmin");
    EXPECT_NEAR(area, 4.0, TOL);
}

// 3. Validation de l'interpolation et de la norme L2 (Cas Polynomial P2)
TEST_F(ScalarFunctionalEvaluatorTest, L2NormP2Projection) {
    auto Q = LagrangeSpace::create(mesh, Order::P2, FieldType::Scalar);
    auto p = Trial(Q, "p");

    // Une fonction quadratique doit être représentée EXACTEMENT en P2
    auto analytical = scalar(x * x + y * y + z * z);
    
    p->interpolate(analytical);

    CfemReal l2Err = 1.0, h1Err = 1.0;
    integrator->L2Norm(p - analytical, &l2Err, &h1Err);

    // L'erreur doit être proche de la précision machine
    EXPECT_NEAR(l2Err, 0.0, TOL);
    EXPECT_NEAR(h1Err, 0.0, TOL);
}

// // 4. Test d'un champ vectoriel (Divergence / Gradient)
// TEST_F(ScalarFunctionalEvaluatorTest, VectorFieldDivergence) {
//     auto V = LagrangeSpace::create(mesh, Order::P2, FieldType::Vector3D);
//     auto u = Trial(V, "u");

//     // u = [x, y, z] -> div(u) = 1 + 1 + 1 = 3
//     auto analytical_u = vector(x, y, z);
//     u->interpolate(analytical_u);

//     // Intégrale de div(u) sur le volume
//     // Théorème de la divergence : int(div(u)) = int(u . n) sur la surface
//     // Ici int(3) dv = 3 * Volume = 3 * 8 = 24
//     auto expr_div = div(u);
//     CfemReal integral_div = integrator->computeIntegral(expr_div);
    
//     EXPECT_NEAR(integral_div, 24.0, 1e-11);
// }

// 5. Test de symétrie / Coordonnées
TEST_F(ScalarFunctionalEvaluatorTest, IntegralLinearFunction) {
    // L'intégrale de 'x' sur [-1, 1] est 0 par symétrie
    auto expr_x = scalar(x);
    CfemReal val = integrator->computeIntegral(expr_x);
    EXPECT_NEAR(val, 0.0, TOL);
    
    // L'intégrale de 'x+1' sur le volume 8 est 8.0
    CfemReal val2 = integrator->computeIntegral(scalar(x + 1.0));
    EXPECT_NEAR(val2, 8.0, TOL);
}

// 6. Test de l'erreur sur une face spécifique (Gradient normal)
TEST_F(ScalarFunctionalEvaluatorTest, SurfaceGradientError) {
    auto Q = LagrangeSpace::create(mesh, Order::P2, FieldType::Scalar);
    auto p = Trial(Q, "p2");
    
    auto analytical = scalar(x * x + y * z + y * y - z * z); 
    p->interpolate(analytical);
    
    // Sur la face x=1, grad(x * x + y * z - z * z) = [2x, z, 2z] = [2, 2y + z, y - 2z]
    // Intégrale sur la face = 88/3
    CfemReal surface_grad_int = integrator->computeIntegral(mag_sqr(grad(p)), "xmax");
    
    EXPECT_NEAR(surface_grad_int, 88.0/3.0, TOL);
}

class SurfaceManifoldTest : public ::testing::Test {
protected:
    void SetUp() override {
        VariableRegistry::clear();
    }

    // Paramètres de la sphère
    const CfemReal x0 = 0.0, y0 = 0.0, z0 = 0.0;
    const CfemReal r0 = 0.5;
    const CfemReal h = 0.03125; // Maillage assez fin pour la convergence
};

/**
 * @brief Vérifie que le flux calculé sur une sphère 2D (Manifold)
 * correspond à la valeur analytique PI/2.
 */
TEST_F(SurfaceManifoldTest, ComputeFluxOnSphere) {
    auto mesh = Mesh::surface_sphere(MeshOrder::Quadratic, h, Vector3D{x0, y0, z0}, r0);
    ASSERT_NE(mesh, nullptr);

    auto V = LagrangeSpace::create(mesh, Order::P2, FieldType::Vector3D);
    auto u = Trial(V, "u");

    // u = [ x² + y² + x + y,  y² + z² + y + z,  x² + z² + x + z ]
    auto analytical_u = vector( x * x + y * y + x + y,
                                y * y + z * z + y + z,
                                x * x + z * z + x + z);

    u->interpolate(analytical_u);

    IntegrationSchemeParams params{};
    ScalarFunctionalEvaluator integrator(mesh, params);

    CfemReal computedFlux = integrator.computeFlux(u, UNIVERSE_ID);

    // (Valeur théorique : PI / 2)
    const CfemReal expectedFlux = M_PI / 2.0;
    EXPECT_NEAR(computedFlux, expectedFlux, 1e-5);
}

/**
 * @brief Vérifie la convergence du gradient de surface (grad_s).
 */
TEST_F(SurfaceManifoldTest, SurfaceGradientError) {
    auto mesh = Mesh::surface_sphere(MeshOrder::Quadratic, h, Vector3D{x0, y0, z0}, r0);
    
    auto Q = LagrangeSpace::create(mesh, Order::P2, FieldType::Scalar);
    auto p = Trial(Q, "p");


    // Champ scalaire analytique : exp(x*y*z)
    auto analytical_p = scalar(exp(x * y * z));
    p->interpolate(analytical_p);

    IntegrationSchemeParams params{};
    ScalarFunctionalEvaluator integrator(mesh, params);

    CfemReal error_grad_s = sqrt(integrator.computeIntegral(mag_sqr(grad_s(p) - grad_s(analytical_p))));
    EXPECT_LT(error_grad_s, 1e-4);
}