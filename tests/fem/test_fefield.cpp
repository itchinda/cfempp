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
#include "Mesh.h"
#include "FEField.h"
#include "MathOperators.h"
#include "FunctionIntegrator.h"
#include "LagrangeSpace.h"
#include "DGSpace.h"

using namespace cfem;
using namespace cfem::sym;
using namespace cfem::sym::symengine;

/**
 * @class FEFieldTest
 * @brief Test fixture for FEField validations.
 * Sets up a standard 3D block mesh and corresponding Finite Element Spaces
 * (Scalar P1 and Vector P2) to be used across multiple tests.
 */
class FEFieldTest : public ::testing::Test 
{
protected:
    std::shared_ptr<Mesh> m_mesh;

    std::shared_ptr<FESpace> m_V_scalar_p1;
    std::shared_ptr<FESpace> m_V_scalar_p2;
    std::shared_ptr<FESpace> m_V_vector_p1;
    std::shared_ptr<FESpace> m_V_vector_p2;
    
    std::shared_ptr<FESpace> m_V_dg_scalar_p0;
    std::shared_ptr<FESpace> m_V_dg_scalar_p1;
    std::shared_ptr<FESpace> m_V_dg_scalar_p2;
    std::shared_ptr<FESpace> m_V_dg_vector_p0;
    std::shared_ptr<FESpace> m_V_dg_vector_p1;
    std::shared_ptr<FESpace> m_V_dg_vector_p2;

    void SetUp() override 
    {
        m_mesh = Mesh::block(MeshOrder::Quadratic, {4, 4, 4}, 
                             Mesh::BlockBounds({0.0, 0.0, 0.0}, {1.0, 1.0, 1.0}), 
                             Mesh::BlockElementType::Tetrahedron);

        // CG Spaces
        m_V_scalar_p1 = LagrangeSpace::create(m_mesh, Order::P1, FieldType::Scalar);
        m_V_scalar_p2 = LagrangeSpace::create(m_mesh, Order::P2, FieldType::Scalar);
        m_V_vector_p1 = LagrangeSpace::create(m_mesh, Order::P1, FieldType::Vector3D);
        m_V_vector_p2 = LagrangeSpace::create(m_mesh, Order::P2, FieldType::Vector3D);
        
        // DG Spaces
        m_V_dg_scalar_p0 = DGSpace::create(m_mesh, Order::P0, FieldType::Scalar);
        m_V_dg_scalar_p1 = DGSpace::create(m_mesh, Order::P1, FieldType::Scalar);
        m_V_dg_scalar_p2 = DGSpace::create(m_mesh, Order::P2, FieldType::Scalar);

        m_V_dg_vector_p0 = DGSpace::create(m_mesh, Order::P0, FieldType::Vector3D);
        m_V_dg_vector_p1 = DGSpace::create(m_mesh, Order::P1, FieldType::Vector3D);
        m_V_dg_vector_p2 = DGSpace::create(m_mesh, Order::P2, FieldType::Vector3D);
    }
};

// ========================================================================
// --- 1. Basic Initialization & API Tests
// ========================================================================

TEST_F(FEFieldTest, InitializationAndProperties) 
{
    auto tempField = FEField::create("Temperature", m_V_scalar_p1);
    
    // Check naming and sym registry
    EXPECT_EQ(tempField->getName(), "Temperature");
    EXPECT_EQ(tempField->toString(), "FEField::Temperature");
    
    // Check dimensional properties
    EXPECT_EQ(tempField->getDimension(), 1);
    
    // Check expression tree properties
    EXPECT_TRUE(tempField->dependsOnFEField());
    EXPECT_TRUE(tempField->dependsOnShape());
    EXPECT_EQ(tempField->getType(), sym::ExpressionType::FEField);
    
    // Check PETSc Vector allocation
    const PetscVector& data = tempField->getData();
    EXPECT_EQ(data.globalSize(), m_V_scalar_p1->getGlobalNumDofs());
}

// ========================================================================
// --- 2. Data Management Tests
// ========================================================================

TEST_F(FEFieldTest, LocalBufferSynchronization) 
{
    auto pressure = FEField::create("Pressure", m_V_scalar_p2);
    Vec petscVec = pressure->getData().getRaw();

    // Directly alter the PETSc backend vector
    PetscScalar testValue = 42.5;
    VecSet(petscVec, testValue);

    // Trigger local buffer synchronization
    pressure->updateLocalBuffer();

    // Verify that getAtNode accesses the correctly synchronized local buffer
    EXPECT_DOUBLE_EQ(pressure->getAtNode(0, 0), testValue);
    EXPECT_DOUBLE_EQ(pressure->getAtNode(m_mesh->getNumVertices() - 1, 0), testValue);
}

TEST_F(FEFieldTest, ReallocationClearsData) 
{
    auto field = FEField::create("Temp", m_V_scalar_p1);
    VecSet(field->getData().getRaw(), 10.0);
    field->updateLocalBuffer();
    
    field->reallocate();
    field->updateLocalBuffer();

    EXPECT_DOUBLE_EQ(field->getAtNode(0, 0), 0.0);
}

// ========================================================================
// --- 1. TESTS DE COHÉRENCE DG (NOMBRE DE DOFS)
// ========================================================================

TEST_F(FEFieldTest, DGSpaceDofConsistency) 
{
    CfemInt nCells = m_mesh->getNumCells();

    // P0 DG: 1 DOF par cellule (Scalaire), 3 DOFs par cellule (Vector3D)
    EXPECT_EQ(m_V_dg_scalar_p0->getGlobalNumDofs(), nCells * 1);
    EXPECT_EQ(m_V_dg_vector_p0->getGlobalNumDofs(), nCells * 3);

    // P1 DG (Tetrahedron): 4 Nodes par cellule
    EXPECT_EQ(m_V_dg_scalar_p1->getGlobalNumDofs(), nCells * 4);
    EXPECT_EQ(m_V_dg_vector_p1->getGlobalNumDofs(), nCells * 4 * 3);

    // P2 DG (Tetrahedron): 10 Nodes par cellule
    EXPECT_EQ(m_V_dg_scalar_p2->getGlobalNumDofs(), nCells * 10);
    EXPECT_EQ(m_V_dg_vector_p2->getGlobalNumDofs(), nCells * 10 * 3);
}

TEST_F(FEFieldTest, CG_Scalar_P1_Exactness) 
{
    auto f = FEField::create("f", m_V_scalar_p1);
    auto analytical = scalar(2.0 * x + y - 0.5 * z);
    f->project(analytical);

    CfemReal l2, h1;
    FunctionIntegrator(m_mesh).L2Norm(f - analytical, &l2, &h1);
    EXPECT_NEAR(l2, 0.0, 1e-12);
    EXPECT_NEAR(h1, 0.0, 1e-11);
}

TEST_F(FEFieldTest, CG_Vector_P1_Exactness) 
{
    auto velocity = FEField::create("Velocity", m_V_vector_p1);
    auto analytical = vector(1.5 * x, -2.0 * y, 0.5 * z);
    velocity->project(analytical);

    CfemReal l2Norm, h1SemiNorm;
    FunctionIntegrator integrator(m_mesh, IntegrationSchemeParams{});
    integrator.L2Norm(velocity - analytical, &l2Norm, &h1SemiNorm);

    EXPECT_NEAR(l2Norm, 0.0, 1e-12) << "P1 GC Vector L2 projection failed.";
    EXPECT_NEAR(h1SemiNorm, 0.0, 1e-11) << "P1 CG Vector H1 projection failed.";
}

TEST_F(FEFieldTest, CG_Scalar_P2_Exactness) 
{
    auto f = FEField::create("f", m_V_scalar_p2);
    auto analytical = scalar(x*x + y*z + 3.0);
    f->project(analytical);

    CfemReal l2, h1;
    FunctionIntegrator(m_mesh).L2Norm(f - analytical, &l2, &h1);
    EXPECT_NEAR(l2, 0.0, 1e-12);
    EXPECT_NEAR(h1, 0.0, 1e-11);
}

TEST_F(FEFieldTest, CG_Vector_P2_Lambda_OpenMP) 
{
    auto u = FEField::create("u", m_V_vector_p2);
    // Test avec une lambda complexe pour stress-tester OpenMP et le mapping P2
    u->project([](const Vector3D& p) -> std::array<CfemReal, 3> {
        return { p.x*p.x, p.x*p.y, p.x*p.z };
    });

    auto analytical = vector(x*x, y*x, x*z);
    CfemReal l2;
    FunctionIntegrator(m_mesh).L2Norm(u - analytical, &l2, nullptr);
    EXPECT_NEAR(l2, 0.0, 1e-11);
    EXPECT_NEAR(l2, 0.0, 1e-11);
}

TEST_F(FEFieldTest, DG_Scalar_P0_Constant) 
{
    auto f = FEField::create("f", m_V_dg_scalar_p0);
    auto analytical = 3.14159;
    f->project(analytical);

    CfemReal l2;
    FunctionIntegrator(m_mesh).L2Norm(f - analytical, &l2, nullptr);
    EXPECT_NEAR(l2, 0.0, 1e-12);
}

TEST_F(FEFieldTest, DG_Vector_P0_Constant) 
{
    auto u = FEField::create("u", m_V_dg_vector_p0);
    auto analytical = vector(1.0, 2.0, 3.0);
    u->project(analytical);

    CfemReal l2;
    FunctionIntegrator(m_mesh).L2Norm(u - analytical, &l2, nullptr);
    EXPECT_NEAR(l2, 0.0, 1e-12);
}

TEST_F(FEFieldTest, DG_Scalar_P1_Exactness) 
{
    auto f = FEField::create("f", m_V_dg_scalar_p1);
    auto analytical = scalar(0.5 * x + y);
    f->project(analytical);

    CfemReal l2, h1;
    FunctionIntegrator(m_mesh).L2Norm(f - analytical, &l2, &h1);
    EXPECT_NEAR(l2, 0.0, 1e-12);
    EXPECT_NEAR(h1, 0.0, 1e-11);
}

TEST_F(FEFieldTest, DG_Vector_P1_Lambda) 
{
    auto u = FEField::create("u", m_V_dg_vector_p1);
    u->project([](const Vector3D& p) -> std::array<CfemReal, 3> {
        return { p.x + 1.0, p.y - 2.0, p.z };
    });

    auto analytical = vector(x + 1.0, y - 2.0, z);
    CfemReal l2;
    FunctionIntegrator(m_mesh).L2Norm(u - analytical, &l2, nullptr);
    EXPECT_NEAR(l2, 0.0, 1e-12);
}

TEST_F(FEFieldTest, DG_Scalar_P2_Exactness) 
{
    auto f = FEField::create("f", m_V_dg_scalar_p2);
    auto analytical = scalar(x*x + y*y + z*z);
    f->project(analytical);

    CfemReal l2, h1;
    FunctionIntegrator(m_mesh).L2Norm(f - analytical, &l2, &h1);
    EXPECT_NEAR(l2, 0.0, 1e-12);
}

TEST_F(FEFieldTest, DG_Vector_P2_Exactness) 
{
    auto u = FEField::create("u", m_V_dg_vector_p2);
    auto analytical = vector(x*x, y*y, z*z);
    u->project(analytical);

    CfemReal l2, h1;
    FunctionIntegrator(m_mesh).L2Norm(u - analytical, &l2, &h1);
    EXPECT_NEAR(l2, 0.0, 1e-12);
}

// ========================================================================
// --- 4. TESTS AVANCÉS (GRADIENTS, MPI SYNC, REALLOCATION)
// ========================================================================

TEST_F(FEFieldTest, GradientEvaluation_P2_CG) 
{
    auto f = FEField::create("f", m_V_scalar_p2);
    auto analytical = scalar(x*x + 2.0*y*y);
    auto analytical_grad = vector(2.0*x, 4.0*y, 0.0);
    
    f->project(analytical);

    // L'intégrateur H1 semi-norm calcule ||grad(f) - grad(analytique)||
    CfemReal h1 = sqrt( FunctionIntegrator(m_mesh).computeIntegral(mag_sqr( grad(f - analytical) ) ) );
    
    EXPECT_NEAR(h1, 0.0, 1e-11) << "Gradient evaluation failed for P2 CG field.";
}

TEST_F(FEFieldTest, MPI_LocalBufferUpdate) 
{
    auto f = FEField::create("f", m_V_scalar_p1);
    Vec v = f->getData().getRaw();
    VecSet(v, 1.234);
    
    // Test de synchronisation du buffer local utilisé par evaluate()
    f->updateLocalBuffer();
    
    // On vérifie une valeur arbitraire dans le buffer local
    EXPECT_DOUBLE_EQ(f->getAtNode(0, 0), 1.234);
}

TEST_F(FEFieldTest, ReallocatePreservesMetadata) 
{
    auto u = FEField::create("u", m_V_vector_p2);
    std::string originalName = u->getName();
    auto originalSpace = u->getSpace();
    
    u->reallocate(); // Simule un changement de maillage ou de distribution MPI
    
    EXPECT_EQ(u->getName(), originalName);
    EXPECT_EQ(u->getSpace(), originalSpace);
    EXPECT_EQ(u->getData().globalSize(), originalSpace->getGlobalNumDofs());
}

TEST_F(FEFieldTest, DGScalarP0Projection) 
{
    // P0 is piecewise constant. It can only exactly represent a constant field.
    auto density = FEField::create("Density", m_V_dg_scalar_p0);
    auto analytical = 1050.0; // e.g., density of a fluid
    
    density->project(analytical);

    CfemReal l2Norm = 1.0;
    FunctionIntegrator integrator(m_mesh, IntegrationSchemeParams{});
    integrator.L2Norm(density - analytical, &l2Norm, nullptr);

    EXPECT_NEAR(l2Norm, 0.0, 1e-12) << "P0 DG L2 projection error is too high. Check cell-center mapping.";
}