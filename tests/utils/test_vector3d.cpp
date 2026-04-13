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
#include "Vector3D.h"

using namespace cfem;

/**
 * @brief Unit tests for the Vector3D struct.
 */

// Tolerance for floating point comparisons
const double TOL = cfem::ZERO_TOLERANCE;

TEST(Vector3DTest, BasicArithmetic) {
    Vector3D v1{1.0, 2.0, 3.0};
    Vector3D v2{4.0, 5.0, 6.0};

    // Test Addition
    Vector3D v_sum = v1 + v2;
    EXPECT_NEAR(v_sum.x, 5.0, TOL);
    EXPECT_NEAR(v_sum.y, 7.0, TOL);
    EXPECT_NEAR(v_sum.z, 9.0, TOL);

    // Test Subtraction
    Vector3D v_diff = v2 - v1;
    EXPECT_NEAR(v_diff.x, 3.0, TOL);
    EXPECT_NEAR(v_diff.y, 3.0, TOL);
    EXPECT_NEAR(v_diff.z, 3.0, TOL);

    // Test Scalar multiplication
    Vector3D v_mult = v1 * 2.0;
    EXPECT_NEAR(v_mult.x, 2.0, TOL);
    EXPECT_NEAR(v_mult.y, 4.0, TOL);
    EXPECT_NEAR(v_mult.z, 6.0, TOL);
    
    // Test Left-side scalar multiplication
    Vector3D v_lmult = 3.0 * v1;
    EXPECT_NEAR(v_lmult.x, 3.0, TOL);
}

TEST(Vector3DTest, VectorCalculus) {
    Vector3D v1{1.0, 0.0, 0.0};
    Vector3D v2{0.0, 1.0, 0.0};

    // Test Dot Product (Orthogonal vectors)
    EXPECT_NEAR(v1.dot(v2), 0.0, TOL);
    
    // Test Dot Product (Parallel)
    EXPECT_NEAR(v1.dot(v1), 1.0, TOL);

    // Test Cross Product (Standard Basis)
    Vector3D v_cross = v1.cross(v2);
    EXPECT_NEAR(v_cross.x, 0.0, TOL);
    EXPECT_NEAR(v_cross.y, 0.0, TOL);
    EXPECT_NEAR(v_cross.z, 1.0, TOL); // x cross y = z

    // Test Norm
    Vector3D v3{3.0, 4.0, 0.0}; // 3-4-5 triangle
    EXPECT_NEAR(v3.norm(), 5.0, TOL);
    EXPECT_NEAR(v3.normSquare(), 25.0, TOL);
}

TEST(Vector3DTest, Normalization) {
    Vector3D v{10.0, 0.0, 0.0};
    Vector3D v_unit = v.normalized();
    
    EXPECT_NEAR(v_unit.norm(), 1.0, TOL);
    EXPECT_NEAR(v_unit.x, 1.0, TOL);

    // Test near-null vector normalization
    Vector3D v_null{1e-18, 1e-18, 1e-18};
    Vector3D v_null_unit = v_null.normalized();
    EXPECT_NEAR(v_null_unit.norm(), 0.0, TOL);
}

TEST(Vector3DTest, Accessors) {
    Vector3D v{1.1, 2.2, 3.3};
    
    // Read
    EXPECT_EQ(v[0], 1.1);
    
    // Write
    v[1] = 5.5;
    EXPECT_EQ(v.y, 5.5);
}

TEST(Vector3DTest, Printing) {
    Vector3D v{1.0, 2.0, 3.0};
    std::stringstream ss;
    ss << v;
    EXPECT_EQ(ss.str(), "(1, 2, 3)");
}