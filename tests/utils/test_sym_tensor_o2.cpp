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
#include "TensorO2.h"

using namespace cfem;

const double TOL = 1e-12;

TEST(SymTensorO2Test, InitializationAndIdentity) {
    auto I = SymTensorO2::Identity();
    EXPECT_DOUBLE_EQ(I[0], 1.0);
    EXPECT_DOUBLE_EQ(I[1], 1.0);
    EXPECT_DOUBLE_EQ(I[2], 1.0);
    EXPECT_DOUBLE_EQ(I[3], 0.0);
    
    EXPECT_DOUBLE_EQ(I.trace(), 3.0);
    EXPECT_DOUBLE_EQ(I.mean(), 1.0);
}

TEST(SymTensorO2Test, Arithmetic) {
    SymTensorO2 A = {1, 2, 3, 4, 5, 6};
    SymTensorO2 B = {6, 5, 4, 3, 2, 1};
    
    auto C = A + B;
    for(int i=0; i<6; ++i) EXPECT_DOUBLE_EQ(C[i], 7.0);
    
    auto D = A * 2.0;
    EXPECT_DOUBLE_EQ(D[5], 12.0);

    auto E = 2.0 * A;
    EXPECT_DOUBLE_EQ(E[5], 12.0);
    
    auto F = -A;
    EXPECT_DOUBLE_EQ(F[0], -1.0);
}

TEST(SymTensorO2Test, TractionVector) {
    // Pure shear stress sigma_xy = 50
    SymTensorO2 sigma = {0, 0, 0, 0, 0, 50.0};
    // Normal in X direction
    Vector3D n = {1.0, 0.0, 0.0};
    
    // Traction t = sigma * n should be {0, 50, 0}
    Vector3D t = sigma * n;
    EXPECT_NEAR(t.x, 0.0, TOL);
    EXPECT_NEAR(t.y, 50.0, TOL);
    EXPECT_NEAR(t.z, 0.0, TOL);
}

TEST(SymTensorO2Test, InvariantsAndPhysics) {
    // Example: Pure hydrostatic pressure
    SymTensorO2 pressure = {100, 100, 100, 0, 0, 0};
    EXPECT_DOUBLE_EQ(pressure.mean(), 100.0);
    EXPECT_NEAR(pressure.vonMises(), 0.0, TOL);
    
    // Example: Pure uniaxial tension in Z
    SymTensorO2 tension = {0, 0, 100, 0, 0, 0};
    // For uniaxial tension, Von Mises stress == tension value
    EXPECT_NEAR(tension.vonMises(), 100.0, TOL);
    
    // Check J2 for uniaxial tension: J2 = (1/3) * sigma^2 = (1/3)*10000 = 3333.333
    EXPECT_NEAR(tension.J2(), 10000.0/3.0, TOL);
}

TEST(SymTensorO2Test, DeviatoricSplit) {
    SymTensorO2 A = {10, 20, 30, 1, 2, 3};
    auto s = A.deviator();
    
    // The trace of any deviatoric tensor must be zero
    EXPECT_NEAR(s.trace(), 0.0, TOL);
    
    // Off-diagonals should remain unchanged
    EXPECT_DOUBLE_EQ(s[3], 1.0);
    EXPECT_DOUBLE_EQ(s[5], 3.0);
}