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
#include "TensorO2.h"

using namespace cfem;

const CfemReal TOL = 1e-12;

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


class SkewSymTensorO2Test : public ::testing::Test {
protected:
    void SetUp() override {
    }
};

// ========================================================================
// Accessors and Skew-symmetry validation
// ========================================================================
TEST_F(SkewSymTensorO2Test, AccessAndSignReversal) {
    SkewSymTensorO2 W;
    W.data[idSSTO2::xy] =  1.5;
    W.data[idSSTO2::xz] = -2.0;
    W.data[idSSTO2::yz] =  3.1;

    // Strictly null diag
    EXPECT_DOUBLE_EQ(W(0, 0), 0.0);
    EXPECT_DOUBLE_EQ(W(1, 1), 0.0);
    EXPECT_DOUBLE_EQ(W(2, 2), 0.0);

    // Upper components
    EXPECT_DOUBLE_EQ(W(0, 1),  1.5);
    EXPECT_DOUBLE_EQ(W(0, 2), -2.0);
    EXPECT_DOUBLE_EQ(W(1, 2),  3.1);

    // Lower components (Sign inversion)
    EXPECT_DOUBLE_EQ(W(1, 0), -1.5);
    EXPECT_DOUBLE_EQ(W(2, 0),  2.0);
    EXPECT_DOUBLE_EQ(W(2, 1), -3.1);
}

// ========================================================================
// Axial vector (w) and conversion
// ========================================================================
TEST_F(SkewSymTensorO2Test, AxialVectorConversion) {
    // Axial vector w = [w_x, w_y, w_z]
    FixedVector<CfemReal, 3> w_in;
    w_in[0] =  4.0;
    w_in[1] = -5.0;
    w_in[2] =  6.0;

    // Conversion Vector -> Tensor
    SkewSymTensorO2 W = SkewSymTensorO2::fromAxialVector(w_in);

    // Conversion check W * v = w x v
    // W_xy = -w_z, W_xz = w_y, W_yz = -w_x
    EXPECT_DOUBLE_EQ(W.data[idSSTO2::xy], -6.0);
    EXPECT_DOUBLE_EQ(W.data[idSSTO2::xz], -5.0);
    EXPECT_DOUBLE_EQ(W.data[idSSTO2::yz], -4.0);

    // Inverse conversion Tensor -> Vector
    FixedVector<CfemReal, 3> w_out = W.axialVector();

    EXPECT_DOUBLE_EQ(w_out[0], w_in[0]);
    EXPECT_DOUBLE_EQ(w_out[1], w_in[1]);
    EXPECT_DOUBLE_EQ(w_out[2], w_in[2]);
}

// ========================================================================
// Transpose
// ========================================================================
TEST_F(SkewSymTensorO2Test, Transpose) {
    SkewSymTensorO2 W;
    W.data[idSSTO2::xy] = 1.0;
    W.data[idSSTO2::xz] = 2.0;
    W.data[idSSTO2::yz] = 3.0;

    SkewSymTensorO2 WT = W.transpose();

    // W^T = -W (skew-symmetric)
    EXPECT_DOUBLE_EQ(WT.data[idSSTO2::xy], -1.0);
    EXPECT_DOUBLE_EQ(WT.data[idSSTO2::xz], -2.0);
    EXPECT_DOUBLE_EQ(WT.data[idSSTO2::yz], -3.0);
}

// ========================================================================
// Math properties (Determinant et and eigen values)
// ========================================================================
TEST_F(SkewSymTensorO2Test, MathProperties) {
    FixedVector<CfemReal, 3> w;
    w[0] = 3.0; w[1] = 4.0; w[2] = 0.0; // Norme = 5.0
    SkewSymTensorO2 W = SkewSymTensorO2::fromAxialVector(w);

    EXPECT_DOUBLE_EQ(det(W), 0.0);

    // The imaginary eigen values magnitude must be equal to the magnitude(norm) of the axial vector
    EXPECT_DOUBLE_EQ(W.imaginaryEigenvalueMagnitude(), 5.0);
}

// ========================================================================
//  Exponential (Rodrigues' formula)
// ========================================================================
TEST_F(SkewSymTensorO2Test, ExponentialRotationZ) {
    // 90 degree (PI/2) rotation around the Z-axis
    FixedVector<CfemReal, 3> w;
    w[0] = 0.0;
    w[1] = 0.0;
    w[2] = CFEM_PI / 2.0; 

    SkewSymTensorO2 W = SkewSymTensorO2::fromAxialVector(w);
    TensorO2 R = W.exp();

    /* Rotation matrix 90 degree around Z-axis :
       [  0  -1   0 ]
       [  1   0   0 ]
       [  0   0   1 ]
    */
    EXPECT_NEAR(R(0, 0),  0.0, TOL);
    EXPECT_NEAR(R(0, 1), -1.0, TOL);
    EXPECT_NEAR(R(0, 2),  0.0, TOL);

    EXPECT_NEAR(R(1, 0),  1.0, TOL);
    EXPECT_NEAR(R(1, 1),  0.0, TOL);
    EXPECT_NEAR(R(1, 2),  0.0, TOL);

    EXPECT_NEAR(R(2, 0),  0.0, TOL);
    EXPECT_NEAR(R(2, 1),  0.0, TOL);
    EXPECT_NEAR(R(2, 2),  1.0, TOL);
}

TEST_F(SkewSymTensorO2Test, ExponentialSmallAngleLimit) {
    // Taylor expansion test for small angle (theta < 1e-14)
    FixedVector<CfemReal, 3> w;
    w[0] = 1e-15; // Too small angles
    w[1] = 0.0;
    w[2] = 0.0;

    SkewSymTensorO2 W = SkewSymTensorO2::fromAxialVector(w);
    TensorO2 R = W.exp();

    // R must be equal to I + W
    EXPECT_NEAR(R(0, 0), 1.0, TOL);
    EXPECT_NEAR(R(1, 1), 1.0, TOL);
    EXPECT_NEAR(R(2, 2), 1.0, TOL);
    
    // Perturbation
    EXPECT_NEAR(R(2, 1),  1e-15, TOL);
    EXPECT_NEAR(R(1, 2), -1e-15, TOL);
}