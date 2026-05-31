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

#include <array>
#include <cmath>
#include <gtest/gtest.h>

#include "Tensor3D.h"

using namespace cfem;

// Standard Matrix substraction to be used as the reference result
Mat3x3 ReferenceSub(const Mat3x3& A, const Mat3x3& B) {
    Mat3x3 C = {{{0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}}};
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
    return C;
}

template<typename Mat>
void ExpectTensorEqMat(const Mat& res, const Mat3x3& expected) {
    constexpr CfemReal tol = 1e-12;
    EXPECT_NEAR(res(0, 0), expected[0][0], tol);
    EXPECT_NEAR(res(0, 1), expected[0][1], tol);
    EXPECT_NEAR(res(0, 2), expected[0][2], tol);
    EXPECT_NEAR(res(1, 0), expected[1][0], tol);
    EXPECT_NEAR(res(1, 1), expected[1][1], tol);
    EXPECT_NEAR(res(1, 2), expected[1][2], tol);
    EXPECT_NEAR(res(2, 0), expected[2][0], tol);
    EXPECT_NEAR(res(2, 1), expected[2][1], tol);
    EXPECT_NEAR(res(2, 2), expected[2][2], tol);
}

// ==============================================================================
// FIXTURE
// ==============================================================================

class TensorOperatorsTest : public ::testing::Test {
protected:
    Tensor3D T;
    SymTensor3D S;
    AxiSymTensor3D A;
    SkewSymTensor3D W;

    void SetUp() override {

        T.xx() = 1.0; T.xy() = 2.0; T.xz() = 3.0;
        T.yx() = 4.0; T.yy() = 5.0; T.yz() = 6.0;
        T.zx() = 7.0; T.zy() = 8.0; T.zz() = 9.0;

        S.xx() = 1.1; S.yy() = 2.2; S.zz() = 3.3;
        S.xy() = 4.4; S.xz() = 5.5; S.yz() = 6.6;

        A.xx() = 7.7; A.yy() = 8.8; A.zz() = 9.9;
        A.xy() = 1.2;

        W.xy() = 2.3; W.xz() = 3.4; W.yz() = 4.5;
    }
};

TEST_F(TensorOperatorsTest, Tensor3D_sub_SymTensor3D) {
    auto res = T - S;
    Mat3x3 expected = ReferenceSub(T.toMat3x3(), S.toMat3x3());
    ExpectTensorEqMat(res, expected);
}

TEST_F(TensorOperatorsTest, Tensor3D_sub_AxiSymTensor3D) {
    auto res = T - A;
    Mat3x3 expected = ReferenceSub(T.toMat3x3(), A.toMat3x3());
    ExpectTensorEqMat(res, expected);
}

TEST_F(TensorOperatorsTest, Tensor3D_sub_SkewSymTensor3D) {
    auto res = T - W;
    Mat3x3 expected = ReferenceSub(T.toMat3x3(), W.toMat3x3());
    ExpectTensorEqMat(res, expected);
}

TEST_F(TensorOperatorsTest, AxiSymTensor3D_sub_Tensor3D) {
    auto res = A - T;
    Mat3x3 expected = ReferenceSub(A.toMat3x3(), T.toMat3x3());
    ExpectTensorEqMat(res, expected);
}

TEST_F(TensorOperatorsTest, AxiSymTensor3D_sub_SymTensor3D) {
    auto res = A - S;
    Mat3x3 expected = ReferenceSub(A.toMat3x3(), S.toMat3x3());
    ExpectTensorEqMat(res, expected);
}

TEST_F(TensorOperatorsTest, AxiSymTensor3D_sub_SkewSymTensor3D) {
    auto res = A - W;
    Mat3x3 expected = ReferenceSub(A.toMat3x3(), W.toMat3x3());
    ExpectTensorEqMat(res, expected);
}

TEST_F(TensorOperatorsTest, SymTensor3D_sub_Tensor3D) {
    auto res = S - T;
    Mat3x3 expected = ReferenceSub(S.toMat3x3(), T.toMat3x3());
    ExpectTensorEqMat(res, expected);
}

TEST_F(TensorOperatorsTest, SymTensor3D_sub_SkewSymTensor3D) {
    auto res = S - W;
    Mat3x3 expected = ReferenceSub(S.toMat3x3(), W.toMat3x3());
    ExpectTensorEqMat(res, expected);
}

TEST_F(TensorOperatorsTest, SymTensor3D_sub_AxiSymTensor3D) {
    auto res = S - A;
    Mat3x3 expected = ReferenceSub(S.toMat3x3(), A.toMat3x3());
    ExpectTensorEqMat(res, expected);
}

TEST_F(TensorOperatorsTest, SkewSymTensor3D_sub_Tensor3D) {
    auto res = W - T;
    Mat3x3 expected = ReferenceSub(W.toMat3x3(), T.toMat3x3());
    ExpectTensorEqMat(res, expected);
}

TEST_F(TensorOperatorsTest, SkewSymTensor3D_sub_SymTensor3D) {
    auto res = W - S;
    Mat3x3 expected = ReferenceSub(W.toMat3x3(), S.toMat3x3());
    ExpectTensorEqMat(res, expected);
}

TEST_F(TensorOperatorsTest, SkewSymTensor3D_sub_AxiSymTensor3D) {
    auto res = W - A;
    Mat3x3 expected = ReferenceSub(W.toMat3x3(), A.toMat3x3());
    ExpectTensorEqMat(res, expected);
}