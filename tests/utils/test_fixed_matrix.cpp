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

#include "FixedMatrix.h"

using namespace cfem;

/**
 * @brief Unit tests for the FixedMatrix template class.
 */

const CfemReal TOL = cfem::ZERO_TOLERANCE;

TEST(FixedMatrixTest, ConstructionAndAccess) {
    // Initializer list construction
    FixedMatrix<CfemReal, 2, 2> mat = {1.0, 2.0, 
                                      3.0, 4.0};
    
    EXPECT_EQ(mat.rows(), 2);
    EXPECT_EQ(mat.cols(), 2);
    EXPECT_DOUBLE_EQ(mat(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(mat(1, 1), 4.0);

    // Dimension pair check
    auto dims = mat.dimensions();
    EXPECT_EQ(dims.first, 2);
    EXPECT_EQ(dims.second, 2);
}

TEST(FixedMatrixTest, BasicOperations) {
    FixedMatrix<CfemReal, 3, 3> I;
    I.setIdentity();
    
    for(CfemInt i = 0; i < 3; ++i) {
        for(CfemInt j = 0; j < 3; ++j) {
            if(i == j) EXPECT_DOUBLE_EQ(I(i,j), 1.0);
            else       EXPECT_DOUBLE_EQ(I(i,j), 0.0);
        }
    }

    FixedMatrix<CfemReal, 2, 2> zero;
    zero.setZero();
    EXPECT_DOUBLE_EQ(zero(1,0), 0.0);
}

TEST(FixedMatrixTest, Arithmetic) {
    FixedMatrix<CfemReal, 2, 2> A = {1, 1, 1, 1};
    FixedMatrix<CfemReal, 2, 2> B = {2, 2, 2, 2};

    // Addition & Subtraction
    auto C = A + B;
    EXPECT_DOUBLE_EQ(C(0,0), 3.0);
    
    auto D = B - A;
    EXPECT_DOUBLE_EQ(D(1,1), 1.0);

    // Scalar mult/div
    auto E = A * 5.0;
    EXPECT_DOUBLE_EQ(E(0,1), 5.0);
    
    auto F = B / 2.0;
    EXPECT_DOUBLE_EQ(F(0,0), 1.0);
}

TEST(FixedMatrixTest, Multiplication) {
    // Matrix-Vector Product
    FixedMatrix<CfemReal, 2, 3> A = {1, 2, 3, 
                                    4, 5, 6};
    FixedVector<CfemReal, 3> v = {1, 0, 2}; // 1*{1,4} + 0*{2,5} + 2*{3,6} = {7, 16}
    
    auto resV = A * v;
    EXPECT_DOUBLE_EQ(resV[0], 7.0);
    EXPECT_DOUBLE_EQ(resV[1], 16.0);

    // Matrix-Matrix Product
    FixedMatrix<CfemReal, 2, 2> M1 = {1, 2, 
                                     3, 4};
    FixedMatrix<CfemReal, 2, 2> M2 = {5, 6, 
                                     7, 8};
    // [1*5+2*7, 1*6+2*8] = [19, 22]
    // [3*5+4*7, 3*6+4*8] = [43, 50]
    auto M3 = M1 * M2;
    EXPECT_DOUBLE_EQ(M3(0,0), 19.0);
    EXPECT_DOUBLE_EQ(M3(1,1), 50.0);
}

TEST(FixedMatrixTest, InPlaceArithmetic)
{
    FixedMatrix<CfemReal, 2, 2> A = {1,2,3,4};

    A += FixedMatrix<CfemReal,2,2>{1,1,1,1};
    EXPECT_DOUBLE_EQ(A(0,0), 2.0);
    EXPECT_DOUBLE_EQ(A(1,1), 5.0);

    A -= FixedMatrix<CfemReal,2,2>{1,1,1,1};
    EXPECT_DOUBLE_EQ(A(0,0), 1.0);
    EXPECT_DOUBLE_EQ(A(1,1), 4.0);

    A *= 2.0;
    EXPECT_DOUBLE_EQ(A(0,1), 4.0);

    A /= 2.0;
    EXPECT_DOUBLE_EQ(A(1,0), 3.0);
}

TEST(FixedMatrixTest, IdentityMultiplication)
{
    FixedMatrix<CfemReal,3,3> I;
    I.setIdentity();

    FixedMatrix<CfemReal,3,3> A = {
        1,2,3, 4,5,6, 7,8,9
    };

    auto R1 = I * A;
    auto R2 = A * I;

    for(CfemInt i=0;i<3;++i)
        for(CfemInt j=0;j<3;++j)
        {
            EXPECT_DOUBLE_EQ(R1(i,j), A(i,j));
            EXPECT_DOUBLE_EQ(R2(i,j), A(i,j));
        }
}

TEST(FixedMatrixTest, Transpose) {
    FixedMatrix<CfemReal, 2, 3> A = {1, 2, 3, 
                                    4, 5, 6};
    auto AT = A.transpose();
    
    EXPECT_EQ(AT.rows(), 3);
    EXPECT_EQ(AT.cols(), 2);
    EXPECT_DOUBLE_EQ(AT(2,0), 3.0);
    EXPECT_DOUBLE_EQ(AT(1,1), 5.0);
}

TEST(FixedMatrixTest, DoubleTranspose)
{
    FixedMatrix<CfemReal,2,3> A = {
        1,2,3, 4,5,6
    };

    auto ATT = A.transpose().transpose();
    for(CfemInt i=0;i<2;++i)
        for(CfemInt j=0;j<3;++j)
            EXPECT_DOUBLE_EQ(ATT(i,j), A(i,j));
}

TEST(FixedMatrixTest, FrobeniusNorm)
{
    FixedMatrix<CfemReal,2,2> A = { 1,2, 3,4 };

    CfemReal expected = std::sqrt(1+4+9+16);
    EXPECT_NEAR(A.frobeniusNorm(), expected, TOL);
}

TEST(FixedMatrixTest, Inverse3x3)
{
    FixedMatrix<CfemReal,3,3> A = {
        4,7,2, 3,6,1, 2,5,1
    };

    auto invA = inverse(A);
    auto I = A * invA;
    for(CfemInt i=0;i<3;++i) {
        for(CfemInt j=0;j<3;++j) {
            if(i==j)
                EXPECT_NEAR(I(i,j), 1.0, 1e-10);
            else
                EXPECT_NEAR(I(i,j), 0.0, 1e-10);
        }
    }
}

TEST(FixedMatrixTest, JacobianTools) {
    // 2x2 Determinant
    FixedMatrix<CfemReal, 2, 2> m2 = {4, 7, 
                                     2, 6};
    EXPECT_NEAR(det(m2), 10.0, TOL); // 24 - 14

    // 3x3 Determinant
    FixedMatrix<CfemReal, 3, 3> m3 = {1, 2, 3, 
                                     0, 4, 5, 
                                     1, 0, 6};
    // 1*(24-0) - 2*(0-5) + 3*(0-4) = 24 + 10 - 12 = 22
    EXPECT_NEAR(det(m3), 22.0, TOL);
}

TEST(FixedMatrixTest, SolveLinearSystem)
{
    FixedMatrix<CfemReal,3,3> A = {
        3,2,-1,
        2,-2,4,
        -1,0.5,-1
    };

    FixedVector<CfemReal,3> b = {1,-2,0};
    FixedVector<CfemReal,3> x;

    solve(A,b,x);

    EXPECT_NEAR(x[0], 1.0, 1e-10);
    EXPECT_NEAR(x[1], -2.0, 1e-10);
    EXPECT_NEAR(x[2], -2.0, 1e-10);
}

TEST(FixedMatrixTest, SolveRequiresPivoting)
{
    FixedMatrix<CfemReal,2,2> A = { 0,1,  1,1 };

    FixedVector<CfemReal,2> b = {1,2};
    FixedVector<CfemReal,2> x;

    solve(A,b,x);

    EXPECT_NEAR(x[0], 1.0, 1e-12);
    EXPECT_NEAR(x[1], 1.0, 1e-12);
}

TEST(FixedMatrixTest, OuterProduct)
{
    FixedVector<CfemReal,3> a = {1,2,3};
    FixedVector<CfemReal,3> b = {4,5,6};

    auto M = outerProduct(a,b);

    EXPECT_DOUBLE_EQ(M(0,0), 4.0);
    EXPECT_DOUBLE_EQ(M(1,2), 12.0);
    EXPECT_DOUBLE_EQ(M(2,1), 15.0);
}

TEST(FixedMatrixTest, BlockExtraction)
{
    FixedMatrix<CfemReal,4,4> A = {
        1,2,3,4,
        5,6,7,8,
        9,10,11,12,
        13,14,15,16
    };

    auto B = A.block<2,2>(1,1);

    EXPECT_DOUBLE_EQ(B(0,0), 6.0);
    EXPECT_DOUBLE_EQ(B(0,1), 7.0);
    EXPECT_DOUBLE_EQ(B(1,0), 10.0);
    EXPECT_DOUBLE_EQ(B(1,1), 11.0);
}

TEST(FixedMatrixTest, NearSingularSolve)
{
    FixedMatrix<CfemReal,2,2> A = {
        1,1,
        1,1.000000001
    };

    FixedVector<CfemReal,2> b = {2,2.000000001};
    FixedVector<CfemReal,2> x;

    solve(A,b,x);

    EXPECT_NEAR(x[0],1.0,1e-6);
    EXPECT_NEAR(x[1],1.0,1e-6);
}