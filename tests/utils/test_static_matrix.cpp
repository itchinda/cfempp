#include <gtest/gtest.h>
#include "StaticMatrix.h"

/**
 * @brief Unit tests for the StaticMatrix template class.
 */

const double TOL = cfem::ZERO_TOLERANCE;

TEST(StaticMatrixTest, ConstructionAndAccess) {
    // Initializer list construction
    StaticMatrix<double, 2, 2> mat = {1.0, 2.0, 
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

TEST(StaticMatrixTest, BasicOperations) {
    StaticMatrix<double, 3, 3> I;
    I.setIdentity();
    
    for(int i = 0; i < 3; ++i) {
        for(int j = 0; j < 3; ++j) {
            if(i == j) EXPECT_DOUBLE_EQ(I(i,j), 1.0);
            else       EXPECT_DOUBLE_EQ(I(i,j), 0.0);
        }
    }

    StaticMatrix<double, 2, 2> zero;
    zero.setZero();
    EXPECT_DOUBLE_EQ(zero(1,0), 0.0);
}

TEST(StaticMatrixTest, Arithmetic) {
    StaticMatrix<double, 2, 2> A = {1, 1, 1, 1};
    StaticMatrix<double, 2, 2> B = {2, 2, 2, 2};

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

TEST(StaticMatrixTest, Multiplication) {
    // Matrix-Vector Product
    StaticMatrix<double, 2, 3> A = {1, 2, 3, 
                                    4, 5, 6};
    StaticVector<double, 3> v = {1, 0, 2}; // 1*{1,4} + 0*{2,5} + 2*{3,6} = {7, 16}
    
    auto resV = A * v;
    EXPECT_DOUBLE_EQ(resV[0], 7.0);
    EXPECT_DOUBLE_EQ(resV[1], 16.0);

    // Matrix-Matrix Product
    StaticMatrix<double, 2, 2> M1 = {1, 2, 
                                     3, 4};
    StaticMatrix<double, 2, 2> M2 = {5, 6, 
                                     7, 8};
    // [1*5+2*7, 1*6+2*8] = [19, 22]
    // [3*5+4*7, 3*6+4*8] = [43, 50]
    auto M3 = M1 * M2;
    EXPECT_DOUBLE_EQ(M3(0,0), 19.0);
    EXPECT_DOUBLE_EQ(M3(1,1), 50.0);
}

TEST(StaticMatrixTest, Transpose) {
    StaticMatrix<double, 2, 3> A = {1, 2, 3, 
                                    4, 5, 6};
    auto AT = A.transpose();
    
    EXPECT_EQ(AT.rows(), 3);
    EXPECT_EQ(AT.cols(), 2);
    EXPECT_DOUBLE_EQ(AT(2,0), 3.0);
    EXPECT_DOUBLE_EQ(AT(1,1), 5.0);
}

TEST(StaticMatrixTest, JacobianTools) {
    // 2x2 Determinant
    StaticMatrix<double, 2, 2> m2 = {4, 7, 
                                     2, 6};
    EXPECT_NEAR(det(m2), 10.0, TOL); // 24 - 14

    // 3x3 Determinant
    StaticMatrix<double, 3, 3> m3 = {1, 2, 3, 
                                     0, 4, 5, 
                                     1, 0, 6};
    // 1*(24-0) - 2*(0-5) + 3*(0-4) = 24 + 10 - 12 = 22
    EXPECT_NEAR(det(m3), 22.0, TOL);
}