#include <gtest/gtest.h>
#include "DynamicMatrix.h"
#include "FixedMatrix.h"
#include "DynamicVector.h"

TEST(DynamicMatrixTest, DefaultConstructor) {
    DynamicMatrix<double> M;
    EXPECT_EQ(M.rows(), 0);
    EXPECT_EQ(M.cols(), 0);
}

TEST(DynamicMatrixTest, SizeConstructor) {
    DynamicMatrix<double> M(3, 4, 1.5);
    EXPECT_EQ(M.rows(), 3);
    EXPECT_EQ(M.cols(), 4);
    EXPECT_DOUBLE_EQ(M(0, 0), 1.5);
    EXPECT_DOUBLE_EQ(M(2, 3), 1.5);
}

TEST(DynamicMatrixTest, InitializerList) {
    // Matrix 2x2 : | 1.0  2.0 |
    //               | 3.0  4.0 |
    DynamicMatrix<double> M(2, 2, {1.0, 2.0, 3.0, 4.0});
    EXPECT_DOUBLE_EQ(M(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(M(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(M(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(M(1, 1), 4.0);
}

TEST(DynamicMatrixTest, ConversionFromFixed) {
    FixedMatrix<double, 2, 2> fixed = {10.0, 20.0, 30.0, 40.0};
    DynamicMatrix<double> dynamic(fixed);
    
    ASSERT_EQ(dynamic.rows(), 2);
    ASSERT_EQ(dynamic.cols(), 2);
    EXPECT_DOUBLE_EQ(dynamic(0, 0), 10.0);
    EXPECT_DOUBLE_EQ(dynamic(0, 1), 20.0);
    EXPECT_DOUBLE_EQ(dynamic(1, 0), 30.0);
    EXPECT_DOUBLE_EQ(dynamic(1, 1), 40.0);
}

// --- Tests Arithmétiques ---

TEST(DynamicMatrixTest, IdentityAndZero) {
    DynamicMatrix<double> M(3, 3);
    M.setIdentity();
    
    for (int i = 0; i < 3; ++i) {
        for (int j = 0; j < 3; ++j) {
            if (i == j) EXPECT_DOUBLE_EQ(M(i, j), 1.0);
            else        EXPECT_DOUBLE_EQ(M(i, j), 0.0);
        }
    }
}

TEST(DynamicMatrixTest, ScalarOperations) {
    DynamicMatrix<double> M(2, 2, 2.0);
    M *= 3.0;
    EXPECT_DOUBLE_EQ(M(0, 0), 6.0);
    EXPECT_DOUBLE_EQ(M(0, 1), 6.0);
    EXPECT_DOUBLE_EQ(M(1, 0), 6.0);
    EXPECT_DOUBLE_EQ(M(1, 1), 6.0);
    
    DynamicMatrix<double> div = M / 2.0;
    EXPECT_DOUBLE_EQ(div(0, 0), 2.0);
    EXPECT_DOUBLE_EQ(div(0, 1), 2.0);
    EXPECT_DOUBLE_EQ(div(1, 0), 2.0);
    EXPECT_DOUBLE_EQ(div(1, 1), 2.0);
}

// --- Tests Algèbre Linéaire ---

TEST(DynamicMatrixTest, MatrixVectorProduct) {
    // M = | 1  2 |   v = | 5 |   res = | 1*5 + 2*6 | = | 17 |
    //     | 3  4 |       | 6 |         | 3*5 + 4*6 |   | 39 |
    DynamicMatrix<double> M(2, 2, {1.0, 2.0, 3.0, 4.0});
    DynamicVector<double> v = {5.0, 6.0};
    
    DynamicVector<double> res = M * v;
    
    ASSERT_EQ(res.size(), 2);
    EXPECT_DOUBLE_EQ(res[0], 17.0);
    EXPECT_DOUBLE_EQ(res[1], 39.0);
}

TEST(DynamicMatrixTest, MatrixMatrixProduct) {
    // A = | 1 2 |   B = | 5 6 |   C = | 1*5+2*7  1*6+2*8 | = | 19 22 |
    //     | 3 4 |       | 7 8 |       | 3*5+4*7  3*6+4*8 |   | 43 50 |
    DynamicMatrix<double> A(2, 2, {1.0, 2.0, 3.0, 4.0});
    DynamicMatrix<double> B(2, 2, {5.0, 6.0, 7.0, 8.0});
    
    DynamicMatrix<double> C = A * B;
    
    EXPECT_DOUBLE_EQ(C(0, 0), 19.0);
    EXPECT_DOUBLE_EQ(C(0, 1), 22.0);
    EXPECT_DOUBLE_EQ(C(1, 0), 43.0);
    EXPECT_DOUBLE_EQ(C(1, 1), 50.0);
}

TEST(DynamicMatrixTest, Transpose) {
    // M = | 1 2 3 |  MT = | 1 4 |
    //     | 4 5 6 |       | 2 5 |
    //                     | 3 6 |
    DynamicMatrix<double> M(2, 3, {1, 2, 3, 4, 5, 6});
    DynamicMatrix<double> MT = M.transpose();
    
    ASSERT_EQ(MT.rows(), 3);
    ASSERT_EQ(MT.cols(), 2);
    EXPECT_DOUBLE_EQ(MT(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(MT(2, 1), 6.0);
}