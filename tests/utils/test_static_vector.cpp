#include <gtest/gtest.h>
#include "StaticVector.h"

/**
 * @brief Unit tests for the StaticVector template class.
 */

const double TOL = cfem::ZERO_TOLERANCE;

TEST(StaticVectorTest, ConstructionAndInitialization) {
    // Default constructor
    StaticVector<double, 4> v1;
    for (int i = 0; i < 4; ++i) {
        EXPECT_NEAR(v1[i], 0.0, TOL);
    }

    // Initializer list
    StaticVector<double, 3> v2 = {1.0, 2.0, 3.0};
    EXPECT_NEAR(v2[0], 1.0, TOL);
    EXPECT_NEAR(v2[1], 2.0, TOL);
    EXPECT_NEAR(v2[2], 3.0, TOL);
}

TEST(StaticVectorTest, BasicArithmetic) {
    StaticVector<double, 2> a = {1.0, 2.0};
    StaticVector<double, 2> b = {3.0, 4.0};

    // Addition
    StaticVector<double, 2> sum = a + b;
    EXPECT_NEAR(sum[0], 4.0, TOL);
    EXPECT_NEAR(sum[1], 6.0, TOL);

    // Subtraction
    StaticVector<double, 2> diff = b - a;
    EXPECT_NEAR(diff[0], 2.0, TOL);
    EXPECT_NEAR(diff[1], 2.0, TOL);

    // Scalar Multiplication (both sides)
    StaticVector<double, 2> mult1 = a * 10.0;
    StaticVector<double, 2> mult2 = 10.0 * a;
    EXPECT_NEAR(mult1[0], 10.0, TOL);
    EXPECT_NEAR(mult2[1], 20.0, TOL);

    // Division
    StaticVector<double, 2> div = b / 2.0;
    EXPECT_NEAR(div[0], 1.5, TOL);
    EXPECT_NEAR(div[1], 2.0, TOL);
}

TEST(StaticVectorTest, CompoundAssignment) {
    StaticVector<double, 3> v = {1.0, 1.0, 1.0};
    StaticVector<double, 3> step = {1.0, 2.0, 3.0};

    v += step; // {2, 3, 4}
    EXPECT_NEAR(v[2], 4.0, TOL);

    v *= 0.5;  // {1, 1.5, 2}
    EXPECT_NEAR(v[1], 1.5, TOL);
}

TEST(StaticVectorTest, LinearAlgebra) {
    StaticVector<double, 3> v1 = {1.0, 2.0, 3.0};
    StaticVector<double, 3> v2 = {0.0, 1.0, 0.5};

    // Dot Product: (1*0) + (2*1) + (3*0.5) = 2 + 1.5 = 3.5
    EXPECT_NEAR(v1.dot(v2), 3.5, TOL);

    // Norm: Vector {3, 4} -> Norm 5
    StaticVector<double, 2> v3 = {3.0, 4.0};
    EXPECT_NEAR(v3.norm(), 5.0, TOL);
    EXPECT_NEAR(v3.normSquare(), 25.0, TOL);
}

TEST(StaticVectorTest, Vector3DHandshake) {
    // Test conversion from Vector3D to StaticVector<T, 3>
    Vector3D geoP = {10.5, 20.5, 30.5};
    StaticVector<double, 3> algebraicV = geoP;
    
    EXPECT_NEAR(algebraicV[0], 10.5, TOL);
    EXPECT_NEAR(algebraicV[2], 30.5, TOL);

    // Test conversion from StaticVector<T, 3> back to Vector3D
    StaticVector<double, 3> sv = {1.0, -1.0, 5.0};
    Vector3D p = sv; // Implicit cast
    
    EXPECT_NEAR(p.x, 1.0, TOL);
    EXPECT_NEAR(p.y, -1.0, TOL);
    EXPECT_NEAR(p.z, 5.0, TOL);
}

TEST(StaticVectorTest, Printing) {
    StaticVector<double, 2> v = {1.1, 2.2};
    std::stringstream ss;
    ss << v;
    EXPECT_EQ(ss.str(), "(1.1, 2.2)");
}