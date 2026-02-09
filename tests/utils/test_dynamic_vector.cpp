#include <gtest/gtest.h>
#include "DynamicVector.h"
#include "FixedVector.h"

using namespace cfem;

TEST(DynamicVectorTest, ConstructorAndSize) {
    DynamicVector<double> v(5, 1.23);
    
    EXPECT_EQ(v.size(), 5);
    for (int i = 0; i < v.size(); ++i) {
        EXPECT_DOUBLE_EQ(v[i], 1.23);
    }
}

TEST(DynamicVectorTest, InitializerList) {
    DynamicVector<double> v = {1.0, 2.0, 3.0};
    
    ASSERT_EQ(v.size(), 3);
    EXPECT_DOUBLE_EQ(v[0], 1.0);
    EXPECT_DOUBLE_EQ(v[1], 2.0);
    EXPECT_DOUBLE_EQ(v[2], 3.0);
}

TEST(DynamicVectorTest, ConversionFromFixed) {
    FixedVector<double, 3> fixed = {10.0, 20.0, 30.0};
    DynamicVector<double> dynamic(fixed);
    
    ASSERT_EQ(dynamic.size(), 3);
    EXPECT_DOUBLE_EQ(dynamic[0], 10.0);
    EXPECT_DOUBLE_EQ(dynamic[2], 30.0);
}

TEST(DynamicVectorTest, AdditionAndSubtraction) {
    DynamicVector<double> v1 = {1.0, 2.0};
    DynamicVector<double> v2 = {3.0, 4.0};
    
    DynamicVector<double> sum = v1 + v2;
    EXPECT_DOUBLE_EQ(sum[0], 4.0);
    EXPECT_DOUBLE_EQ(sum[1], 6.0);
    
    v1 += v2;
    EXPECT_DOUBLE_EQ(v1[0], 4.0);
    
    DynamicVector<double> diff = v1 - v2; // (4,6) - (3,4) = (1,2)
    EXPECT_DOUBLE_EQ(diff[0], 1.0);
    EXPECT_DOUBLE_EQ(diff[1], 2.0);
}

TEST(DynamicVectorTest, ScalarOperations) {
    DynamicVector<double> v = {1.0, 2.0, 3.0};
    
    DynamicVector<double> res = v * 2.0;
    EXPECT_DOUBLE_EQ(res[0], 2.0);
    EXPECT_DOUBLE_EQ(res[1], 4.0);
    EXPECT_DOUBLE_EQ(res[2], 6.0);

    res = 2.0*v;
    EXPECT_DOUBLE_EQ(res[0], 2.0);
    EXPECT_DOUBLE_EQ(res[1], 4.0);
    EXPECT_DOUBLE_EQ(res[2], 6.0);

    res = v/2.0;
    EXPECT_DOUBLE_EQ(res[0], 0.5);
    EXPECT_DOUBLE_EQ(res[1], 1.0);
    EXPECT_DOUBLE_EQ(res[2], 1.5);
    
    res /= 2.0;
    EXPECT_DOUBLE_EQ(res[0], 0.25);
    EXPECT_DOUBLE_EQ(res[1], 0.5);
    EXPECT_DOUBLE_EQ(res[2], 0.75);

#ifndef NDEBUG
    EXPECT_THROW(res/cfem::ZERO_TOLERANCE, cfem::AssertionException);
#endif
}

// --- Tests d'Algèbre Linéaire ---

TEST(DynamicVectorTest, LinearAlgebra) {
    DynamicVector<double> v1 = {1.0, 0.0};
    DynamicVector<double> v2 = {0.0, 1.0};
    
    // Produit scalaire de vecteurs orthogonaux
    EXPECT_DOUBLE_EQ(v1.dot(v2), 0.0);
    
    // Norme de (3, 4)
    DynamicVector<double> v3 = {3.0, 4.0};
    EXPECT_DOUBLE_EQ(v3.normSquare(), 25.0);
    EXPECT_DOUBLE_EQ(v3.norm(), 5.0);
}

// --- Test des Assertions (Optionnel) ---

TEST(DynamicVectorTest, AssertionOnSizeMismatch) {
    DynamicVector<double> v1(3);
    DynamicVector<double> v2(2);
    
#ifndef NDEBUG
    EXPECT_THROW(v1 + v2, cfem::AssertionException); 
    EXPECT_THROW(v1 - v2, cfem::AssertionException); 
#endif


}