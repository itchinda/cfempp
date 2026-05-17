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
#include "DynamicMatrix.h"
#include "FixedMatrix.h"
#include "DynamicVector.h"

using namespace cfem;

TEST(DynamicMatrixTest, DefaultConstructor) {
    DynamicMatrix<double> M;
    EXPECT_EQ(M.n_rows(), 0);
    EXPECT_EQ(M.n_cols(), 0);
}

TEST(DynamicMatrixTest, SizeConstructor) {
    DynamicMatrix<double> M(3, 4, 1.5);
    EXPECT_EQ(M.n_rows(), 3);
    EXPECT_EQ(M.n_cols(), 4);
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
    
    ASSERT_EQ(dynamic.n_rows(), 2);
    ASSERT_EQ(dynamic.n_cols(), 2);
    EXPECT_DOUBLE_EQ(dynamic(0, 0), 10.0);
    EXPECT_DOUBLE_EQ(dynamic(0, 1), 20.0);
    EXPECT_DOUBLE_EQ(dynamic(1, 0), 30.0);
    EXPECT_DOUBLE_EQ(dynamic(1, 1), 40.0);
}

// --- Tests Arithmétiques ---

TEST(DynamicMatrixTest, IdentityAndZero) {
    DynamicMatrix<double> M(3, 3);
    M.setIdentity();
    
    for (CfemInt i = 0; i < 3; ++i) {
        for (CfemInt j = 0; j < 3; ++j) {
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
    EXPECT_DOUBLE_EQ(div(0, 0), 3.0);
    EXPECT_DOUBLE_EQ(div(0, 1), 3.0);
    EXPECT_DOUBLE_EQ(div(1, 0), 3.0);
    EXPECT_DOUBLE_EQ(div(1, 1), 3.0);
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
    
    ASSERT_EQ(MT.n_rows(), 3);
    ASSERT_EQ(MT.n_cols(), 2);
    EXPECT_DOUBLE_EQ(MT(0, 0), 1.0);
    EXPECT_DOUBLE_EQ(MT(2, 1), 6.0);
}

TEST(DynamicMatrixTest, SpanRowAccess) {
    DynamicMatrix<double> M(3, 3, {1, 2, 3, 
                                   4, 5, 6, 
                                   7, 8, 9});
    
    // Extraction d'une ligne via std::span
    auto row1 = M[1]; // Doit pointer sur {4, 5, 6}
    ASSERT_EQ(row1.size(), 3);
    EXPECT_DOUBLE_EQ(row1[0], 4.0);
    EXPECT_DOUBLE_EQ(row1[2], 6.0);

    // Modification via le span (doit modifier la matrice)
    row1[1] = 50.0;
    EXPECT_DOUBLE_EQ(M(1, 1), 50.0);
}

TEST(DynamicMatrixTest, MatrixViewAndSubView) {
    // Matrice 4x4
    DynamicMatrix<double> M(4, 4, { 1,  2,  3,  4,
                                    5,  6,  7,  8,
                                    9, 10, 11, 12,
                                   13, 14, 15, 16});

    // Extraction d'un bloc 2x2 au centre (lignes 1-2, cols 1-2)
    auto view = M.view(1, 1, 2, 2);

    ASSERT_EQ(view.n_rows(), 2);
    ASSERT_EQ(view.n_cols(), 2);
    EXPECT_EQ(view.stride(), 4); // Le stride doit rester celui de la matrice parent (4)

    // Vérification des valeurs du bloc central
    EXPECT_DOUBLE_EQ(view(0, 0), 6.0);
    EXPECT_DOUBLE_EQ(view(0, 1), 7.0);
    EXPECT_DOUBLE_EQ(view(1, 0), 10.0);
    EXPECT_DOUBLE_EQ(view(1, 1), 11.0);

    // Modification via la vue
    view(0, 0) = -6.0;
    EXPECT_DOUBLE_EQ(M(1, 1), -6.0); // La matrice d'origine doit être impactée
}

TEST(DynamicMatrixTest, AdditionAndSubtraction) {
    DynamicMatrix<double> A(2, 2, {1.0, 2.0, 3.0, 4.0});
    DynamicMatrix<double> B(2, 2, {4.0, 3.0, 2.0, 1.0});

    DynamicMatrix<double> sum = A + B;
    EXPECT_DOUBLE_EQ(sum(0, 0), 5.0);
    EXPECT_DOUBLE_EQ(sum(1, 1), 5.0);

    A -= B; // A devient {-3.0, -1.0, 1.0, 3.0}
    EXPECT_DOUBLE_EQ(A(0, 0), -3.0);
    EXPECT_DOUBLE_EQ(A(0, 1), -1.0);
    EXPECT_DOUBLE_EQ(A(1, 0), 1.0);
    EXPECT_DOUBLE_EQ(A(1, 1), 3.0);
}

TEST(DynamicMatrixTest, CopyAndMoveSemantics) {
    DynamicMatrix<double> orig(2, 2, {1, 2, 3, 4});

    // --- Copy Constructor ---
    DynamicMatrix<double> copy(orig);
    copy(0, 0) = 99.0;
    EXPECT_DOUBLE_EQ(orig(0, 0), 1.0);  // Orig doit être inchangé (Deep Copy)
    EXPECT_DOUBLE_EQ(copy(0, 0), 99.0);

    // --- Move Constructor ---
    DynamicMatrix<double> moved(std::move(orig));
    EXPECT_EQ(orig.n_rows(), 0);        // L'original doit être vidé
    EXPECT_EQ(orig.n_cols(), 0);
    EXPECT_EQ(orig.size(), 0);
    EXPECT_EQ(moved.n_rows(), 2);       // Le nouveau possède les données
    EXPECT_DOUBLE_EQ(moved(0, 1), 2.0);
}

TEST(DynamicMatrixTest, ReshapeAndAssign) {
    DynamicMatrix<double> M(2, 3, 1.0); // 2x3 remplie de 1.0
    EXPECT_EQ(M.size(), 6);

    // Reshape sans changer la mémoire totale
    M.reshape(3, 2); 
    EXPECT_EQ(M.n_rows(), 3);
    EXPECT_EQ(M.n_cols(), 2);

    // Assign (change la taille et réinitialise)
    M.assign(4, 4, 7.0);
    EXPECT_EQ(M.n_rows(), 4);
    EXPECT_EQ(M.n_cols(), 4);
    EXPECT_EQ(M.size(), 16);
    EXPECT_DOUBLE_EQ(M(3, 3), 7.0);
}

TEST(DynamicMatrixTest, ConstMatrixView) {
    const DynamicMatrix<double> M(2, 3, 
        {1, 2, 3,
         4, 5, 6}
    );
    auto view = M.view();
    static_assert( std::is_same_v< decltype(view),  ConstMatrixView<double> > );

    EXPECT_EQ(view.n_rows(), 2);
    EXPECT_EQ(view.n_cols(), 3);

    EXPECT_DOUBLE_EQ(view(1, 2), 6.0);

    auto row = view.row(0);
    EXPECT_DOUBLE_EQ(row[1], 2.0);
}

TEST(DynamicMatrixTest, MatrixViewContiguous) {
    DynamicMatrix<double> M(4, 4, 1.0);

    auto full = M.view();
    EXPECT_TRUE(full.isContiguous());

    auto sub = M.view(1, 1, 2, 2);

    // stride = 4, cols = 2
    EXPECT_FALSE(sub.isContiguous());
}

TEST(DynamicMatrixTest, EmptyMatrix) {
    DynamicMatrix<double> M;

    EXPECT_TRUE(M.empty());
    M.assign(2, 2, 1.0);
    EXPECT_FALSE(M.empty());

    M.clear();

    EXPECT_TRUE(M.empty());
}

TEST(DynamicMatrixTest, SameShape) {
    DynamicMatrix<double> A(2, 3);
    DynamicMatrix<double> B(2, 3);
    DynamicMatrix<double> C(3, 2);

    EXPECT_TRUE(A.sameShape(B));
    EXPECT_FALSE(A.sameShape(C));
}

TEST(DynamicMatrixTest, IsSquare) {
    DynamicMatrix<double> A(3, 3);
    DynamicMatrix<double> B(2, 3);

    EXPECT_TRUE(A.isSquare());
    EXPECT_FALSE(B.isSquare());
}

TEST(DynamicMatrixTest, NestedSubView) {
    DynamicMatrix<double> M(5, 5, 
        { 1,  2,  3,  4,  5,
          6,  7,  8,  9, 10,
         11, 12, 13, 14, 15,
         16, 17, 18, 19, 20,
         21, 22, 23, 24, 25 }
    );

    auto v1 = M.view(1, 1, 3, 3);
    auto v2 = v1.subView(1, 1, 2, 2);

    EXPECT_EQ(v2.n_rows(), 2);
    EXPECT_EQ(v2.n_cols(), 2);

    EXPECT_DOUBLE_EQ(v2(0,0), 13.0);
    EXPECT_DOUBLE_EQ(v2(1,1), 19.0);

    v2(0,0) = -13.0;

    EXPECT_DOUBLE_EQ(M(2,2), -13.0);
}

TEST(DynamicMatrixTest, TransposeAliasing) {
    DynamicMatrix<double> M(2, 2, 
        { 1, 2,
          3, 4 }
    );

    M.transpose(M);

    EXPECT_DOUBLE_EQ(M(0,0), 1.0);
    EXPECT_DOUBLE_EQ(M(0,1), 3.0);
    EXPECT_DOUBLE_EQ(M(1,0), 2.0);
    EXPECT_DOUBLE_EQ(M(1,1), 4.0);
}

TEST(DynamicMatrixTest, ReshapePreservesData) {
    DynamicMatrix<double> M(2, 3, 
        { 1, 2, 3,
          4, 5, 6 }
    );

    M.reshape(3, 2);

    EXPECT_DOUBLE_EQ(M(0,0), 1.0);
    EXPECT_DOUBLE_EQ(M(0,1), 2.0);
    EXPECT_DOUBLE_EQ(M(1,0), 3.0);
    EXPECT_DOUBLE_EQ(M(2,1), 6.0);
}

TEST(DynamicMatrixTest, MoveAssignment) {
    DynamicMatrix<double> A(2, 2,
       { 1, 2,
         3, 4 }
    );

    DynamicMatrix<double> B;

    B = std::move(A);

    EXPECT_EQ(A.size(), 0);

    EXPECT_EQ(B.n_rows(), 2);
    EXPECT_EQ(B.n_cols(), 2);

    EXPECT_DOUBLE_EQ(B(1,1), 4.0);
}

TEST(DynamicMatrixTest, ConstRowAccess) {
    const DynamicMatrix<double> M(2, 2, 
        { 1, 2,
          3, 4 }
    );

    auto row = M[1];
    static_assert( std::is_same_v< decltype(row), std::span<const double> > );
    EXPECT_DOUBLE_EQ(row[0], 3.0);
    EXPECT_DOUBLE_EQ(row[1], 4.0);
}

TEST(DynamicMatrixTest, MatrixViewStrideAccess) {
    DynamicMatrix<double> M(4, 4,
        { 1,  2,  3,  4,
          5,  6,  7,  8,
          9, 10, 11, 12,
         13, 14, 15, 16 }
    );

    auto sub = M.view(1, 1, 2, 2);

    EXPECT_EQ(sub.stride(), 4);

    EXPECT_DOUBLE_EQ(sub(0,0), 6.0);
    EXPECT_DOUBLE_EQ(sub(1,0), 10.0);
}

TEST(DynamicMatrixTest, ZeroSizedMatrix) {
    DynamicMatrix<double> M(0, 0);

    EXPECT_TRUE(M.empty());
    EXPECT_EQ(M.size(), 0);

    auto v = M.view();

    EXPECT_EQ(v.n_rows(), 0);
    EXPECT_EQ(v.n_cols(), 0);
}