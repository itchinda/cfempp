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
#include "Vector3D.h"
#include "TensorO2.h"
#include "MathOperators.h"

using namespace cfem;
using namespace cfem::sym;
using namespace cfem::sym::cfemvars;

class AdvancedExpressionTest : public ::testing::Test {
protected:
    EvalContext ctx;
    Vector3D testPoint{1.0, 2.0, 3.0};
    CfemInt testCellId = 42;

    void SetUp() override {
        ctx.physPos = testPoint;
        ctx.currentCellId = testCellId;
    }

    // Helper to evaluate at testPoint and return a full vector
    std::vector<CfemReal> evaluateFull(std::shared_ptr<const CfemExpression> expr) {
        std::vector<CfemReal> out(expr->getDimension(), 0.0);
        expr->evaluate(testCellId, {0,0,0}, out, ctx);
        return out;
    }
    
    // Quick single-component fetcher
    CfemReal quickEval(std::shared_ptr<const CfemExpression> expr, int comp = 0) {
        return evaluateFull(expr)[comp];
    }
};

// --- 1. Evaluation Context Propagation ---
TEST_F(AdvancedExpressionTest, ContextIsPropagated) {
    // Ensure that LambdaExpression correctly receives the context and cell ID
    auto context_aware_lambda = scalar([](const Vector3D& pos) {
        return pos.x + pos.y + pos.z;
    });

    CfemReal expected = testPoint.x + testPoint.y + testPoint.z;
    EXPECT_DOUBLE_EQ(quickEval(context_aware_lambda), expected);
}

// --- 2. Advanced Binary Operator Algebra ---
TEST_F(AdvancedExpressionTest, FullAlgebraicSuite) {
    auto f = scalar(x * x); // at x=1 -> 1
    auto g = scalar(y);     // at y=2 -> 2

    auto add_expr = f + g;
    auto sub_expr = f - g;
    auto mul_expr = f * g;
    auto div_expr = f / g;

    EXPECT_DOUBLE_EQ(quickEval(add_expr), 1.0 + 2.0);
    EXPECT_DOUBLE_EQ(quickEval(sub_expr), 1.0 - 2.0);
    EXPECT_DOUBLE_EQ(quickEval(mul_expr), 1.0 * 2.0);
    EXPECT_DOUBLE_EQ(quickEval(div_expr), 1.0 / 2.0);
}

// --- 3. Vector Math (Dot and Cross Products) ---
TEST_F(AdvancedExpressionTest, VectorProducts) {
    auto u = vector(x, y, z); // [1, 2, 3]
    auto v = vector(z, x, y); // [3, 1, 2]

    auto d_prod = dot(u, v);
    auto c_prod = cross(u, v);

    // Note: This assumes BinaryExpression::evaluate handles MathOp::Dot and Cross
    // u.v = 1*3 + 2*1 + 3*2 = 3 + 2 + 6 = 11
    EXPECT_DOUBLE_EQ(quickEval(d_prod), 11.0);

    // u x v = [(2*2 - 3*1), (3*3 - 1*2), (1*1 - 2*3)] = [1, 7, -5]
    auto cross_res = evaluateFull(c_prod);
    ASSERT_EQ(cross_res.size(), 3);
    EXPECT_DOUBLE_EQ(cross_res[0], 1.0);
    EXPECT_DOUBLE_EQ(cross_res[1], 7.0);
    EXPECT_DOUBLE_EQ(cross_res[2], -5.0);
}

// --- 4. Deep Tree Gradient (Chain & Sum Rules) ---
TEST_F(AdvancedExpressionTest, DeepGradientTree) {
    // h = x^2 * y + z^3
    auto f1 = scalar(x * x);
    auto f2 = scalar(y);
    auto f3 = scalar(z * z * z);
    
    auto h = (f1 * f2) + f3;
    
    // grad(h) = [2xy, x^2, 3z^2]
    auto gh = grad(h);
    
    EXPECT_EQ(gh->getDimension(), 3);
    auto res = evaluateFull(gh);
    
    EXPECT_DOUBLE_EQ(res[0], 2.0 * testPoint.x * testPoint.y); // 2 * 1 * 2 = 4
    EXPECT_DOUBLE_EQ(res[1], testPoint.x * testPoint.x);       // 1^2 = 1
    EXPECT_DOUBLE_EQ(res[2], 3.0 * testPoint.z * testPoint.z); // 3 * 3^2 = 27
}

#ifndef NDEBUG
    TEST_F(AdvancedExpressionTest, SpanBufferOverflowProtection) {
        auto t = tensor([](const Vector3D&){ return TensorO2(); }); // Dim = 9

        // Create a span that is purposefully too small (size 3 instead of 9)
        std::vector<CfemReal> tiny_buffer(3, 0.0);
        std::span<CfemReal> bad_span(tiny_buffer.data(), tiny_buffer.size());

        // This should trigger your CFEM_ASSERT: "LambdaExpression: tensor output span too small."
        // EXPECT_DEATH checks that the program aborts with a specific message.
        EXPECT_THROW( (t->evaluate(testCellId, {0,0,0}, bad_span, ctx)),cfem::AssertionException);
    }
#endif

TEST_F(AdvancedExpressionTest, InvalidGradientDimensionAssertion) {
    auto bad_grad = scalar([](const Vector3D& X){ return X.x; }); // Dim 1

    // Vérifie simplement que le code jette une exception
    EXPECT_ANY_THROW({
        auto invalid_vec = vector([](const Vector3D& X){ return X; }, bad_grad);
    });
}

TEST_F(AdvancedExpressionTest, TensorRowMajorMapping) {
    // We want to ensure that LambdaExpression maps a 3x3 TensorO2 to a 9-element 
    // span in the exact order the DomainIntegrator expects (usually row-major).
    auto analytical_grad_u = tensor([](const Vector3D& X) -> TensorO2 {
        // T = [ 1, 2, 3 ]
        //     [ 4, 5, 6 ]
        //     [ 7, 8, 9 ]
        return TensorO2({1., 2., 3., 4., 5., 6., 7., 8., 9.});
    });

    auto res = evaluateFull(analytical_grad_u);
    ASSERT_EQ(res.size(), 9);
    
    // Check linear memory layout
    for(int i = 0; i < 9; ++i) {
        EXPECT_DOUBLE_EQ(res[i], static_cast<double>(i + 1));
    }
}