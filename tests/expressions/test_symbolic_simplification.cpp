#include <gtest/gtest.h>
#include "MathOperators.h" // Your symbolic headers

using namespace cfem::sym;

class SymbolicSimplificationTest : public ::testing::Test {
protected:
    // Helper to check if expression string matches expected
    void ExpectExpr(const std::shared_ptr<const CfemExpression>& e, const std::string& expected) {
        EXPECT_EQ(e->toString(), expected);
    }
};

// --- 1. Basic Identity & Scale Tests ---
TEST_F(SymbolicSimplificationTest, BasicIdentities) {
    auto u = scalar(x);
    
    // u + u -> 2 * u
    ExpectExpr(u + u, "2 * x");
    
    // 2*u + u -> 3 * u
    ExpectExpr((2.0 * u) + u, "3 * x");
    
    // u + (-u) -> 0
    // (Assuming operator- is implemented as + Scaled(-1))
    ExpectExpr(u - u, "0");
}

// --- 2. The "Simple Example" (Cancellation) ---
TEST_F(SymbolicSimplificationTest, CancellationAndReduction) {
    auto X = scalar(x);
    auto s = sin(X);
    auto c = cos(X);

    // (s + c) + (s - c) -> 2 * s
    // Here, c and -c should cancel out via your peek logic
    ExpectExpr((s + c) + (s - c), "2 * sin(x)");

    // Final target: (s + c) + (s - c) - s -> sin(x)
    ExpectExpr(((s + c) + (s - c)) - s, "sin(x)");
}

// --- 3. Nested Peek Logic (Scrambled Terms) ---
TEST_F(SymbolicSimplificationTest, NestedScrambledTerms) {
    auto u = exp(scalar(x));
    auto v = scalar(y);
    auto w = scalar(z);

    // u + (v + u) -> 2*u + v
    ExpectExpr(u + (v + u), "2 * exp(x) + y");

    // (u + v) + (w - u) -> v + w
    // The 'u' and '-u' are buried inside binary nodes
    ExpectExpr((u + v) + (w - u), "y + z");
}

// --- 4. Unary and Complex Structural Tests ---
TEST_F(SymbolicSimplificationTest, UnaryStructuralIdentity) {
    auto X = scalar(x);
    auto f = sin(exp(X));
    
    // sin(exp(x)) + sin(exp(x)) -> 2 * sin(exp(x))
    ExpectExpr(f + f, "2 * sin(exp(x))");
    
    // 3 * f - f -> 2 * f
    ExpectExpr((3.0 * f) - f, "2 * sin(exp(x))");
}

// --- 5. Stress Test (Zero Merging) ---
TEST_F(SymbolicSimplificationTest, ZeroPropagation) {
    auto X = scalar(x);
    auto complexZero = (sin(X) + cos(X)) - (cos(X) + sin(X));
    
    // Everything should vanish
    EXPECT_TRUE(complexZero->toString() == "0");
}

// // --- 6. Vector Component Matching ---
// TEST_F(SymbolicSimplificationTest, VectorIdentities) {
//     // Assuming you have a vector(u, v, w) constructor
//     auto vecA = vector(x, y, 0.0);
//     auto vecB = vector(x, y, 0.0);

//     // vecA + vecB -> 2 * Vector(u, v, 0)
//     // This tests your isIdentical recursive component check
//     ExpectExpr(vecA + vecB, "2 * [x, y, 0]"); 
// }

TEST_F(SymbolicSimplificationTest, VariableInteraction) {
    CfemReal val = 5.0;
    auto p = std::make_shared<VariableExpression>("p", &val, 1);
    auto X = scalar(x);

    // p + x should not simplify (they are different types)
    auto expr = p + X;
    ExpectExpr(expr, "p + x");

    // p - p should simplify to 0
    // (This requires isIdentical to check the m_dataPtr)
    ExpectExpr(p - p, "0");
}