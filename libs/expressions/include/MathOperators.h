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

#pragma once

#include <memory>
#include "CfemExpression.h"
#include "ConstantExpression.h"
#include "Argument.h"
// #include "FEField.h"
// #include "FESpace.h"

#include "symengine/SymEngineSymbols.h"

namespace cfem
{
    class FEField;
    class FESpace;
}
namespace cfem::sym
{

    
    /**
     * @name Fundamental Linear Algebra Operations
     * @{
     */
    std::shared_ptr<const CfemExpression> identityO2TensorExpr();
    std::shared_ptr<const CfemExpression> transpose(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> trace(std::shared_ptr<const CfemExpression> u);
    /** @} */

    /**
     * @name Products and Multilinear Operations
     * @{
     */
    std::shared_ptr<const CfemExpression> dot(std::shared_ptr<const CfemExpression> lhs, 
                                              std::shared_ptr<const CfemExpression> rhs);
    std::shared_ptr<const CfemExpression> inner(std::shared_ptr<const CfemExpression> lhs, 
                                                std::shared_ptr<const CfemExpression> rhs);

    std::shared_ptr<const CfemExpression> cross(std::shared_ptr<const CfemExpression> lhs, 
                                                std::shared_ptr<const CfemExpression> rhs);
    /** @} */

    /**
     * @name Pointwise Mathematical Functions
     * @{
     */
    std::shared_ptr<const CfemExpression> atan2(std::shared_ptr<const CfemExpression> x, 
                                                std::shared_ptr<const CfemExpression> y);
    std::shared_ptr<const CfemExpression> pow(std::shared_ptr<const CfemExpression> base, 
                                              std::shared_ptr<const CfemExpression> expr);

    std::shared_ptr<const CfemExpression> abs(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> sqrt(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> cbrt(std::shared_ptr<const CfemExpression> u);

    std::shared_ptr<const CfemExpression> sign(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> heaviside(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> sin(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> cos(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> tan(std::shared_ptr<const CfemExpression> u);

    std::shared_ptr<const CfemExpression> cot(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> sec(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> csc(std::shared_ptr<const CfemExpression> u);

    std::shared_ptr<const CfemExpression> asin(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> acos(std::shared_ptr<const CfemExpression> u);

    std::shared_ptr<const CfemExpression> atan(std::shared_ptr<const CfemExpression> u);

    std::shared_ptr<const CfemExpression> acot(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> asec(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> acsc(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> sinh(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> cosh(std::shared_ptr<const CfemExpression> u);

    std::shared_ptr<const CfemExpression> tanh(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> asinh(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> acosh(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> atanh(std::shared_ptr<const CfemExpression> u);

    std::shared_ptr<const CfemExpression> exp(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> log(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> erf(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> erfc(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression>        tgamma(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> zeta(std::shared_ptr<const CfemExpression> u);
    /** @} */

    // SymEngine:
    /**
     * @name Common math functions
     *  @{
     */
    SymEngine::Expression sin(const SymEngine::Expression &e);
    SymEngine::Expression cos(const SymEngine::Expression &e);
    SymEngine::Expression tan(const SymEngine::Expression &e);
    SymEngine::Expression cot(const SymEngine::Expression &e);
    SymEngine::Expression sec(const SymEngine::Expression &e);
    SymEngine::Expression csc(const SymEngine::Expression &e);

    SymEngine::Expression asin(const SymEngine::Expression &e);
    SymEngine::Expression acos(const SymEngine::Expression &e);
    SymEngine::Expression atan(const SymEngine::Expression &e);
    SymEngine::Expression acot(const SymEngine::Expression &e);
    SymEngine::Expression asec(const SymEngine::Expression &e);
    SymEngine::Expression acsc(const SymEngine::Expression &e);

    SymEngine::Expression sinh(const SymEngine::Expression &e);
    SymEngine::Expression cosh(const SymEngine::Expression &e);
    SymEngine::Expression tanh(const SymEngine::Expression &e);
    SymEngine::Expression coth(const SymEngine::Expression &e);
    SymEngine::Expression asinh(const SymEngine::Expression &e);
    SymEngine::Expression acosh(const SymEngine::Expression &e);
    SymEngine::Expression atanh(const SymEngine::Expression &e);

    SymEngine::Expression exp(const SymEngine::Expression &e);
    SymEngine::Expression log(const SymEngine::Expression &e);
    SymEngine::Expression ln (const SymEngine::Expression &e);

    SymEngine::Expression sqrt(const SymEngine::Expression &e);
    SymEngine::Expression cbrt(const SymEngine::Expression &e);
    SymEngine::Expression pow(const SymEngine::Expression &b, const SymEngine::Expression &e);

    SymEngine::Expression abs(const SymEngine::Expression &e);

    SymEngine::Expression max(const SymEngine::Expression &a, const SymEngine::Expression &b);
    SymEngine::Expression min(const SymEngine::Expression &a, const SymEngine::Expression &b);

    SymEngine::Expression gamma(const SymEngine::Expression &e);
    SymEngine::Expression erf(const SymEngine::Expression &e);
    SymEngine::Expression erfc(const SymEngine::Expression &e);
    SymEngine::Expression zeta(const SymEngine::Expression &e);
    /** @} */

    /**
     * @name Basic Arithmetic Operators
     * @{
     */
    std::shared_ptr<const CfemExpression> operator-(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression> operator-(std::shared_ptr<const CfemExpression> lhs, 
                                                    std::shared_ptr<const CfemExpression> rhs);

    std::shared_ptr<const CfemExpression> operator+(std::shared_ptr<const CfemExpression> lhs, 
                                                    std::shared_ptr<const CfemExpression> rhs);

    std::shared_ptr<const CfemExpression> operator*(std::shared_ptr<const CfemExpression> lhs, 
                                                    std::shared_ptr<const CfemExpression> rhs);
    std::shared_ptr<const CfemExpression> operator/(std::shared_ptr<const CfemExpression> lhs, 
                                                    std::shared_ptr<const CfemExpression> rhs);
    /** @} */

    namespace converter{
        template<typename T>
        inline std::shared_ptr<const CfemExpression> to_expr(const T& val) {
            return std::make_shared<const ConstantExpression>(val);
        }
    }

    #define CFEM_BINARY_OP_OVERLOADS(op) \
    inline std::shared_ptr<const CfemExpression> operator op (std::shared_ptr<const CfemExpression> lhs, CfemReal rhs) { return lhs op converter::to_expr(rhs); } \
    inline std::shared_ptr<const CfemExpression> operator op (CfemReal lhs, std::shared_ptr<const CfemExpression> rhs) { return converter::to_expr(lhs) op rhs; } \
    inline std::shared_ptr<const CfemExpression> operator op (std::shared_ptr<const CfemExpression> lhs, const Vector3D& rhs) { return lhs op converter::to_expr(rhs); } \
    inline std::shared_ptr<const CfemExpression> operator op (const Vector3D& lhs, std::shared_ptr<const CfemExpression> rhs) { return converter::to_expr(lhs) op rhs; } \
    inline std::shared_ptr<const CfemExpression> operator op (std::shared_ptr<const CfemExpression> lhs, const TensorO2& rhs) { return lhs op converter::to_expr(rhs); } \
    inline std::shared_ptr<const CfemExpression> operator op (const TensorO2& lhs, std::shared_ptr<const CfemExpression> rhs) { return converter::to_expr(lhs) op rhs; }

    CFEM_BINARY_OP_OVERLOADS(+)
    CFEM_BINARY_OP_OVERLOADS(-)
    CFEM_BINARY_OP_OVERLOADS(*)
    CFEM_BINARY_OP_OVERLOADS(/)

    #undef CFEM_BINARY_OP_OVERLOADS


    #define CFEM_GENERATE_BINARY_FUNCTION(func) \
    inline std::shared_ptr<const CfemExpression> func(std::shared_ptr<const CfemExpression> y, CfemReal x) { return func(y, converter::to_expr(x)); } \
    inline std::shared_ptr<const CfemExpression> func(CfemReal y, std::shared_ptr<const CfemExpression> x) { return func(converter::to_expr(y), x); }

    CFEM_GENERATE_BINARY_FUNCTION(atan2)
    CFEM_GENERATE_BINARY_FUNCTION(pow)
    CFEM_GENERATE_BINARY_FUNCTION(max)
    CFEM_GENERATE_BINARY_FUNCTION(min)
    CFEM_GENERATE_BINARY_FUNCTION(dot)
    CFEM_GENERATE_BINARY_FUNCTION(inner)
    #undef CFEM_GENERATE_BINARY_FUNCTION


    /**
     * @name Scalar Expression Factories (Dim 1)
     * @brief Helpers to wrap native C++ lambdas or symbolic variables into 1D CFEM expressions.
     * @{
     */

    /**
     * @brief Wrap a pure spatial scalar lambda f(X) as a CfemExpression.
     * @tparam Func Deduced type of the lambda/functor.
     * @param f Lambda function returning a CfemReal from a Vector3D point.
     * @param gradient (Optional) Expression representing the analytical gradient of f (dimension 3).
     */
    template <typename Func>
    requires std::is_invocable_r_v<CfemReal, Func, const Vector3D&>
    inline std::shared_ptr<const CfemExpression> scalar(
                                                    Func&& f,
                                                    std::shared_ptr<const CfemExpression> gradient = nullptr
                                                );

    /**
     * @brief Wrap a spatio-temporal scalar lambda f(X, t) as a CfemExpression.
     * @tparam Func Deduced type of the lambda/functor.
     * @param f Lambda function returning a CfemReal from a Vector3D and a CfemReal (time).
     * @param gradient (Optional) Expression representing the analytical gradient of f (dimension 3).
     * @param timeDeriv (Optional) Expression representing the time derivative of f (dimension 1).
     */
    template <typename Func>
    requires std::is_invocable_r_v<CfemReal, Func, const Vector3D&, CfemReal>
    inline std::shared_ptr<const CfemExpression> scalar(
                                                    Func&& f,
                                                    std::shared_ptr<const CfemExpression> gradient,
                                                    std::shared_ptr<const CfemExpression> timeDeriv
                                                );
 
    /** 
     * @brief Casts or promotes an existing generic expression into a strictly scalar one. 
     */
    std::shared_ptr<const CfemExpression> scalar(std::shared_ptr<const CfemExpression> expr);

    /** 
     * @brief Generates a constant scalar expression.
     */
    std::shared_ptr<const CfemExpression> constant(CfemReal c);

    /**
     * @brief Convert a SymEngine scalar expression into a CFEM expression.
     * @param ex The scalar symbolic expression
     * @return Shared pointer to a ScalarSymbolicExpression
     */
    std::shared_ptr<const CfemExpression> scalar(const SymEngine::Expression &ex);
    /** @} */

    
    /**
     * @name Vector Expression Factories (Dim 3)
     * @brief Helpers to wrap lambdas or components into 3D CFEM expressions.
     * @{
     */

    /**
     * @brief Wrap a pure spatial vector lambda f(X) as a CfemExpression.
     * @param f Lambda function returning a Vector3D from a Vector3D.
     * @param gradient (Optional) Expression representing the analytical Jacobian of f (dimension 9).
     */
    template <typename Func>
    requires std::is_invocable_r_v<Vector3D, Func, const Vector3D&>
    inline std::shared_ptr<const CfemExpression> vector(
                                                    Func&& f,
                                                    std::shared_ptr<const CfemExpression> gradient = nullptr
                                                );
                            
    /**
     * @brief Wrap a spatio-temporal vector lambda f(X, t) as a CfemExpression.
     * @param f Lambda function returning a Vector3D from a Vector3D and a CfemReal (time).
     * @param gradient (Optional) Expression representing the analytical Jacobian of f (dimension 9).
     * @param timeDeriv (Optional) Expression representing the time derivative of f (dimension 3).
     */
    template <typename Func>
    requires std::is_invocable_r_v<Vector3D, Func, const Vector3D&, CfemReal>
    inline std::shared_ptr<const CfemExpression> vector(
                                                    Func&& f,
                                                    std::shared_ptr<const CfemExpression> gradient,
                                                    std::shared_ptr<const CfemExpression> timeDeriv
                                                );

    /**
     * @brief Casts or promotes an existing expression into a vector expression.
     */
    std::shared_ptr<const CfemExpression> vector(std::shared_ptr<const CfemExpression> expr);
    
    /** 
     * @brief Assemble three independent scalar expressions to form a 3D vector.
     */
    std::shared_ptr<const CfemExpression> vector(std::shared_ptr<const CfemExpression> expr_x,
                                                 std::shared_ptr<const CfemExpression> expr_y,
                                                 std::shared_ptr<const CfemExpression> expr_z);                                      

    /** 
     * @brief Generates a constant vector expression.
     */
    std::shared_ptr<const CfemExpression> constant(const Vector3D& constantVector);

    /**
     * @brief Convert three SymEngine expressions into a vector CFEM expression.
     * @param ex Expression for x-component
     * @param ey Expression for y-component
     * @param ez Expression for z-component
     * @return Shared pointer to a CfemExpression (VectorSymbolicExpression)
     */
    std::shared_ptr<const CfemExpression> vector(const SymEngine::Expression &ex,
                                                 const SymEngine::Expression &ey,
                                                 const SymEngine::Expression &ez);

    /** @} */

    /**
     * @name Tensor Expression Factories (Dim 9)
     * @brief Helpers to wrap lambdas or components into 2nd-order tensor CFEM expressions.
     * @{
     */

    /**
     * @brief Wrap a pure spatial matrix lambda f(X) as a CfemExpression.
     * @param f Lambda returning a 2nd-order tensor (TensorO2) from a Vector3D.
     * @param gradient (Optional) Expression of the 3rd-order gradient tensor (dimension 27).
     */
    template <typename Func>
    requires std::is_invocable_r_v<TensorO2, Func, const Vector3D&>
    inline std::shared_ptr<const CfemExpression> tensor(
                                                    Func&& f,
                                                    std::shared_ptr<const CfemExpression> gradient = nullptr
                                                );

    /**
     * @brief Wrap a spatio-temporal matrix lambda f(X, t) as a CfemExpression.
     * @param f Lambda returning a TensorO2 from a Vector3D and a timestamp.
     * @param gradient (Optional) Expression of the 3rd-order gradient tensor (dimension 27).
     * @param timeDeriv (Optional) Expression of the matrix time derivative (dimension 9).
     */
    template <typename Func>
    requires std::is_invocable_r_v<TensorO2, Func, const Vector3D&, CfemReal>
    inline std::shared_ptr<const CfemExpression> tensor(
                                                    Func&& f,
                                                    std::shared_ptr<const CfemExpression> gradient,
                                                    std::shared_ptr<const CfemExpression> timeDeriv
                                                );
                                        
    /** 
     * @brief Casts or promotes an existing expression into a tensor expression. 
     */
    std::shared_ptr<const CfemExpression> tensor(std::shared_ptr<const CfemExpression> expr);
    
    /** 
     * @brief Assemble 9 independent scalar expressions into a full 2nd-order tensor. 
     */
    std::shared_ptr<const CfemExpression> tensor(std::shared_ptr<const CfemExpression> expr_xx,
                                                 std::shared_ptr<const CfemExpression> expr_xy,
                                                 std::shared_ptr<const CfemExpression> expr_xz,
                                                 std::shared_ptr<const CfemExpression> expr_yx,
                                                 std::shared_ptr<const CfemExpression> expr_yy,
                                                 std::shared_ptr<const CfemExpression> expr_yz,
                                                 std::shared_ptr<const CfemExpression> expr_zx,
                                                 std::shared_ptr<const CfemExpression> expr_zy,
                                                 std::shared_ptr<const CfemExpression> expr_zz);

    /** 
     * @brief Build a symmetric tensor using 6 independent components (implicit Voigt notation). 
     */
    std::shared_ptr<const CfemExpression> symmetricTensor(std::shared_ptr<const CfemExpression> xx, 
                                                          std::shared_ptr<const CfemExpression> yy, 
                                                          std::shared_ptr<const CfemExpression> zz,
                                                          std::shared_ptr<const CfemExpression> yz, 
                                                          std::shared_ptr<const CfemExpression> xz, 
                                                          std::shared_ptr<const CfemExpression> xy);

    /** 
     * @brief Create a diagonal tensor using scalar expressions for the main diagonal. 
     */
    std::shared_ptr<const CfemExpression> diag(std::shared_ptr<const CfemExpression> a,
                                               std::shared_ptr<const CfemExpression> b,
                                               std::shared_ptr<const CfemExpression> c);
    
    /** 
     * @brief Create a constant diagonal tensor. 
     */
    std::shared_ptr<const CfemExpression> diag(CfemReal a, CfemReal b, CfemReal c);

    /** 
     * @brief Returns the Identity tensor expression. 
     */
    std::shared_ptr<const CfemExpression> identity();

    /** 
     * @brief Generates a constant vector expression.
     */
    std::shared_ptr<const CfemExpression> constant(const TensorO2& constantTensor);
    /** @} */

    /**
     * @name Differential Calculus and Operators
     * @{
     */
    
    /**
     * @brief Compute the gradient of a CFEM expression.
     * @param e Input expression. Scalar (dim 1) -> Gradient (dim 3). Vector (dim 3) -> Jacobian (dim 9).
     * @return Pointer to the differentiated expression.
     */
    std::shared_ptr<const CfemExpression> grad(const std::shared_ptr<const CfemExpression> &e);

    /**
     * @brief Compute the divergence of a vector or tensor field.
     * @param e Vector expression (dim 3) or tensor expression (dim 9).
     * @return Divergence (dim 1 for a vector, dim 3 for a tensor).
     */
    std::shared_ptr<const CfemExpression> div(const std::shared_ptr<const CfemExpression> &e);

    /**
     * @brief Compute the magnitude squared of a CFEM expression.
     * @param a Scalar or vector expression.
     * @return Pointer to the squared magnitude expression (scalar).
     */
    std::shared_ptr<const CfemExpression> mag_sqr(const std::shared_ptr<const CfemExpression> a);
    /** @} */

    /**
     * @name Variables and FE Trial/Test Functions
     * @{
     */
    
    /**
     * @brief Defines a physical or state variable (like space or time).
     * @param name Name of the variable for the symbolic engine.
     * @param dimIndex Associated spatial dimension index.
     */
    std::shared_ptr<const CfemExpression> variable(const std::string name, CfemInt dimIndex);

    /**
     * @brief Generate a Trial Function for variational formulations.
     * @param V The FE space the field belongs to.
     * @param name Symbolic variable representation (defaults to "u").
     */
    std::shared_ptr<FEField> Trial(std::shared_ptr<FESpace> V, const std::string& name = "u");
    
    /**
     * @brief Generate a Test Function for variational formulations.
     * @param V The associated FE space.
     * @param name Symbolic variable representation (defaults to "v").
     */
    std::shared_ptr<const TestFunction>  Test (std::shared_ptr<FESpace> V, const std::string& name = "v");
    /** @} */
}

namespace cfem::sym::utils {
    std::shared_ptr<const CfemExpression> extractComponent(std::shared_ptr<const CfemExpression> expr, int idx);

    template <typename Func>
    std::shared_ptr<const CfemExpression> foldUnary(std::shared_ptr<const CfemExpression>   u, 
                                                    MathOp                                  op, 
                                                    Func                                    func); 
    template <typename Func>
    std::shared_ptr<const CfemExpression> foldBinary(std::shared_ptr<const CfemExpression> lhs, 
                                                     std::shared_ptr<const CfemExpression> rhs, 
                                                     MathOp                                op, 
                                                     CfemInt                               outDim,
                                                     Func                                  func);
}

#include "MathOperators.tpp"

/**
 * @namespace cfem::sym::cfemvars
 * @brief Global pre-instantiated spacetime variables for a natural mathematical DSL.
 */
namespace cfem::sym::cfemvars{
    inline auto x = variable("x", 0);
    inline auto y = variable("y", 1);
    inline auto z = variable("z", 2);
    inline auto t = variable("t", 3);
}
