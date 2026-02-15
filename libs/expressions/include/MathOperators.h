#pragma once

#include <memory>
#include "CfemExpression.h"
#include "ConstantExpression.h"

#include "symengine/SymEngineSymbols.h"


namespace cfem::sym
{
    inline std::shared_ptr<const CfemExpression> transpose(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> dot(std::shared_ptr<const CfemExpression> lhs, 
                                                 std::shared_ptr<const CfemExpression> rhs);
    inline std::shared_ptr<const CfemExpression> inner(std::shared_ptr<const CfemExpression> lhs, 
                                                   std::shared_ptr<const CfemExpression> rhs);

    inline std::shared_ptr<const CfemExpression> cross(std::shared_ptr<const CfemExpression> lhs, 
                                                   std::shared_ptr<const CfemExpression> rhs);
    inline std::shared_ptr<const CfemExpression> atan2(std::shared_ptr<const CfemExpression> x, 
                                                       std::shared_ptr<const CfemExpression> y);
    inline std::shared_ptr<const CfemExpression> pow(std::shared_ptr<const CfemExpression> base, 
                                                       std::shared_ptr<const CfemExpression> expr);

    inline std::shared_ptr<const CfemExpression> abs(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> sqrt(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> cbrt(std::shared_ptr<const CfemExpression> u);

    inline std::shared_ptr<const CfemExpression> sign(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> heaviside(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> sin(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> cos(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> tan(std::shared_ptr<const CfemExpression> u);

    inline std::shared_ptr<const CfemExpression> cot(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> sec(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> csc(std::shared_ptr<const CfemExpression> u);

    inline std::shared_ptr<const CfemExpression> asin(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> acos(std::shared_ptr<const CfemExpression> u);

    inline std::shared_ptr<const CfemExpression> atan(std::shared_ptr<const CfemExpression> u);

    inline std::shared_ptr<const CfemExpression> acot(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> asec(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> acsc(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> sinh(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> cosh(std::shared_ptr<const CfemExpression> u);

    inline std::shared_ptr<const CfemExpression> tanh(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> asinh(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> acosh(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> atanh(std::shared_ptr<const CfemExpression> u);

    inline std::shared_ptr<const CfemExpression> exp(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> log(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> erf(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> erfc(std::shared_ptr<const CfemExpression> u);
    std::shared_ptr<const CfemExpression>        tgamma(std::shared_ptr<const CfemExpression> u);

    #if __cplusplus >= 201703L
    inline std::shared_ptr<const CfemExpression> zeta(std::shared_ptr<const CfemExpression> u);
    #endif

    // ==========================================
    // Opérateurs de Base
    // ==========================================
    inline std::shared_ptr<const CfemExpression> operator-(std::shared_ptr<const CfemExpression> u);
    inline std::shared_ptr<const CfemExpression> operator-(std::shared_ptr<const CfemExpression> lhs, 
                                                           std::shared_ptr<const CfemExpression> rhs);

    inline std::shared_ptr<const CfemExpression> operator+(std::shared_ptr<const CfemExpression> lhs, 
                                                           std::shared_ptr<const CfemExpression> rhs);

    inline std::shared_ptr<const CfemExpression> operator*(std::shared_ptr<const CfemExpression> lhs, 
                                                        std::shared_ptr<const CfemExpression> rhs);
    inline std::shared_ptr<const CfemExpression> operator/(std::shared_ptr<const CfemExpression> lhs, 
                                                           std::shared_ptr<const CfemExpression> rhs);
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
     * @brief Wrap a scalar lambda as a CfemExpression.
     * @param f Lambda function returning a CfemReal value from a 3D point (Vector3D).
     * @param gradient (Optional) An expression representing the analytical gradient of f. 
     *                  Must have a dimension of 3 (Vector).
     * @return A shared pointer to a CfemExpression of dimension 1.
     * @note If 'gradient' is not provided, any attempt to call grad() on this expression 
     *       (e.g., during H1-seminorm computation) will raise a runtime exception.
     * @throw CFEM_ASSERT if gradient is provided and its dimension is not 3.
     */
    std::shared_ptr<const CfemExpression> scalar(
                                            const std::function<CfemReal(const Vector3D &)> &f, 
                                            std::shared_ptr<const CfemExpression> gradient = nullptr);

    /**
     * @brief Wrap a vector lambda as a CfemExpression.
     * @param f Lambda function returning a Vector3D from a 3D point.
     * @param gradient (Optional) An expression representing the analytical Jacobian of f. 
     *                  Must have a dimension of 9 (3x3 Tensor).
     * @return A shared pointer to a CfemExpression of dimension 3.
     * @note If 'gradient' is not provided, operations requiring the derivative 
     *       (H1 semi-norm, divergence, etc.) will raise a runtime exception.
     * @throw CFEM_ASSERT if gradient is provided and its dimension is not 9.
     */
    std::shared_ptr<const CfemExpression> vector(
                                            const std::function<Vector3D(const Vector3D &)> &f,
                                            std::shared_ptr<const CfemExpression> gradient = nullptr);
    
    /**
     * @brief Wrap a 2nd-order tensor lambda as a CfemExpression.
     * @param f Lambda returning a 9-component array (3x3 tensor) from a 3D point.
     * @param gradient (Optional) Expression representing the 3rd-order gradient tensor.
     *         Must have a dimension of 27.
     * @return A shared pointer to a CfemExpression of dimension 9.
     * @note If 'gradient' is not provided, operations requiring the derivative 
     *       (H1 semi-norm, divergence, etc.) will raise a runtime exception.
     * @throw CFEM_ASSERT if gradient is provided and its dimension is not 27.
     */
    std::shared_ptr<const CfemExpression> tensor(
                                            const std::function<TensorO2(const Vector3D &)> &f,
                                            std::shared_ptr<const CfemExpression> gradient = nullptr);

    /**
     * @brief Compute the gradient of a CFEM expression.
     * @param e Shared pointer to the input expression (scalar or vector)
     * @return Shared pointer to the gradient (vector or tensor)
     */
    std::shared_ptr<const CfemExpression> grad(const std::shared_ptr<const CfemExpression> &e);

    /**
     * @brief Compute the divergence of a CFEM expression.
     * @param e Shared pointer to the input vector/tensor expression
     * @return Shared pointer to the scalar divergence
     */
    std::shared_ptr<const CfemExpression> div(const std::shared_ptr<const CfemExpression> &e);

    /**
     * @brief Compute the magnitude squared of a CFEM expression.
     * @param a Shared pointer to a scalar or vector expression
     * @return Shared pointer to the squared magnitude expression
     */
    std::shared_ptr<const CfemExpression> mag_sqr(const std::shared_ptr<const CfemExpression> a);

    // ---------------------------------------------------
    // SymEngine Expression Helpers
    // ---------------------------------------------------
    /**
     * @brief Convert a SymEngine scalar expression into a CFEM expression.
     * @param ex The scalar symbolic expression
     * @return Shared pointer to a CfemExpression (ScalarSymbolicExpression)
     */
    std::shared_ptr<const CfemExpression> scalar(const Expression &ex);

    /**
     * @brief Convert three SymEngine expressions into a vector CFEM expression.
     * @param ex Expression for x-component
     * @param ey Expression for y-component
     * @param ez Expression for z-component
     * @return Shared pointer to a CfemExpression (VectorSymbolicExpression)
     */
    std::shared_ptr<const CfemExpression> vector(const Expression &ex,
                                                 const Expression &ey,
                                                 const Expression &ez);

}

#include "MathOperators.tpp"
