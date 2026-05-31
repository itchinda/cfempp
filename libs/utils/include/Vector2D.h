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

/**
 * @file Vector2D.h
 * @author Itchinda
 * @brief Lightweight 3D vector class for Finite Element geometric operations.
 * @version 0.1
 * @date 2026-01-27
 *
 * @copyright Copyright (c) 2026
 *
 */

#pragma once

#include "cfemutils.h"
#include <cmath>
#include <ostream>
#include <algorithm>

namespace cfem
{

    /**
     * @struct Vector2D
     * @brief Represents a vector or point in 2D Euclidean space.
     * 
     * This struct provides the basic linear algebra and vector calculus tools 
     * required for mesh transformations, Jacobian calculations, and normal extraction.
     * It is highly optimized, utilizing `constexpr` for compile-time evaluation where possible.
     * 
     * @par Example Usage:
     * \code{.cpp}
     * cfem::Vector2D a(1.0, 0.0);
     * cfem::Vector2D b(0.0, 1.0);
     * cfem::CfemReal normal_z = a.cross(b).normalized(); // Returns (0, 0, 1)
     * // cfem::CfemReal cross_prod = a.cross(b); // Returns 1.0 (scalar representing the Z-component)
     * \endcode
     */
    struct Vector2D
    {
        CfemReal x = 0.0; ///< X-coordinate
        CfemReal y = 0.0; ///< Y-coordinate

        static constexpr CfemInt size = 2;

        constexpr Vector2D() noexcept = default;
        constexpr Vector2D( CfemReal x_, CfemReal y_) noexcept
            : x(x_), y(y_)
        {
        }

        // Basic Arithmetic Operators

        /** @brief Vector addition. */
        [[nodiscard]] constexpr Vector2D operator+(const Vector2D &other) const noexcept
        {
            return {x + other.x, y + other.y};
        }

        /** 
         * @brief Vector subtraction. 
         */
        [[nodiscard]] constexpr Vector2D operator-(const Vector2D &other) const noexcept
        {
            return {x - other.x, y - other.y};
        }

        /** @brief Right-side scalar multiplication. */
        [[nodiscard]] constexpr Vector2D operator*(CfemReal scalar) const noexcept
        {
            return {x * scalar, y * scalar};
        }

        [[nodiscard]] constexpr Vector2D operator-() const noexcept
        {
            return {-x, -y};
        }

        /**
         * @brief Right-side scalar division.
         * @param scalar The denominator.
         * @return A new vector divided by the scalar. 
         * @note If the scalar is near zero (below cfem::ZERO_TOLERANCE), a CFEM_ERROR is logged, 
         * but the division proceeds (which may result in Inf or NaN values).
         */
        [[nodiscard]] inline Vector2D operator/(CfemReal scalar) const
        {
            CFEM_ASSERT(std::abs(scalar) > ZERO_TOLERANCE, "Division by zero in Vector2D::operator/");
            
            CfemReal inv = 1.0 / scalar;
            return {x * inv, y * inv};
        }

        // --- Compound Assignment Operators ---

        constexpr Vector2D &operator+=(const Vector2D &other) noexcept
        {
            x += other.x;
            y += other.y;
            return *this;
        }

        constexpr Vector2D &operator-=(const Vector2D &other) noexcept
        {
            x -= other.x;
            y -= other.y;
            return *this;
        }

        /**
         * @brief In place multiplication by a scalar
         * 
         * @param scalar 
         * @return Vector2D& 
         */
        constexpr Vector2D &operator*=(CfemReal scalar) noexcept
        {
            x *= scalar;
            y *= scalar;
            return *this;
        }

        /**
         * @brief In-place Scalar division.
         * @param scalar The denominator.
         * @note If the scalar is near zero (below cfem::ZERO_TOLERANCE), a CFEM_ERROR is logged, 
         * but the division proceeds (which may result in Inf or NaN values).
         */
        constexpr Vector2D &operator/=(CfemReal scalar)
        {
            CFEM_ASSERT(std::abs(scalar) > ZERO_TOLERANCE, "Division by zero in Vect::operator/=");

            CfemReal inv = 1.0 / scalar;
            x *= inv;
            y *= inv;
            return *this;
        }

        /**
         * @brief Evaluates strict bitwise equality.
         * @param other The vector to compare with.
         * @return True if all components match exactly, false otherwise.
         * @warning Do NOT use this for geometric or physical checks due to floating-point 
         * inaccuracies. Use isApprox() instead for tolerance-based equality.
         */
        [[nodiscard]] constexpr bool operator==(const Vector2D& other) const noexcept {
            return x == other.x && y == other.y;
        }

        [[nodiscard]] bool operator!=(const Vector2D &other) const noexcept
        {
            return !(*this == other);
        }

        /**
         * @brief Checks if two vectors are geometrically identical up to a numerical tolerance.
         * 
         * This method relies on the global NumericPolicy to safely compare floating-point 
         * coordinates. It dynamically computes a local scale for each axis to prevent 
         * floating-point traps when dealing with very large or very small coordinates.
         * 
         * @param other The vector to compare with.
         * @return True if the vectors are identical within the defined tolerance, false otherwise.
         */
        [[nodiscard]] bool isApprox(const Vector2D& other) const noexcept {
            CfemReal scale_x = std::max(std::abs(x), std::abs(other.x));
            CfemReal scale_y = std::max(std::abs(y), std::abs(other.y));

            return NumericPolicy::areEqual(x, other.x, scale_x) &&
                   NumericPolicy::areEqual(y, other.y, scale_y);
        }

        // Accessors
        /** 
         * @brief Read-only access by index (0=x, 1=y, 2=z). 
         */
        [[nodiscard]] inline CfemReal operator[](CfemInt i) const {
            CFEM_ASSERT(0 <= i && i < 2, "Index out of bound in Vector2D::operator[] const.");
            return (i == 0) ? x : y;
        }

        /** 
         * @brief Read-write access by index (0=x, 1=y, 2=z). 
         */
        [[nodiscard]] inline CfemReal& operator[](CfemInt i) {
            CFEM_ASSERT(0 <= i && i < 2, "Index out of bound in Vector2D::operator[].");
            return (i == 0) ? x : y;
        }

        [[nodiscard]] constexpr const CfemReal* data() const noexcept
        {
            return &x;
        }

        [[nodiscard]] constexpr CfemReal* data() noexcept
        {
            return &x;
        }

        // Vector Calculus

        /** 
         * @brief Computes the dot (inner) product.
         * @return a.x*b.x + a.y*b.y
         */
        [[nodiscard]] constexpr CfemReal dot(const Vector2D &other) const noexcept
        {
            return x * other.x + y * other.y;
        }

        /** 
         * @brief Computes the cross product.
         * @return A scalar representing the z-component of a vector perpendicular to both vectors
         */
        [[nodiscard]] constexpr CfemReal cross(const Vector2D &other) const noexcept
        {
            return x * other.y - y * other.x;
        }

        /** @brief Returns the squared Euclidean norm (faster than norm()). */
        [[nodiscard]] constexpr CfemReal normSquare() const noexcept
        {
            return x * x + y * y;
        }

        /** @brief Returns the Euclidean norm (length) of the vector. */
        [[nodiscard]] CfemReal norm() const noexcept
        {
            return std::sqrt(normSquare());
        }

        /** 
         * @brief Returns a unit vector in the same direction.
         * Returns {0,0,0} if the vector is near-null.
         */
        [[nodiscard]] Vector2D normalized() const noexcept
        {
            CfemReal n = norm();
            if (n < cfem::ZERO_TOLERANCE)
                return {0, 0};
            return (*this) / n;
        }

        Vector2D &normalize() noexcept
        {
            CfemReal n = norm();
            if (n > cfem::ZERO_TOLERANCE)
            {
                *this /= n;
            }
            return *this;
        }

        /** 
         * @brief Returns the squared distance to another point. (Faster, no sqrt)
         */
        [[nodiscard]] constexpr CfemReal distanceSquaredTo(const Vector2D& other) const noexcept {
            return (*this - other).normSquare();
        }

        /** 
         * @brief Returns the Euclidean distance to another point. 
         */
        [[nodiscard]] CfemReal distanceTo(const Vector2D& other) const noexcept {
            return std::sqrt(distanceSquaredTo(other));
        }

        [[nodiscard]] constexpr bool isClose(const Vector2D &other, CfemReal tol = cfem::ZERO_TOLERANCE) const noexcept
        {
            return (std::abs(x - other.x) < tol &&
                    std::abs(y - other.y) < tol);
        }

        /** 
         * @brief Computes the angle (in radians) between this vector and another.
         * Safe against floating point inaccuracies.
         */
        [[nodiscard]] CfemReal angleWith(const Vector2D& other) const noexcept {
            return std::acos(cosAngleWith(other));
        }

        /**
         * @brief Computes the cosine of the angle between this vector and another.
         * @param other The second vector.
         * @return The cosine of the angle in the range [-1.0, 1.0].
         * @note Faster than angleWith() and numerically stable.
         */
        [[nodiscard]] CfemReal cosAngleWith(const Vector2D& other) const noexcept {
            CfemReal denom = this->norm() * other.norm();
            
            if (denom < cfem::ZERO_TOLERANCE) return 0.0;
            
            return std::clamp(this->dot(other) / denom, -1.0, 1.0); // std::clamp to avoid floating point errors (the result should stay between -1 and 1)
        }

        /** 
         * @brief Computes the orthogonal projection of this vector onto a target direction.
         * @param target The vector defining the projection direction.
         * @return The projected Vector2D. Returns {0,0,0} if the target is a null vector.
         */
        [[nodiscard]] constexpr Vector2D projectOnto(const Vector2D& target) const {
            CfemReal targetNormSq = target.normSquare();
            if (targetNormSq < cfem::ZERO_TOLERANCE) return {0, 0};
            return target * (this->dot(target) / targetNormSq);
        }

        /** 
         * @brief Returns the component of this vector perpendicular to 'target'. 
         */
        [[nodiscard]] constexpr Vector2D perpendicularTo(const Vector2D& target) const {
            return *this - projectOnto(target);
        }

        /** 
         * @brief Returns a vector containing the minimum components of both vectors. 
         */
        [[nodiscard]] constexpr Vector2D cwiseMin(const Vector2D& other) const {
            return {std::min(x, other.x), std::min(y, other.y)};
        }

        /** 
         * @brief Returns a vector containing the maximum components of both vectors. 
         */
        [[nodiscard]] constexpr Vector2D cwiseMax(const Vector2D& other) const {
            return {std::max(x, other.x), std::max(y, other.y)};
        }
        
        /** 
         * @brief Returns the vector resulting from component-wise multiplication. 
         */
        [[nodiscard]] constexpr Vector2D cwiseProduct(const Vector2D& other) const {
            return {x * other.x, y * other.y};
        }
    };

    // Non-member Helper Functions

    /** @brief Left-side scalar multiplication (e.g., 2.0 * vector). */
    [[nodiscard]] constexpr Vector2D operator*(CfemReal scalar, const Vector2D &p) noexcept
    {
        return {p.x * scalar, p.y * scalar};
    }

    /** 
     * @brief Returns the Euclidean norm (length) of the vector. 
     */
    [[nodiscard]] inline CfemReal norm(const Vector2D &v) noexcept
    {
        return v.norm();
    }

    /** 
     * @brief Functional style dot product. 
     */
    [[nodiscard]] constexpr CfemReal dot(const Vector2D &a, const Vector2D &b) noexcept
    {
        return a.dot(b);
    }

    /** @brief Functional style cross product. */
    [[nodiscard]] constexpr CfemReal cross(const Vector2D &a, const Vector2D &b) noexcept
    {
        return a.cross(b);
    }

    /** @brief Overloads the output stream operator for easy printing/logging.
     *  @param os Output stream (e.g., std::cout or a file stream).
     *  @param p The Vector2D to print.
     *  @return Reference to the stream for chaining.
     */
    inline std::ostream &operator<<(std::ostream &os, const Vector2D &p)
    {
        return os << "(" << p.x << ", " << p.y << ", " << ")";
    }

} // namespace
