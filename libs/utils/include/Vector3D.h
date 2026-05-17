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
 * @file Vector3D.h
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
     * @struct Vector3D
     * @brief Represents a vector or point in 3D Euclidean space.
     * 
     * This struct provides the basic linear algebra and vector calculus tools 
     * required for mesh transformations, Jacobian calculations, and normal extraction.
     * It is highly optimized, utilizing `constexpr` for compile-time evaluation where possible.
     * 
     * @par Example Usage:
     * \code{.cpp}
     * cfem::Vector3D a(1.0, 0.0, 0.0);
     * cfem::Vector3D b(0.0, 1.0, 0.0);
     * cfem::Vector3D normal = a.cross(b).normalized(); // Returns (0, 0, 1)
     * CfemReal area = a.cross(b).norm();               // Returns 1.0
     * \endcode
     */
    struct Vector3D
    {
        CfemReal x = 0.0; ///< X-coordinate
        CfemReal y = 0.0; ///< Y-coordinate
        CfemReal z = 0.0; ///< Z-coordinate

        static constexpr CfemInt size = 3;

        constexpr Vector3D() noexcept = default;

        constexpr Vector3D( CfemReal x_, CfemReal y_, CfemReal z_) noexcept
            : x(x_), y(y_), z(z_) 
        {
        }

        // Basic Arithmetic Operators

        /** @brief Vector addition. */
        [[nodiscard]] constexpr Vector3D operator+(const Vector3D &other) const noexcept
        {
            return {x + other.x, y + other.y, z + other.z};
        }

        /** 
         * @brief Vector subtraction. 
         */
        [[nodiscard]] constexpr Vector3D operator-(const Vector3D &other) const noexcept
        {
            return {x - other.x, y - other.y, z - other.z};
        }

        /** @brief Right-side scalar multiplication. */
        [[nodiscard]] constexpr Vector3D operator*(CfemReal scalar) const noexcept
        {
            return {x * scalar, y * scalar, z * scalar};
        }

        /**
         * @brief Right-side scalar division.
         * @param scalar The denominator.
         * @return A new vector divided by the scalar. 
         * @note If the scalar is near zero (below cfem::ZERO_TOLERANCE), a CFEM_ERROR is logged, 
         * but the division proceeds (which may result in Inf or NaN values).
         */
        [[nodiscard]] inline Vector3D operator/(CfemReal scalar) const
        {
            CFEM_ASSERT(std::abs(scalar) > ZERO_TOLERANCE, "Division by zero in Vector3D::operator/");
            
            CfemReal inv = 1.0 / scalar;
            return {x * inv, y * inv, z * inv};
        }

        [[nodiscard]] constexpr Vector3D operator-() const noexcept
        {
            return {-x, -y, -z};
        }

        // Compound Assignment Operators

        constexpr Vector3D &operator+=(const Vector3D &other) noexcept
        {
            x += other.x;
            y += other.y;
            z += other.z;
            return *this;
        }

        constexpr Vector3D &operator-=(const Vector3D &other) noexcept
        {
            x -= other.x;
            y -= other.y;
            z -= other.z;
            return *this;
        }

        /**
         * @brief In place multiplication by a scalar
         * 
         * @param scalar 
         * @return Vector3D& 
         */
        constexpr Vector3D &operator*=(CfemReal scalar) noexcept
        {
            x *= scalar;
            y *= scalar;
            z *= scalar;
            return *this;
        }

        /**
         * @brief In-place Scalar division.
         * @param scalar The denominator.
         * @note If the scalar is near zero (below cfem::ZERO_TOLERANCE), a CFEM_ERROR is logged, 
         * but the division proceeds (which may result in Inf or NaN values).
         */
        inline Vector3D &operator/=(CfemReal scalar)
        {
            CFEM_ASSERT(std::abs(scalar) > ZERO_TOLERANCE, "Division by zero in Vector3D::operator/=");
            
            CfemReal inv = 1.0 / scalar;
            x *= inv;
            y *= inv;
            z *= inv;
            return *this;
        }

        /**
         * @brief Evaluates strict bitwise equality.
         * @param other The vector to compare with.
         * @return True if all components match exactly, false otherwise.
         * @warning Do NOT use this for geometric or physical checks due to floating-point 
         * inaccuracies. Use isApprox() instead for tolerance-based equality.
         */
        [[nodiscard]] constexpr bool operator==(const Vector3D& other) const noexcept {
            return x == other.x && y == other.y && z == other.z;
        }

        [[nodiscard]] constexpr bool operator!=(const Vector3D &other) const noexcept
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
        [[nodiscard]] constexpr bool isApprox(const Vector3D& other) const noexcept {
            CfemReal scale_x = std::max(std::abs(x), std::abs(other.x));
            CfemReal scale_y = std::max(std::abs(y), std::abs(other.y));
            CfemReal scale_z = std::max(std::abs(z), std::abs(other.z));

            return NumericPolicy::areEqual(x, other.x, scale_x) &&
                   NumericPolicy::areEqual(y, other.y, scale_y) &&
                   NumericPolicy::areEqual(z, other.z, scale_z);
        }

        // Accessors

        /** @brief Read-only access by index (0=x, 1=y, 2=z). */
        [[nodiscard]] constexpr CfemReal operator[](CfemInt i) const
        {
            CFEM_ASSERT( 0 <= i && i < 3, "Index out of bound in Vector3D::operator[] const.");
            return (i == 0) ? x : ((i == 1) ? y : z);
        }

        /** 
         * @brief Read-write access by index (0=x, 1=y, 2=z). 
         */
        [[nodiscard]] constexpr CfemReal &operator[](CfemInt i)
        {
            CFEM_ASSERT( 0 <= i && i < 3, "Index out of bound in Vector3D::operator[].");
            return (i == 0) ? x : ((i == 1) ? y : z);
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

        /** * @brief Computes the dot (inner) product.
         * @return a.x*b.x + a.y*b.y + a.z*b.z
         */
        [[nodiscard]] constexpr CfemReal dot(const Vector3D &other) const noexcept
        {
            return x * other.x + y * other.y + z * other.z;
        }

        /** * @brief Computes the cross product.
         * @return A vector perpendicular to both this and other.
         */
        [[nodiscard]] constexpr Vector3D cross(const Vector3D &other) const noexcept
        {
            return {
                y * other.z - z * other.y,
                z * other.x - x * other.z,
                x * other.y - y * other.x};
        }

        /** 
         * @brief Computes the scalar triple product: this . (b x c).
         * Extremely useful for computing 3D element volumes and Jacobian determinants.
         */
        [[nodiscard]] constexpr CfemReal scalarTripleProduct(const Vector3D& b, const Vector3D& c) const {
            return this->dot(b.cross(c));
        }

        /** @brief Returns the squared Euclidean norm (faster than norm()). */
        [[nodiscard]] constexpr CfemReal normSquare() const noexcept
        {
            return x * x + y * y + z * z;
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
        [[nodiscard]] Vector3D normalized() const noexcept
        {
            CfemReal n = norm();
            if (n < cfem::ZERO_TOLERANCE)
                return {0, 0, 0};
            return (*this) / n;
        }

        Vector3D &normalize() noexcept
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
        [[nodiscard]] constexpr CfemReal distanceSquaredTo(const Vector3D& other) const {
            return (*this - other).normSquare();
        }

        /** 
         * @brief Returns the Euclidean distance to another point. 
         */
        [[nodiscard]] CfemReal distanceTo(const Vector3D& other) const {
            return std::sqrt(distanceSquaredTo(other));
        }

        [[nodiscard]] constexpr bool isClose(const Vector3D &other, CfemReal tol = cfem::ZERO_TOLERANCE) const noexcept
        {
            return (std::abs(x - other.x) < tol &&
                    std::abs(y - other.y) < tol &&
                    std::abs(z - other.z) < tol);
        }

        /** 
         * @brief Computes the angle (in radians) between this vector and another.
         * Safe against floating point inaccuracies.
         */
        [[nodiscard]] CfemReal angleWith(const Vector3D& other) const {
            return std::acos(cosAngleWith(other));
        }

        /**
         * @brief Computes the cosine of the angle between this vector and another.
         * @param other The second vector.
         * @return The cosine of the angle in the range [-1.0, 1.0].
         * @note Faster than angleWith() and numerically stable.
         */
        [[nodiscard]] CfemReal cosAngleWith(const Vector3D& other) const {
            CfemReal denom = this->norm() * other.norm();
            
            if (denom < cfem::ZERO_TOLERANCE) return 0.0;
            
            return std::clamp(this->dot(other) / denom, -1.0, 1.0); // std::clamp to avoid floating point errors (the result should stay between -1 and 1)
        }

        /** 
         * @brief Computes the orthogonal projection of this vector onto a target direction.
         * @param target The vector defining the projection direction.
         * @return The projected Vector3D. Returns {0,0,0} if the target is a null vector.
         */
        [[nodiscard]] constexpr Vector3D projectOnto(const Vector3D& target) const {
            CfemReal targetNormSq = target.normSquare();
            if (targetNormSq < cfem::ZERO_TOLERANCE) return {0, 0, 0};
            return target * (this->dot(target) / targetNormSq);
        }

        /** 
         * @brief Returns the component of this vector perpendicular to 'target'. 
         */
        [[nodiscard]] constexpr Vector3D perpendicularTo(const Vector3D& target) const {
            return *this - projectOnto(target);
        }

        /** 
         * @brief Returns a vector containing the minimum components of both vectors. 
         */
        [[nodiscard]] constexpr Vector3D cwiseMin(const Vector3D& other) const {
            return {std::min(x, other.x), std::min(y, other.y), std::min(z, other.z)};
        }

        /** 
         * @brief Returns a vector containing the maximum components of both vectors. 
         */
        [[nodiscard]] constexpr Vector3D cwiseMax(const Vector3D& other) const {
            return {std::max(x, other.x), std::max(y, other.y), std::max(z, other.z)};
        }
        
        /** 
         * @brief Returns the vector resulting from component-wise multiplication. 
         */
        [[nodiscard]] constexpr Vector3D cwiseProduct(const Vector3D& other) const {
            return {x * other.x, y * other.y, z * other.z};
        }
    };

    // Non-member Helper Functions

    /** @brief Left-side scalar multiplication (e.g., 2.0 * vector). */
    [[nodiscard]] constexpr Vector3D operator*(CfemReal scalar, const Vector3D &p) noexcept
    {
        return {p.x * scalar, p.y * scalar, p.z * scalar};
    }

    /** @brief Returns the Euclidean norm (length) of the vector. */
    [[nodiscard]] inline CfemReal norm(const Vector3D &v) noexcept
    {
        return v.norm();
    }

    /** @brief Functional style dot product. */
    [[nodiscard]] constexpr CfemReal dot(const Vector3D &a, const Vector3D &b) noexcept
    {
        return a.dot(b);
    }

    /** @brief Functional style cross product. */
    [[nodiscard]] constexpr Vector3D cross(const Vector3D &a, const Vector3D &b) noexcept
    {
        return a.cross(b);
    }

    /** @brief Overloads the output stream operator for easy printing/logging.
     *  @param os Output stream (e.g., std::cout or a file stream).
     *  @param p The Vector3D to print.
     *  @return Reference to the stream for chaining.
     */
    inline std::ostream &operator<<(std::ostream &os, const Vector3D &p)
    {
        return os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
    }

} // namespace