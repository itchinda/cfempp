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

#ifndef CFEM_VECTOR3D
#define CFEM_VECTOR3D

#include "cfemutils.h"
#include <cmath>
#include <sstream>

/**
 * @struct Vector3D
 * @brief Represents a vector or point in 3D Euclidean space.
 * * This struct provides the basic linear algebra and vector calculus tools 
 * required for mesh transformations, Jacobian calculations, and normal extraction.
 */
struct Vector3D {
    CfemReal x = 0.0; ///< X-coordinate
    CfemReal y = 0.0; ///< Y-coordinate
    CfemReal z = 0.0; ///< Z-coordinate

    static constexpr int size = 3;

    // --- Basic Arithmetic Operators ---

    /** @brief Vector addition. */
    Vector3D operator+(const Vector3D& other) const {
        return {x + other.x, y + other.y, z + other.z};
    }

    /** @brief Vector subtraction. */
    Vector3D operator-(const Vector3D& other) const {
        return {x - other.x, y - other.y, z - other.z};
    }

    /** @brief Right-side scalar multiplication. */
    Vector3D operator*(CfemReal scalar) const {
        return {x * scalar, y * scalar, z * scalar};
    }

    /** @brief Scalar division with zero-tolerance check. */
    Vector3D operator/(CfemReal scalar) const {
        if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
            CFEM_ERROR("Division by zero in Vector3D operator/");
        return {x / scalar, y / scalar, z / scalar};
    }

    // --- Compound Assignment Operators ---

    Vector3D& operator+=(const Vector3D& other) {
        x += other.x; y += other.y; z += other.z;
        return *this;
    }

    Vector3D& operator-=(const Vector3D& other) {
        x -= other.x; y -= other.y; z -= other.z;
        return *this;
    }

    Vector3D& operator*=(CfemReal scalar) {
        x *= scalar; y *= scalar; z *= scalar;
        return *this;
    }

    Vector3D& operator/=(CfemReal scalar) {
        if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
            CFEM_ERROR("Division by zero in Vector3D operator/=");
        x /= scalar; y /= scalar; z /= scalar;
        return *this;
    }
    
    // --- Accessors ---

    /** @brief Read-only access by index (0=x, 1=y, 2=z). */
    CfemReal operator[](int i) const {
        if (i == 0) return x;
        if (i == 1) return y;
        if (i == 2) return z;
        CFEM_ERROR("Index out of bound in Vector3D::operator[] const.");
        return 0.0; 
    }

    /** @brief Read-write access by index (0=x, 1=y, 2=z). */
    CfemReal& operator[](int i) {
        if (i == 0) return x;
        if (i == 1) return y;
        if (i == 2) return z;
        CFEM_ERROR("Index out of bound in Vector3D::operator[].");
        return x; 
    }
    
    // --- Vector Calculus ---

    /** * @brief Computes the dot (inner) product.
     * @return a.x*b.x + a.y*b.y + a.z*b.z
     */
    CfemReal dot(const Vector3D& other) const {
        return x * other.x + y * other.y + z * other.z;
    }

    /** * @brief Computes the cross product. 
     * @return A vector perpendicular to both this and other.
     */
    Vector3D cross(const Vector3D& other) const {
        return {
            y * other.z - z * other.y,
            z * other.x - x * other.z,
            x * other.y - y * other.x
        };
    }
    
    /** @brief Returns the squared Euclidean norm (faster than norm()). */
    CfemReal normSquare() const {
        return x * x + y * y + z * z;
    }

    /** @brief Returns the Euclidean norm (length) of the vector. */
    CfemReal norm() const {
        return std::sqrt(normSquare());
    }

    /** * @brief Returns a unit vector in the same direction.
     * Returns {0,0,0} if the vector is near-null.
     */
    Vector3D normalized() const {
        CfemReal n = norm();
        if (n < cfem::ZERO_TOLERANCE) return {0, 0, 0};
        return (*this) / n;
    }

    const Vector3D& normalize() {
        CfemReal n = norm();
        if (n > cfem::ZERO_TOLERANCE) {
            x /= n; y /= n; z /= n;
        }
        return *this;
    }

    bool isClose(const Vector3D& other, CfemReal tol = cfem::ZERO_TOLERANCE) const {
        return (std::abs(x - other.x) < tol && 
                std::abs(y - other.y) < tol && 
                std::abs(z - other.z) < tol);
    }
};

// --- Non-member Helper Functions ---

/** @brief Left-side scalar multiplication (e.g., 2.0 * vector). */
inline Vector3D operator*(CfemReal scalar, const Vector3D& p) {
    return {p.x * scalar, p.y * scalar, p.z * scalar};
}

/** @brief Returns the Euclidean norm (length) of the vector. */
inline CfemReal norm(const Vector3D& v)  {
    return v.norm();
}

/** @brief Functional style dot product. */
inline CfemReal dot(const Vector3D& a, const Vector3D& b) {
    return a.dot(b);
}

/** @brief Functional style cross product. */
inline Vector3D cross(const Vector3D& a, const Vector3D& b) {
    return a.cross(b);
}

/** @brief Overloads the output stream operator for easy printing/logging.
 *  @param os Output stream (e.g., std::cout or a file stream).
 *  @param p The Vector3D to print.
 *  @return Reference to the stream for chaining.
 */
inline std::ostream& operator<<(std::ostream& os, const Vector3D& p) {
    return os << "(" << p.x << ", " << p.y << ", " << p.z << ")";
}


#endif // CFEM_VECTOR3D