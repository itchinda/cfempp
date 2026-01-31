/**
 * @file StaticVector.h
 * @author Itchinda
 * @brief Fixed-size stack-allocated vector for Finite Element local operations.
 * @version 0.1
 * @date 2026-01-27
 */

#ifndef CFEM_STATIC_VECTOR
#define CFEM_STATIC_VECTOR

#include "cfemutils.h"
#include "Vector3D.h"
#include <iostream>
#include <cmath>
#include <initializer_list>

/**
 * @brief A fixed-size column vector of size N.
 * @tparam T The numeric type (e.g., CfemReal).
 * @tparam N Dimension of the vector.
 * * This class is intended for small, element-level algebraic operations 
 * (local residuals, shape function values, etc.) where stack allocation
 * and compile-time size checking are critical for performance.
 */
template <typename T, int N>
struct StaticVector {
    /** @brief Contiguous data storage on the stack. */
    T data[N] = {0};

    // --- Constructors ---

    /** @brief Default constructor (initializes to zero). */
    StaticVector() = default;

    /** @brief Initializer list constructor for syntax: StaticVector v = {1, 2, 3}; */
    StaticVector(std::initializer_list<T> list) {
        CFEM_ASSERT(list.size() == N);
        int i = 0;
        for (auto it = list.begin(); it != list.end() && i < N; ++it, ++i) {
            data[i] = *it;
        }
    }

    /** * @brief Conversion Constructor from Vector3D. 
     * Only available when N == 3 via static_assert.
     */
    StaticVector(const Vector3D& v) {
        static_assert(N == 3, "StaticVector size must be 3 to initialize from Vector3D.");
        data[0] = v.x;
        data[1] = v.y;
        data[2] = v.z;
    }

    // --- Type Conversion (Casting) ---

    /** * @brief Implicit conversion to Vector3D.
     * Allows seamless transition between algebraic and geometric logic.
     */
    operator Vector3D() const {
        static_assert(N == 3, "StaticVector size must be 3 to convert to Vector3D.");
        return { data[0], data[1], data[2] };
    }

    // --- Accessors ---

    /** @brief Access element with [] notation. */
    inline T& operator[](int i) {
        CFEM_ASSERT(i >= 0 && i < N);
        return data[i];
    }

    /** @brief Read-only access element with [] notation. */
    inline const T& operator[](int i) const {
        CFEM_ASSERT(i >= 0 && i < N);
        return data[i];
    }

    /** @brief Returns the fixed size of the vector. */
    static constexpr int size() noexcept { return N; }

    // --- Basic Math ---

    /** @brief Resets all components to zero. */
    void setZero() {
        for (int i = 0; i < N; ++i) data[i] = static_cast<T>(0);
    }

    // --- Arithmetic Operators ---

    StaticVector operator+(const StaticVector& other) const {
        StaticVector result;
        for (int i = 0; i < N; ++i) result.data[i] = data[i] + other.data[i];
        return result;
    }

    StaticVector operator-(const StaticVector& other) const {
        StaticVector result;
        for (int i = 0; i < N; ++i) result.data[i] = data[i] - other.data[i];
        return result;
    }

    StaticVector operator*(T scalar) const {
        StaticVector result;
        for (int i = 0; i < N; ++i) result.data[i] = data[i] * scalar;
        return result;
    }

    StaticVector operator/(T scalar) const {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        StaticVector result;
        for (int i = 0; i < N; ++i) result.data[i] = data[i] / scalar;
        return result;
    }

    // --- Compound Operators ---

    StaticVector& operator+=(const StaticVector& other) {
        for (int i = 0; i < N; ++i) data[i] += other.data[i];
        return *this;
    }

    StaticVector& operator-=(const StaticVector& other) {
        for (int i = 0; i < N; ++i) data[i] -= other.data[i];
        return *this;
    }

    StaticVector& operator*=(T scalar) {
        for (int i = 0; i < N; ++i) data[i] *= scalar;
        return *this;
    }

    StaticVector& operator/=(T scalar) {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        for (int i = 0; i < N; ++i) data[i] /= scalar;
        return *this;
    }

    // --- Linear Algebra Tools ---

    /** @brief Computes the dot (inner) product. */
    T dot(const StaticVector<T, N>& other) const {
        T result = 0;
        for (int i = 0; i < N; ++i) result += data[i] * other.data[i];
        return result;
    }

    /** @brief Returns squared Euclidean norm. */
    T normSquare() const {
        return dot(*this);
    }

    /** @brief Returns Euclidean norm (length). */
    T norm() const {
        return std::sqrt(normSquare());
    }
};

// --- Non-member Helper Functions ---

/** @brief Left-side scalar multiplication (e.g., scalar * vector). */
template <typename T, int N>
inline StaticVector<T, N> operator*(T scalar, const StaticVector<T, N>& v) {
    return v * scalar;
}

/** @brief Formats the vector for standard output. */
template <typename T, int N>
inline std::ostream& operator<<(std::ostream& os, const StaticVector<T, N>& v) {
    os << "(";
    for (int i = 0; i < N; ++i) {
        os << v[i] << (i == N - 1 ? "" : ", ");
    }
    os << ")";
    return os;
}

#endif // CFEM_STATIC_VECTOR