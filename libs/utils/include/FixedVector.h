/*!
 * @file FixedVector.h
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


namespace cfem{


/*!
 * @struct FixedVector
 * @brief A fixed-size column vector of size N.
 * @tparam T The numeric type (e.g., CfemReal, double).
 * @tparam N Dimension of the vector (must be > 0).
 * @note This class is intended for small, element-level algebraic operations 
 *       (local residuals, shape function values, etc.) where stack allocation
 *       and compile-time size checking are critical for performance.
 */
template <typename T, int N>
struct FixedVector {
    
    /*! 
     * @brief Contiguous data storage on the stack.
     */
    T data[N] = {0};

    /*! 
     * @brief Default constructor. Initializes all elements to zero.
     */
    FixedVector() = default;

    /*!
     * @brief Initializer list constructor.
     * @param list A list of values (e.g., {1.0, 2.0, 3.0}).
     * @throw CFEM_ASSERT if list.size() != N.
     */
    FixedVector(std::initializer_list<T> list) {
        CFEM_ASSERT(list.size() == N);
        int i = 0;
        for (auto it = list.begin(); it != list.end() && i < N; ++it, ++i) {
            data[i] = *it;
        }
    }

    /*! 
     * @brief Conversion Constructor from Vector3D. 
     * @note Only available when N == 3.
     * @param v Geometric Vector3D to convert from.
     */
    FixedVector(const Vector3D& v) {
        static_assert(N == 3, "FixedVector size must be 3 to initialize from Vector3D.");
        data[0] = v.x;
        data[1] = v.y;
        data[2] = v.z;
    }

    /*! 
     * @brief Implicit conversion to Vector3D.
     * @note Only available when N == 3.
     * @return A Vector3D representation of the first three components.
     */
    operator Vector3D() const {
        static_assert(N == 3, "FixedVector size must be 3 to convert to Vector3D.");
        return { data[0], data[1], data[2] };
    }

    // --- Accessors ---

    /*! 
     * @brief Access element with [] notation. 
     * @param i Index of the element (0-indexed).
     * @return Reference to the i-th element.
     */
    inline T& operator[](int i) {
        CFEM_ASSERT(i >= 0 && i < N);
        return data[i];
    }

    /*! 
     * @brief Read-only access element with [] notation. 
     * @param i Index of the element (0-indexed).
     * @return Constant reference to the i-th element.
     */
    inline const T& operator[](int i) const {
        CFEM_ASSERT(i >= 0 && i < N);
        return data[i];
    }

    /*! 
     * @brief Returns the fixed size of the vector. 
     * @return Dimension N.
     */
    static constexpr int size() noexcept { return N; }

    /*! 
     * @brief Resets all components to zero.
     */
    void setZero() {
        for (int i = 0; i < N; ++i) data[i] = static_cast<T>(0);
    }

    // --- Arithmetic Operators ---

    /*! 
     * @brief Component-wise addition. 
     */
    FixedVector operator+(const FixedVector& other) const {
        FixedVector result;
        for (int i = 0; i < N; ++i) result.data[i] = data[i] + other.data[i];
        return result;
    }

    /*!
     * @brief Component-wise subtraction. 
     */
    FixedVector operator-(const FixedVector& other) const {
        FixedVector result;
        for (int i = 0; i < N; ++i) result.data[i] = data[i] - other.data[i];
        return result;
    }

    /*!
     * @brief Scalar multiplication. 
     */
    FixedVector operator*(T scalar) const {
        FixedVector result;
        for (int i = 0; i < N; ++i) result.data[i] = data[i] * scalar;
        return result;
    }

    /*! 
     * @brief Scalar division. 
     * @throw CFEM_ASSERT if scalar is near zero.
     */
    FixedVector operator/(T scalar) const {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        FixedVector result;
        for (int i = 0; i < N; ++i) result.data[i] = data[i] / scalar;
        return result;
    }

    /*! 
     * @brief In-place addition. 
     */
    FixedVector& operator+=(const FixedVector& other) {
        for (int i = 0; i < N; ++i) data[i] += other.data[i];
        return *this;
    }

    /*! 
     * @brief In-place subtraction. 
     */
    FixedVector& operator-=(const FixedVector& other) {
        for (int i = 0; i < N; ++i) data[i] -= other.data[i];
        return *this;
    }

    /*! 
     * @brief In-place scalar multiplication. 
     */
    FixedVector& operator*=(T scalar) {
        for (int i = 0; i < N; ++i) data[i] *= scalar;
        return *this;
    }

    /*! 
     * @brief In-place scalar division. 
     */
    FixedVector& operator/=(T scalar) {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        for (int i = 0; i < N; ++i) data[i] /= scalar;
        return *this;
    }

    // --- Linear Algebra Tools ---

    /*! 
     * @brief Computes the dot (inner) product between two vectors. 
     * @param other The vector to dot with.
     * @return The scalar result \f$ \sum v_i \cdot w_i \f$.
     */
    T dot(const FixedVector<T, N>& other) const {
        T result = 0;
        for (int i = 0; i < N; ++i) result += data[i] * other.data[i];
        return result;
    }

    /*! 
     * @brief Returns squared Euclidean norm. 
     * @note Faster than norm() as it avoids std::sqrt.
     */
    T normSquare() const {
        return dot(*this);
    }

    /*! 
     * @brief Returns Euclidean norm (L2 length). 
     */
    T norm() const {
        return std::sqrt(normSquare());
    }
};

/*! 
 * @brief Left-side scalar multiplication. 
 * @relates FixedVector
 */
template <typename T, int N>
inline FixedVector<T, N> operator*(T scalar, const FixedVector<T, N>& v) {
    return v * scalar;
}

/*! 
 * @brief Formats the vector for standard output. 
 * @relates FixedVector
 */
template <typename T, int N>
inline std::ostream& operator<<(std::ostream& os, const FixedVector<T, N>& v) {
    os << "(";
    for (int i = 0; i < N; ++i) {
        os << v[i] << (i == N - 1 ? "" : ", ");
    }
    os << ")";
    return os;
}


} // namespace
#endif // CFEM_STATIC_VECTOR