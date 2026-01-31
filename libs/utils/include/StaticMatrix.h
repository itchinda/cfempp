/**
 * @file StaticMatrix.h
 * @author Itchinda
 * @brief Fixed-size dense matrix for local element operations.
 * @date 2026-01-27
 */

#ifndef CFEM_STATIC_MATRIX
#define CFEM_STATIC_MATRIX

#include "cfemutils.h"
#include "StaticVector.h"
#include <iostream>
#include <utility>


/**
 * @brief A stack-allocated, fixed-size dense matrix.
 * @tparam T The numeric type (e.g., CfemReal).
 * @tparam R Number of rows.
 * @tparam C Number of columns.
 */
template <typename T, int R, int C>
struct StaticMatrix {
    T data[R * C] = {0};

    StaticMatrix() = default;

    StaticMatrix(std::initializer_list<T> list) {
        CFEM_ASSERT(list.size() == R * C);
        int i = 0;
        for (auto val : list) data[i++] = val;
    }

    // --- Accessors ---

    inline T& operator()(int r, int c) {
        CFEM_ASSERT(r < R && c < C);
        return data[r * C + c];
    }

    inline const T& operator()(int r, int c) const {
        CFEM_ASSERT(r < R && c < C);
        return data[r * C + c];
    }

    static constexpr int rows() noexcept { return R; }
    static constexpr int cols() noexcept { return C; }
    static constexpr std::pair<int,int> dimensions() noexcept { return {R, C}; }

    // --- Basic Operations ---

    void setZero() {
        for (int i = 0; i < R * C; ++i) data[i] = static_cast<T>(0);
    }

    void setIdentity() {
        static_assert(R == C, "Identity matrix must be square.");
        setZero();
        for (int i = 0; i < R; ++i) (*this)(i, i) = static_cast<T>(1);
    }

    // --- Arithmetic ---

    StaticMatrix<T, R, C> operator*(T scalar) const {
        StaticMatrix<T, R, C> res;
        for (int i = 0; i < R * C; ++i) res.data[i] = data[i] * scalar;
        return res;
    }

    StaticMatrix<T, R, C> operator/(T scalar) const {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        StaticMatrix<T, R, C> res;
        for (int i = 0; i < R * C; ++i) res.data[i] = data[i] / scalar;
        return res;
    }

    StaticMatrix<T, R, C> operator+(const StaticMatrix& other) const {
        StaticMatrix<T, R, C> res;
        for (int i = 0; i < R * C; ++i) res.data[i] = data[i] + other.data[i];
        return res;
    }

    StaticMatrix<T, R, C> operator-(const StaticMatrix& other) const {
        StaticMatrix<T, R, C> res;
        for (int i = 0; i < R * C; ++i) res.data[i] = data[i] - other.data[i];
        return res;
    }

    StaticMatrix<T, R, C>& operator*=(T scalar) {

        for (int i = 0; i < R * C; ++i) data[i] *= scalar;
        return *this;
    }

    StaticMatrix<T, R, C>& operator/=(T scalar) {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        for (int i = 0; i < R * C; ++i) data[i] /= scalar;
        return *this;
    }

    StaticMatrix<T, R, C>& operator+=(const StaticMatrix& other) {
        for (int i = 0; i < R * C; ++i) data[i] += other.data[i];
        return *this;
    }

    StaticMatrix<T, R, C>& operator-=(const StaticMatrix& other) {
        for (int i = 0; i < R * C; ++i) data[i] -= other.data[i];
        return *this;
    }


    // --- Transpose ---

    /** @brief Returns a new matrix that is the transpose of this. */
    StaticMatrix<T, C, R> transpose() const {
        StaticMatrix<T, C, R> res;
        for (int i = 0; i < R; ++i)
            for (int j = 0; j < C; ++j)
                res(j, i) = (*this)(i, j);
        return res;
    }
};

// --- Multiplications (Non-member for template flexibility) ---

/** @brief Matrix-Vector Product: res(R) = Mat(R,C) * Vec(C) */
template <typename T, int R, int C>
inline StaticVector<T, R> operator*(const StaticMatrix<T, R, C>& m, const StaticVector<T, C>& v) {
    StaticVector<T, R> res;
    for (int i = 0; i < R; ++i) {
        T sum = 0;
        for (int j = 0; j < C; ++j) {
            sum += m(i, j) * v[j];
        }
        res[i] = sum;
    }
    return res;
}

/** @brief Matrix-Matrix Product: res(R, Inner2) = A(R, Inner1) * B(Inner1, Inner2) */
template <typename T, int R, int Inner, int C2>
inline StaticMatrix<T, R, C2> operator*(const StaticMatrix<T, R, Inner>& A, const StaticMatrix<T, Inner, C2>& B) {
    StaticMatrix<T, R, C2> res;
    res.setZero();
    for (int i = 0; i < R; ++i) {
        for (int k = 0; k < Inner; ++k) {
            T val = A(i, k);
            for (int j = 0; j < C2; ++j) {
                res(i, j) += val * B(k, j);
            }
        }
    }
    return res;
}

/** @brief Scalar multiplication (left side) */
template <typename T, int R, int C>
inline StaticMatrix<T, R, C> operator*(T scalar, const StaticMatrix<T, R, C>& m) {
    return m * scalar;
}

// --- Specialized Tools for Jacobians ---

/** @brief Determinant for 2x2. */
template <typename T>
inline T det(const StaticMatrix<T, 2, 2>& m) {
    return m(0,0)*m(1,1) - m(0,1)*m(1,0);
}

/** @brief Determinant for 3x3. */
template <typename T>
inline T det(const StaticMatrix<T, 3, 3>& m) {
    return m(0,0)*(m(1,1)*m(2,2) - m(1,2)*m(2,1)) -
           m(0,1)*(m(1,0)*m(2,2) - m(1,2)*m(2,0)) +
           m(0,2)*(m(1,0)*m(2,1) - m(1,1)*m(2,0));
}

// --- Output Stream ---

template <typename T, int R, int C>
inline std::ostream& operator<<(std::ostream& os, const StaticMatrix<T, R, C>& m) {
    for (int i = 0; i < R; ++i) {
        os << "| ";
        for (int j = 0; j < C; ++j) {
            os << m(i, j) << (j == C - 1 ? "" : ", ");
        }
        os << " |\n";
    }
    return os;
}

#endif