/**
 * @file FixedMatrix.h
 * @author Itchinda
 * @brief Fixed-size dense matrix for local element operations.
 * @date 2026-01-27
 */

#ifndef CFEM_STATIC_MATRIX
#define CFEM_STATIC_MATRIX

#include "cfemutils.h"
#include "FixedVector.h"
#include <iostream>
#include <utility>


namespace cfem {

/**
 * @struct FixedMatrix
 * @brief A stack-allocated, fixed-size dense matrix.
 *        Storage is row-major by default due to the indexing \f$ i \cdot C + j \f$.
 *        Optimized for small local matrices like Jacobians, local stiffness matrices (\f$ K_e \f$),
 *        and strain-displacement matrices (\f$ B \f$).
 * @tparam T The numeric type (e.g., CfemReal).
 * @tparam R Number of rows.
 * @tparam C Number of columns.
 */
template <typename T, int R, int C>
struct FixedMatrix {
    /** 
     * @brief Flat array storage to ensure memory contiguity. 
     */
    T data[R * C] = {0};

    /** 
     * @brief Default constructor. Initializes all entries to zero. 
     */
    FixedMatrix() = default;

    /** 
     * @brief Initializer list constructor.
     * @param list Values provided in row-major order.
     * @throw CFEM_ASSERT if list size does not match R*C.
     */
    FixedMatrix(std::initializer_list<T> list) {
        CFEM_ASSERT(list.size() == R * C);
        int i = 0;
        for (auto val : list) data[i++] = val;
    }

    /** 
     * @brief Function call operator for element access (read-write).
     * @param r Row index.
     * @param c Column index.
     * @return Reference to the element at (r, c).
     */
    inline T& operator()(int r, int c) {
        CFEM_ASSERT(r < R && c < C);
        return data[r * C + c];
    }

    /** 
     * @brief Function call operator for element access (read-only).
     * @param r Row index.
     * @param c Column index.
     * @return Constant reference to the element at (r, c).
     */
    inline const T& operator()(int r, int c) const {
        CFEM_ASSERT(r < R && c < C);
        return data[r * C + c];
    }

    /** 
     * @brief Returns number of rows. 
     */
    static constexpr int rows() noexcept { return R; }

    /**
     * @brief Returns number of columns. 
     */
    static constexpr int cols() noexcept { return C; }

    /** 
     * @brief Returns dimensions as a std::pair {Rows, Cols}. 
     */
    static constexpr std::pair<int,int> dimensions() noexcept { return {R, C}; }

    // --- Basic Operations ---

    /** @brief Resets all matrix entries to zero. */
    void setZero() {
        for (int i = 0; i < R * C; ++i) data[i] = static_cast<T>(0);
    }

    /** 
     * @brief Sets the matrix to Identity (\f$ I \f$).
     * @note Only available for square matrices.
     */
    void setIdentity() {
        static_assert(R == C, "Identity matrix must be square.");
        setZero();
        for (int i = 0; i < R; ++i) (*this)(i, i) = static_cast<T>(1);
    }

    // --- Arithmetic ---

    /** 
     * @brief Matrix-scalar multiplication. 
     */
    FixedMatrix<T, R, C> operator*(T scalar) const {
        FixedMatrix<T, R, C> res;
        for (int i = 0; i < R * C; ++i) res.data[i] = data[i] * scalar;
        return res;
    }

    /** 
     * @brief Matrix-scalar division. 
     */
    FixedMatrix<T, R, C> operator/(T scalar) const {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        FixedMatrix<T, R, C> res;
        for (int i = 0; i < R * C; ++i) res.data[i] = data[i] / scalar;
        return res;
    }

    /**
     * @brief Matrix-matrix addition.
     */
    FixedMatrix<T, R, C> operator+(const FixedMatrix& other) const {
        FixedMatrix<T, R, C> res;
        for (int i = 0; i < R * C; ++i) res.data[i] = data[i] + other.data[i];
        return res;
    }

    /** 
     * @brief Matrix-matrix subtraction.
     */
    FixedMatrix<T, R, C> operator-(const FixedMatrix& other) const {
        FixedMatrix<T, R, C> res;
        for (int i = 0; i < R * C; ++i) res.data[i] = data[i] - other.data[i];
        return res;
    }

    /** 
     * @brief In-place scalar multiplication. 
     */
    FixedMatrix<T, R, C>& operator*=(T scalar) {
        for (int i = 0; i < R * C; ++i) data[i] *= scalar;
        return *this;
    }

    /** 
     * @brief In-place scalar division.
     */
    FixedMatrix<T, R, C>& operator/=(T scalar) {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        for (int i = 0; i < R * C; ++i) data[i] /= scalar;
        return *this;
    }

    /** 
     * @brief In-place matrix addition. 
     */
    FixedMatrix<T, R, C>& operator+=(const FixedMatrix& other) {
        for (int i = 0; i < R * C; ++i) data[i] += other.data[i];
        return *this;
    }

    /** 
     * @brief In-place matrix subtraction. 
     */
    FixedMatrix<T, R, C>& operator-=(const FixedMatrix& other) {
        for (int i = 0; i < R * C; ++i) data[i] -= other.data[i];
        return *this;
    }

    /** 
     * @brief Computes the transpose of the matrix.
     * @return A new FixedMatrix of dimensions C x R.
     */
    FixedMatrix<T, C, R> transpose() const {
        FixedMatrix<T, C, R> res;
        for (int i = 0; i < R; ++i)
            for (int j = 0; j < C; ++j)
                res(j, i) = (*this)(i, j);
        return res;
    }
};

/** 
 * @brief Matrix-Vector Product: \f$ y = A \cdot x \f$
 * @relates FixedMatrix
 */
template <typename T, int R, int C>
inline FixedVector<T, R> operator*(const FixedMatrix<T, R, C>& m, const FixedVector<T, C>& v) {
    FixedVector<T, R> res;
    for (int i = 0; i < R; ++i) {
        T sum = 0;
        for (int j = 0; j < C; ++j) {
            sum += m(i, j) * v[j];
        }
        res[i] = sum;
    }
    return res;
}

/** 
 * @brief Matrix-Matrix Product: \f$ C = A \cdot B \f$
 *        Optimized with ikj loop order for better cache performance.
 * @relates FixedMatrix
 */
template <typename T, int R, int Inner, int C2>
inline FixedMatrix<T, R, C2> operator*(const FixedMatrix<T, R, Inner>& A, const FixedMatrix<T, Inner, C2>& B) {
    FixedMatrix<T, R, C2> res;
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

/** 
 * @brief 3x3 Matrix-Vector3D Product: \f$ C = A \cdot B \f$
 *        Optimized with ikj loop order for better cache performance.
 * @relates Vector3D
 */
inline Vector3D operator*(const FixedMatrix<CfemReal, 3, 3>& A, const Vector3D& B) {
    return { 
            A(0,0)*B.x + A(0,1)*B.y + A(0,2)*B.z, 
            A(1,0)*B.x + A(1,1)*B.y + A(1,2)*B.z, 
            A(2,0)*B.x + A(2,1)*B.y + A(2,2)*B.z
        };
}

/** @brief Left-side scalar multiplication. @relates FixedMatrix */
template <typename T, int R, int C>
inline FixedMatrix<T, R, C> operator*(T scalar, const FixedMatrix<T, R, C>& m) {
    return m * scalar;
}

// --- Specialized Tools for Jacobians ---

/** 
 * @brief Determinant for 2x2 matrix.
 * \f$ det(M) = m_{00}m_{11} - m_{01}m_{10} \f$
 */
template <typename T>
inline T det(const FixedMatrix<T, 2, 2>& m) {
    return m(0,0)*m(1,1) - m(0,1)*m(1,0);
}

/** 
 * @brief Determinant for 3x3 matrix using Sarrus-like expansion.
 */
template <typename T>
inline T det(const FixedMatrix<T, 3, 3>& m) {
    return m(0,0)*(m(1,1)*m(2,2) - m(1,2)*m(2,1)) -
           m(0,1)*(m(1,0)*m(2,2) - m(1,2)*m(2,0)) +
           m(0,2)*(m(1,0)*m(2,1) - m(1,1)*m(2,0));
}


/** 
 * @brief Calcule l'inverse d'une matrice 3x3 en utilisant la comatrice.
 * @param m La matrice à inverser.
 * @param determinant Le déterminant déjà calculé (pour éviter de le recalculer).
 * @return La matrice inverse.
 * @throw CFEM_ASSERT si le déterminant est proche de zéro.
 */
template <typename T>
inline FixedMatrix<T, 3, 3> inverse(const FixedMatrix<T, 3, 3>& m, T determinant) {
    CFEM_ASSERT(std::abs(determinant) > cfem::ZERO_TOLERANCE && "Tentative d'inversion d'une matrice singulière (J nul ou élément retourné).");

    FixedMatrix<T, 3, 3> res;
    T invDet = static_cast<T>(1.0) / determinant;

    // Calcul des cofacteurs et transposition directe (Matrice Adjointe)
    res(0, 0) = (m(1, 1) * m(2, 2) - m(1, 2) * m(2, 1)) * invDet;
    res(0, 1) = (m(0, 2) * m(2, 1) - m(0, 1) * m(2, 2)) * invDet;
    res(0, 2) = (m(0, 1) * m(1, 2) - m(0, 2) * m(1, 1)) * invDet;

    res(1, 0) = (m(1, 2) * m(2, 0) - m(1, 0) * m(2, 2)) * invDet;
    res(1, 1) = (m(0, 0) * m(2, 2) - m(0, 2) * m(2, 0)) * invDet;
    res(1, 2) = (m(1, 0) * m(0, 2) - m(0, 0) * m(1, 2)) * invDet;

    res(2, 0) = (m(1, 0) * m(2, 1) - m(1, 1) * m(2, 0)) * invDet;
    res(2, 1) = (m(2, 0) * m(0, 1) - m(0, 0) * m(2, 1)) * invDet;
    res(2, 2) = (m(0, 0) * m(1, 1) - m(1, 0) * m(0, 1)) * invDet;

    return res;
}

/** 
 * @brief Prints the matrix in a readable grid format. 
 */
template <typename T, int R, int C>
inline std::ostream& operator<<(std::ostream& os, const FixedMatrix<T, R, C>& m) {
    for (int i = 0; i < R; ++i) {
        os << "| ";
        for (int j = 0; j < C; ++j) {
            os << m(i, j) << (j == C - 1 ? "" : ", ");
        }
        os << " |\n";
    }
    return os;
}

} // namespace

#endif // CFEM_STATIC_MATRIX