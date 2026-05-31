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
 * @file FixedMatrix.h
 * @author Itchinda
 * @brief Fixed-size dense matrix for local element operations.
 * @date 2026-01-27
 */

#ifndef CFEM_STATIC_MATRIX
#define CFEM_STATIC_MATRIX

#include <iostream>
#include <utility>

#include "cfemutils.h"
#include "FixedVector.h"

namespace cfem
{

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
    template <typename T, CfemInt R, CfemInt C>
    struct FixedMatrix
    {
        /**
         * @brief Flat array storage to ensure memory contiguity.
         */
        T data[R * C] = {};

        /**
         * @brief Default constructor. Initializes all entries to zero.
         */
        constexpr FixedMatrix() noexcept = default;

        /**
         * @brief Initializer list constructor.
         * @param list Values provided in row-major order.
         * @throw CFEM_ASSERT if list size does not match R*C.
         */
        constexpr FixedMatrix(std::initializer_list<T> list) noexcept(false)
        {
            CFEM_ASSERT(list.size() == R * C, "Invalid initializer_list size.");
            CfemInt i = 0;
            for (const auto &val : list)
                data[i++] = val;
        }

        /**
         * @brief Function call operator for element access (read-write).
         * @param r Row index.
         * @param c Column index.
         * @return Reference to the element at (r, c).
         */
        [[nodiscard]] constexpr T &operator()(CfemInt r, CfemInt c)
        {
            CFEM_ASSERT(r >= 0 && r < R && c >= 0 && c < C, "Matrix Index out of bounds.");
            return data[r * C + c];
        }

        /**
         * @brief Function call operator for element access (read-only).
         * @param r Row index.
         * @param c Column index.
         * @return Constant reference to the element at (r, c).
         */
        [[nodiscard]] constexpr inline const T &operator()(CfemInt r, CfemInt c) const
        {
            CFEM_ASSERT(r >= 0 && r < R && c >= 0 && c < C, "Matrix Index out of bounds.");
            return data[r * C + c];
        }

        [[nodiscard]] constexpr T &at(CfemInt r, CfemInt c)
        {
            return operator()(r, c);
        }

        [[nodiscard]] constexpr const T &at(CfemInt r, CfemInt c) const
        {
            return operator()(r, c);
        }

        /**
         * @brief Returns number of rows.
         */
        [[nodiscard]] static constexpr CfemInt rows() noexcept { return R; }

        /**
         * @brief Returns number of columns.
         */
        [[nodiscard]] static constexpr CfemInt cols() noexcept { return C; }

        /**
         * @brief Returns dimensions as a std::pair {Rows, Cols}.
         */
        [[nodiscard]] static constexpr std::pair<CfemInt, CfemInt> dimensions() noexcept { return {R, C}; }

        // --- Basic Operations ---

        /** @brief Resets all matrix entries to zero. */
        constexpr void setZero() noexcept
        {
            std::fill_n(data, R * C, T{});
        }

        /**
         * @brief Sets the matrix to Identity (\f$ I \f$).
         * @note Only available for square matrices.
         */
        constexpr void setIdentity() noexcept
        {
            CFEM_STATIC_ASSERT(R == C, "Identity matrix must be square.");
            setZero();
            for (CfemInt i = 0; i < R; ++i)
                data[i * C + i] = static_cast<T>(1);
        }

        // --- Arithmetic ---

        /**
         * @brief Matrix-scalar multiplication.
         */
        [[nodiscard]] constexpr FixedMatrix<T, R, C> operator*(T scalar) const noexcept
        {
            FixedMatrix<T, R, C> res;
            for (CfemInt i = 0; i < R * C; ++i)
                res.data[i] = data[i] * scalar;
            return res;
        }

        /**
         * @brief Matrix-scalar division.
         */
        [[nodiscard]] constexpr FixedMatrix<T, R, C> operator/(T scalar) const noexcept
        {
            if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
                CFEM_ERROR << "Division by zero in FixedMatrix operator/";
            FixedMatrix<T, R, C> res;
            for (CfemInt i = 0; i < R * C; ++i)
                res.data[i] = data[i] / scalar;
            return res;
        }

        /**
         * @brief Matrix-matrix addition.
         */
        [[nodiscard]] constexpr FixedMatrix<T, R, C> operator+(const FixedMatrix &other) const noexcept
        {
            FixedMatrix<T, R, C> res;
            for (CfemInt i = 0; i < R * C; ++i)
                res.data[i] = data[i] + other.data[i];
            return res;
        }

        /**
         * @brief Matrix-matrix subtraction.
         */
        [[nodiscard]] constexpr FixedMatrix<T, R, C> operator-(const FixedMatrix &other) const noexcept
        {
            FixedMatrix<T, R, C> res;
            for (CfemInt i = 0; i < R * C; ++i)
                res.data[i] = data[i] - other.data[i];
            return res;
        }

        /**
         * @brief In-place scalar multiplication.
         */
        constexpr FixedMatrix<T, R, C> &operator*=(T scalar) noexcept
        {
            for (CfemInt i = 0; i < R * C; ++i)
                data[i] *= scalar;
            return *this;
        }

        /**
         * @brief In-place scalar division.
         */
        constexpr FixedMatrix<T, R, C> &operator/=(T scalar) noexcept
        {
            if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
                CFEM_ERROR << "Division by zero in FixedMatrix operator/";
            for (CfemInt i = 0; i < R * C; ++i)
                data[i] /= scalar;
            return *this;
        }

        /**
         * @brief In-place matrix addition.
         */
        constexpr FixedMatrix<T, R, C> &operator+=(const FixedMatrix &other) noexcept
        {
            for (CfemInt i = 0; i < R * C; ++i)
                data[i] += other.data[i];
            return *this;
        }

        /**
         * @brief In-place matrix subtraction.
         */
        constexpr FixedMatrix<T, R, C> &operator-=(const FixedMatrix &other) noexcept
        {
            for (CfemInt i = 0; i < R * C; ++i)
                data[i] -= other.data[i];
            return *this;
        }

        template <CfemInt r, CfemInt c>
        [[nodiscard]] constexpr FixedMatrix<T, r, c> block(CfemInt start_row, CfemInt start_col) const 
        {
            CFEM_ASSERT(start_row >= 0 && start_row + r <= R && 
                        start_col >= 0 && start_col + c <= C && 
                        "Block extraction out of bounds!");

            FixedMatrix<T, r, c> res{};

            for (CfemInt i = 0; i < r; ++i) {
                // Précalcul des offsets de ligne (mise en registre)
                const CfemInt parent_row_offset = (start_row + i) * C;
                const CfemInt child_row_offset = i * c;

                for (CfemInt j = 0; j < c; ++j) {
                    res.data[child_row_offset + j] = this->data[parent_row_offset + start_col + j];
                }
            }

            return res;
        }

        /**
         * @brief Computes the transpose of the matrix.
         * @return A new FixedMatrix of dimensions C x R.
         */
        [[nodiscard]] constexpr FixedMatrix<T, C, R> transpose() const noexcept
        {
            FixedMatrix<T, C, R> res;
            for (CfemInt i = 0; i < R; ++i)
            {
                CfemInt offset{i * C};
                for (CfemInt j = 0; j < C; ++j)
                {
                    res.data[j * R + i] = data[offset + j];
                }
            }

            return res;
        }

        /**
         * @brief Computes the trace of the square matrix.
         * @return The trace of the FixedMatrix.
         */
        [[nodiscard]] constexpr T trace() const noexcept
        {
            CFEM_STATIC_ASSERT(R == C, "FixedMatrix::trace() is only defined for square matrices.");
            
            T res = T{};
            for (CfemInt i = 0; i < R; ++i)
            {
                res += data[i * C + i];
            }
            return res;
        }

        /**
         * @brief Compute the double contraction product A:B
         * * @param B The other matrix.
         * @return The scalar result of the double dot product.
         */
        [[nodiscard]] constexpr T doubleDot(const FixedMatrix<T, R, C>& B) const noexcept
        {
            T res = T{};
            
            for (CfemInt i = 0; i < R * C; ++i)
            {
                res += data[i] * B.data[i];
            }
            return res;
        }

        /**
         * @brief Frobenius norm of the matrix.
         */
        [[nodiscard]] constexpr T frobeniusNorm() const noexcept
        {
            T sum = T{};
            for (CfemInt i = 0; i < R * C; ++i)
                sum += data[i] * data[i];
            return std::sqrt(sum);
        }
    };

    /**
     * @brief Matrix-Vector Product: \f$ y = A \cdot x \f$
     * @relates FixedMatrix
     */
    template <typename T, CfemInt R, CfemInt C>
    [[nodiscard]] constexpr FixedVector<T, R> operator*(const FixedMatrix<T, R, C> &m,
                                                        const FixedVector<T, C> &v) noexcept
    {
        FixedVector<T, R> res{};

        for (CfemInt i = 0; i < R; ++i)
        {
            const T *row = &m.data[i * C]; // pointer to the first element in row i
            T sum = T{};
            for (CfemInt j = 0; j < C; ++j)
            {
                sum += row[j] * v.data[j];
            }
            res.data[i] = sum;
        }

        return res;
    }

    /**
     * @brief Matrix-Matrix Product: \f$ C = A \cdot B \f$
     *        Optimized with ikj loop order for better cache performance.
     * @relates FixedMatrix
     */
    template <typename T, CfemInt R, CfemInt Inner, CfemInt C2>
    [[nodiscard]] constexpr FixedMatrix<T, R, C2> operator*(const FixedMatrix<T, R, Inner> &A,
                                                            const FixedMatrix<T, Inner, C2> &B) noexcept
    {
        FixedMatrix<T, R, C2> res{};

        for (CfemInt i = 0; i < R; ++i) {
            T *res_row_i = &res.data[i * C2];

            for (CfemInt k = 0; k < Inner; ++k){
                const T a_ik = A.data[i * Inner + k];
                const T *b_row_k = &B.data[k * C2];

                for (CfemInt j = 0; j < C2; ++j){
                    res_row_i[j] += a_ik * b_row_k[j];
                }
            }
        }

        return res;
    }

    /**
     * @brief Specialized 3x3 Matrix-Vector3D Product: \f$ C = A \cdot B \f$
     * @relates Vector3D
     */
    [[nodiscard]] constexpr inline Vector3D operator*(const FixedMatrix<CfemReal, 3, 3> &A, 
                                                      const Vector3D &B) noexcept
    {
        return {
            A.data[0] * B.x + A.data[1] * B.y + A.data[2] * B.z,
            A.data[3] * B.x + A.data[4] * B.y + A.data[5] * B.z,
            A.data[6] * B.x + A.data[7] * B.y + A.data[8] * B.z};
    }

    /** @brief Left-side scalar multiplication. @relates FixedMatrix */
    template <typename T, CfemInt R, CfemInt C>
    [[nodiscard]] constexpr FixedMatrix<T, R, C> operator*(T scalar, const FixedMatrix<T, R, C> &m) noexcept
    {
        return m * scalar;
    }

    /**
     * @brief Determinant for 2x2 matrix.
     * \f$ det(M) = m_{00}m_{11} - m_{01}m_{10} \f$
     */
    template <typename T>
    [[nodiscard]] constexpr T det(const FixedMatrix<T, 2, 2> &m) noexcept
    {
        return m.data[0] * m.data[3] - m.data[1] * m.data[2];
    }

    /**
     * @brief Determinant for 3x3 matrix using Sarrus-like expansion.
     */
    template <typename T>
    [[nodiscard]] constexpr T det(const FixedMatrix<T, 3, 3> &m) noexcept
    {
        const T m00 = m.data[0], m01 = m.data[1], m02 = m.data[2];
        const T m10 = m.data[3], m11 = m.data[4], m12 = m.data[5];
        const T m20 = m.data[6], m21 = m.data[7], m22 = m.data[8];

        return m00 * (m11 * m22 - m12 * m21) - 
               m01 * (m10 * m22 - m12 * m20) + 
               m02 * (m10 * m21 - m11 * m20);
    }

    /**
     * @brief Computes the inverse of a 3x3 matrix using the adjugate method.
     *
     * @details
     * This overload assumes the determinant is already known and avoids
     * recomputing it. This is useful in performance-critical code paths
     * where the determinant has already been evaluated (e.g. Jacobian-based FEM kernels).
     *
     * No numerical regularization is applied.
     *
     * @warning The caller is responsible for ensuring the matrix is not singular
     * or ill-conditioned.
     *
     * @param m Input matrix.
     * @param determinant Precomputed determinant of @p m.
     * @return Inverse matrix.
     */
    template <typename T>
    [[nodiscard]] constexpr inline FixedMatrix<T, 3, 3> inverse(const FixedMatrix<T, 3, 3>& m, T determinant)
    {
        CFEM_ASSERT(std::abs(determinant) > cfem::ZERO_TOLERANCE, "Singular 3x3 matrix/Jacobian detected.");

        FixedMatrix<T, 3, 3> res;
        T invDet = T(1.0) / determinant;

        const T m00 = m.data[0], m01 = m.data[1], m02 = m.data[2];
        const T m10 = m.data[3], m11 = m.data[4], m12 = m.data[5];
        const T m20 = m.data[6], m21 = m.data[7], m22 = m.data[8];

        // Calcul des cofacteurs et transposition directe (Matrice Adjointe)
        res.data[0] = (m11 * m22 - m12 * m21) * invDet;
        res.data[1] = (m02 * m21 - m01 * m22) * invDet;
        res.data[2] = (m01 * m12 - m02 * m11) * invDet;

        res.data[3] = (m12 * m20 - m10 * m22) * invDet;
        res.data[4] = (m00 * m22 - m02 * m20) * invDet;
        res.data[5] = (m10 * m02 - m00 * m12) * invDet;

        res.data[6] = (m10 * m21 - m11 * m20) * invDet;
        res.data[7] = (m20 * m01 - m00 * m21) * invDet;
        res.data[8] = (m00 * m11 - m10 * m01) * invDet;

        return res;
    }

    /**
     * @brief Computes the inverse of a 3x3 matrix.
     *
     * @details
     * Convenience overload that computes the determinant internally
     * before calling the optimized inverse(m, determinant) implementation.
     *
     * @param m Input matrix.
     * @return Inverse matrix.
     */
    template <typename T>
    [[nodiscard]] constexpr inline FixedMatrix<T, 3, 3> inverse(const FixedMatrix<T, 3, 3> &m)
    {
        return inverse(m,  det(m));
    }

    template <typename T, CfemInt N>
    [[nodiscard]] static constexpr FixedMatrix<T, N, N> identity() noexcept
    {
        FixedMatrix<T, N, N> m;
        m.setIdentity();
        return m;
    }

    template <typename T, CfemInt N>
    [[nodiscard]] static constexpr FixedMatrix<T, N, N> zero() noexcept
    {
        FixedMatrix<T, N, N> m{};
        return m;
    }

    template <typename T, CfemInt N>
    [[nodiscard]] static constexpr FixedMatrix<T, N, N> outerProduct(const FixedVector<T, N> &v1,
                                                                     const FixedVector<T, N> &v2) noexcept
    {
        FixedMatrix<T, N, N> res{};

        for (CfemInt i = 0; i < N; ++i) {
            CfemInt row_i_offset = i * N;
            const T v1_i = v1.data[i];
            for (CfemInt j = 0; j < N; ++j) {
                res.data[row_i_offset + j] = v1_i * v2.data[j];
            }
        }
        return res;
    }

    /**
     * @brief Prints the matrix in a readable grid format.
     */
    template <typename T, CfemInt R, CfemInt C>
    inline std::ostream &operator<<(std::ostream &os, const FixedMatrix<T, R, C> &m)
    {
        for (CfemInt i = 0; i < R; ++i)
        {
            os << "| ";
            for (CfemInt j = 0; j < C; ++j)
            {
                os << m.data[i * C + j] << (j == C - 1 ? "" : ", ");
            }
            os << " |\n";
        }
        return os;
    }


    template <typename T, CfemInt N>
    void solve(const FixedMatrix<T, N, N>& A, const FixedVector<T, N>& b, FixedVector<T, N>& x)
    {
        FixedMatrix<T, N, N> LU = A;  // Copy for in-place factorization

        CfemInt P[N]; // Permutation Vector
        for (CfemInt i = 0; i < N; ++i) {
            P[i] = i;
        }

        // PLU factorization (with partial pivoting)
        for (CfemInt i = 0; i < N; ++i) {
            // Search for pivot row
            CfemInt pivot_row = i;
            const CfemInt row_i_offset = i * N;
            T max_val = std::abs(LU.data[row_i_offset + i]);  // (i,i) component
            
            for (CfemInt j = i + 1; j < N; ++j) {
                T current_abs = std::abs(LU.data[j * N + i]);
                if (current_abs > max_val) {
                    max_val = current_abs;
                    pivot_row = j;
                }
            }

            // CFEM_ERROR si max_val < cfem::ZERO_TOLERANCE (Matrice singulière)
            if (max_val < cfem::ZERO_TOLERANCE) {
                CFEM_ERROR << "LU factorization failed. The matrix is singular.";
            }

            // Rows permutation (if required)
            if (pivot_row != i) {
                const CfemInt pivot_offset = pivot_row * N; // Précalcul
                for (CfemInt k = 0; k < N; ++k) {
                    T tmp = LU.data[row_i_offset + k];  
                    LU.data[row_i_offset + k] = LU.data[pivot_offset + k]; 
                    LU.data[pivot_offset + k] = tmp;
                }
                CfemInt tmp_p = P[i];
                P[i] = P[pivot_row];
                P[pivot_row] = tmp_p;
            }

            // Gauss elimination (In-place LU computation)
            const T pivot_inv = T(1.0) / LU.data[row_i_offset + i]; 
            
            if constexpr (N > 4000){
                #pragma omp parallel for schedule(static)
                for (CfemInt j = i + 1; j < N; ++j) {
                    const CfemInt row_j_offset = j * N;
                    LU.data[row_j_offset + i] *= pivot_inv;
                    
                    const T lji = LU.data[row_j_offset + i];
                    #pragma omp simd
                    for (CfemInt k = i + 1; k < N; ++k) {
                        LU.data[row_j_offset + k] -= lji * LU.data[row_i_offset + k]; // Update U
                    }
                }
            } else {
                for (CfemInt j = i + 1; j < N; ++j) {
                    const CfemInt row_j_offset = j * N;
                    LU.data[row_j_offset + i] *= pivot_inv; 

                    const T lji = LU.data[row_j_offset + i];
                    #pragma omp simd
                    for (CfemInt k = i + 1; k < N; ++k) {
                        LU.data[row_j_offset + k] -= lji * LU.data[row_i_offset + k]; // Update U
                    }
                }

            }
        }

        // Descending step : L * y = P * b
        FixedVector<T, N> y{};
        for (CfemInt i = 0; i < N; ++i) {
            const CfemInt row_i_offset = i * N;
            y[i] = b[P[i]]; // Applying Permutation
            for (CfemInt j = 0; j < i; ++j) {
                y[i] -= LU.data[row_i_offset + j] * y[j];
            }
        }

        // Ascending step : U * x = y
        for (CfemInt i = N; i-- > 0;) {
            const CfemInt row_i_offset = i * N;
            x[i] = y[i];
            for (CfemInt j = i + 1; j < N; ++j) {
                x[i] -= LU.data[row_i_offset + j]  * x[j];
            }
            x[i] *= (T(1.0) / LU.data[row_i_offset + i]); 
        }
    }

} // namespace

#endif // CFEM_STATIC_MATRIX