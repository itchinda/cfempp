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
 * @file DynamicMatrix.h
 * @author Itchinda
 * @brief Dynamic-size dense matrix for general Finite Element operations.
 * @version 0.1
 * @date 2026-01-27
 */

#ifndef CFEM_DYNAMIC_MATRIX
#define CFEM_DYNAMIC_MATRIX

#include "cfemutils.h"
#include "FixedMatrix.h"
#include "DynamicVector.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <span>

namespace cfem
{

    template <typename T>
    class MatrixView
    {
    private:
        T *m_data = nullptr;

        CfemInt m_rows = 0;
        CfemInt m_cols = 0;

        // Distance entre deux lignes successives
        // (en nombre d’éléments)
        CfemInt m_stride = 0;

    public:
        using value_type = T;

        constexpr MatrixView() noexcept = default;

        constexpr MatrixView( T *data,
                              CfemInt rows,
                              CfemInt cols,
                              CfemInt stride) noexcept
            : m_data(data),
              m_rows(rows),
              m_cols(cols),
              m_stride(stride)
        { }

        // Cas matrice compacte contiguë
        constexpr MatrixView( T *data,
                              CfemInt rows,
                              CfemInt cols) noexcept
            : m_data(data),
              m_rows(rows),
              m_cols(cols),
              m_stride(cols)
        {}

        constexpr CfemInt n_rows() const noexcept { return m_rows; }
        constexpr CfemInt n_cols() const noexcept { return m_cols; }
        constexpr CfemInt stride() const noexcept { return m_stride; }
        constexpr size_t  size()   const noexcept { return static_cast<size_t>(m_rows) * static_cast<size_t>(m_cols); }
        constexpr bool    empty()  const noexcept { return size() == 0; }

        constexpr T *data() noexcept{ return m_data; }
        constexpr const T *data() const noexcept{ return m_data; }

        constexpr T *rowData(CfemInt r)
        {    
            CFEM_ASSERT(r >= 0 && r < m_rows);   
             return m_data + r * m_stride;
        }

        constexpr const T *rowData(CfemInt r) const
        {
            CFEM_ASSERT(r >= 0 && r < m_rows);
            return m_data + r * m_stride;
        }

        std::span<T> operator[](CfemInt r)
        {
            CFEM_ASSERT(r >= 0 && r < m_rows);
            return std::span<T>( rowData(r), static_cast<size_t>(m_cols) );
        }

        std::span<const T> operator[](CfemInt r) const
        {
            CFEM_ASSERT(r >= 0 && r < m_rows);
            return std::span<const T>( rowData(r), static_cast<size_t>(m_cols) );
        }

        constexpr T &operator()(CfemInt r, CfemInt c)
        {
            CFEM_ASSERT(r >= 0 && r < m_rows);
            CFEM_ASSERT(c >= 0 && c < m_cols);

            return m_data[r * m_stride + c];
        }

        constexpr const T &operator()(CfemInt r, CfemInt c) const
        {
            CFEM_ASSERT(r >= 0 && r < m_rows);
            CFEM_ASSERT(c >= 0 && c < m_cols);

            return m_data[r * m_stride + c];
        }

        std::span<T> row(CfemInt r) noexcept {
            return (*this)[r];
        }

        std::span<const T> row(CfemInt r) const noexcept {
            return (*this)[r];
        }

        void fill(const T &value)
        {
            for (CfemInt i = 0; i < m_rows; ++i){
                auto row = (*this)[i];
                std::fill( row.begin(), row.end(), value);
            }
        }

        constexpr MatrixView subView( CfemInt row0,
                                      CfemInt col0,
                                      CfemInt subRows,
                                      CfemInt subCols) const
        {
            CFEM_ASSERT(row0 >= 0);
            CFEM_ASSERT(col0 >= 0);
            CFEM_ASSERT(subRows >= 0);
            CFEM_ASSERT(subCols >= 0);
            CFEM_ASSERT(row0 + subRows <= m_rows);
            CFEM_ASSERT(col0 + subCols <= m_cols);
            CFEM_ASSERT(m_data != nullptr || subRows * subCols == 0);

            return MatrixView( m_data + row0 * m_stride + col0,
                               subRows,
                               subCols,
                               m_stride);
        }

        constexpr bool isContiguous() const noexcept {
            return m_stride == m_cols;
        }

        [[nodiscard]] constexpr bool sameShape( const MatrixView& other) const noexcept {
            return m_rows == other.m_rows && m_cols == other.m_cols;
        }
    };

    template <typename T>
    using ConstMatrixView = MatrixView<const T>;

    /**
     * @class DynamicMatrix
     * @brief A heap-allocated, dynamically-sized dense matrix.
     * * Storage is row-major: index = r * cols + c.
     * @tparam T The numeric type (e.g., CfemReal, double).
     */
    template <typename T>
    class DynamicMatrix
    {
    private:
        std::vector<T> m_data;
        CfemInt m_rows = 0;
        CfemInt m_cols = 0;

    public:
        // --- Constructors ---

        /**
         * @brief Default constructor (empty matrix).
         */
        DynamicMatrix() = default;

        /**
         * @brief Size constructor.
         * @param r Number of rows.
         * @param c Number of columns.
         * @param init Initial value for all elements.
         */
        DynamicMatrix(CfemInt r, CfemInt c, T init = static_cast<T>(0))
            : m_rows(r), m_cols(c), m_data(static_cast<size_t>(r * c), init) {}

        /**
         * @brief Initializer list constructor (row-major).
         */
        DynamicMatrix(CfemInt r, CfemInt c, std::initializer_list<T> list)
            : m_rows(r), m_cols(c), m_data(list)
        {
            CFEM_ASSERT(m_data.size() == static_cast<size_t>(r * c));
        }

        /**
         * @brief Conversion constructor from FixedMatrix.
         */
        template <CfemInt R, CfemInt C>
        DynamicMatrix(const FixedMatrix<T, R, C> &fixed)
            : m_rows(R), m_cols(C), m_data(std::begin(fixed.data), std::end(fixed.data)) {}

        // Move constructor
        DynamicMatrix(DynamicMatrix &&other) noexcept
            : m_data(std::move(other.m_data)), m_rows(other.m_rows), m_cols(other.m_cols)
        {
            other.m_rows = 0;
            other.m_cols = 0;
        }

        // move assignment
        DynamicMatrix &operator=(DynamicMatrix &&other) noexcept
        {
            if (this != &other)
            {
                m_data = std::move(other.m_data);
                m_rows = other.m_rows;
                m_cols = other.m_cols;
                other.m_rows = 0;
                other.m_cols = 0;
            }
            return *this;
        }

        DynamicMatrix(const DynamicMatrix &other) = default;
        DynamicMatrix &operator=(const DynamicMatrix &other) = default;

        [[nodiscard]] constexpr inline CfemInt n_rows() const noexcept { return m_rows; }
        [[nodiscard]] constexpr inline CfemInt n_cols() const noexcept { return m_cols; }
        [[nodiscard]] constexpr inline size_t  capacity() const noexcept { return m_data.capacity(); }
        [[nodiscard]] constexpr inline size_t  size() const noexcept { return m_rows * m_cols; }

        void clear()
        {
            m_data.clear();
            m_rows = 0;
            m_cols = 0;
        }

        /**
         * @brief Resizes the matrix. 
         * @warning This operation preserves the underlying 1D memory but DESTROYS the 2D layout 
         * if the number of columns changes. Old data should be considered invalidated.
         */
        inline void resize(CfemInt r, CfemInt c)
        {
            CFEM_ASSERT(r >= 0, "rows size should be a positive integer");
            CFEM_ASSERT(c >= 0, "columns size should be a positive integer");
            m_rows = r;
            m_cols = c;
            m_data.resize(static_cast<size_t>(r * c));
        }

        inline void reserve(CfemInt r, CfemInt c)
        {
            CFEM_ASSERT(r >= 0, "rows size should be a positive integer");
            CFEM_ASSERT(c >= 0, "columns size should be a positive integer");

            m_data.reserve(static_cast<size_t>(r * c));
        }

        inline void reshape(CfemInt r, CfemInt c)
        {
            CFEM_ASSERT(r >= 0);
            CFEM_ASSERT(c >= 0);
            CFEM_ASSERT( static_cast<size_t>(r) * static_cast<size_t>(c) == size() );
            m_rows = r;
            m_cols = c;
        }

        /**
         * @brief Redimensionne et initialise la matrice.
         * @param rows Nombre de lignes
         * @param cols Nombre de colonnes
         * @param value Valeur d'initialisation (par défaut 0)
         */
        inline void assign(CfemInt r, CfemInt c, const T &value = T())
        {
            CFEM_ASSERT(r >= 0, "rows size should be a positive integer");
            CFEM_ASSERT(c >= 0, "columns size should be a positive integer");

            m_rows = r;
            m_cols = c;
            m_data.assign(static_cast<size_t>(r * c), value);
        }

        /**
         * @brief Resets all entries to zero.
         */
        void setZero() { std::fill(m_data.begin(), m_data.end(), static_cast<T>(0)); }

        /**
         * @brief Sets a square matrix to identity.
         */
        void setIdentity()
        {
            CFEM_ASSERT(m_rows == m_cols);
            setZero();
            for (CfemInt i = 0; i < m_rows; ++i)
                (*this)(i, i) = static_cast<T>(1);
        }

        [[nodiscard]] T *data() noexcept { return m_data.data(); }
        [[nodiscard]] const T *data() const noexcept { return m_data.data(); }

        /** @brief Itérateur sur le début des données. */
        auto begin() noexcept { return m_data.begin(); }
        auto end() noexcept { return m_data.end(); }

        /** @brief Itérateur constant sur le début des données. */
        auto begin() const noexcept { return m_data.begin(); }
        auto end() const noexcept { return m_data.end(); }

        /** @brief Itérateur constant explicite (cbegin/cend). */
        auto cbegin() const noexcept { return m_data.cbegin(); }
        auto cend() const noexcept { return m_data.cend(); }

        inline std::span<T> operator[](CfemInt r)
        {
            CFEM_ASSERT(r >= 0 && r < m_rows);
            return std::span<T>(m_data.data() + r * m_cols, static_cast<size_t>(m_cols));
        }

        [[nodiscard]] inline std::span<const T> operator[](CfemInt r) const
        {
            CFEM_ASSERT(r >= 0 && r < m_rows);
            return std::span<const T>(m_data.data() + r * m_cols, static_cast<size_t>(m_cols));
        }

        std::span<T> row(CfemInt r) {
            CFEM_ASSERT(r >= 0 && r < m_rows);
            return (*this)[r];
        }

        [[nodiscard]] std::span<const T> row(CfemInt r) const {
            CFEM_ASSERT(r >= 0 && r < m_rows);
            return (*this)[r];
        }

        /**
         * @brief Element access (read-write).
         */
        inline T &operator()(CfemInt r, CfemInt c)
        {
            CFEM_ASSERT(r >= 0 && r < m_rows && c >= 0 && c < m_cols);
            return m_data[static_cast<size_t>(r * m_cols + c)];
        }

        /**
         * @brief Element access (read-only).
         */
        [[nodiscard]] inline const T &operator()(CfemInt r, CfemInt c) const
        {
            CFEM_ASSERT(r >= 0 && r < m_rows && c >= 0 && c < m_cols);
            return m_data[static_cast<size_t>(r * m_cols + c)];
        }

        [[nodiscard]] MatrixView<T> view() {
            return MatrixView<T>(data(), n_rows(), n_cols(), n_cols() );
        }

        [[nodiscard]] ConstMatrixView<T> view() const {
            return ConstMatrixView<T>(data(), n_rows(), n_cols(), n_cols() );
        }

        [[nodiscard]] MatrixView<T> view(CfemInt row0,
                                         CfemInt col0,
                                         CfemInt subRows,
                                         CfemInt subCols ) {
            return view().subView( row0, col0, subRows, subCols );
        }

        [[nodiscard]] ConstMatrixView<T> view( CfemInt row0,
                                               CfemInt col0,
                                               CfemInt subRows,
                                               CfemInt subCols ) const {
            return view().subView( row0, col0, subRows, subCols );
        }

        [[nodiscard]] bool empty() const noexcept { return m_data.empty(); }
        [[nodiscard]] bool sameShape(const DynamicMatrix& other) const noexcept {
            return m_rows == other.m_rows && m_cols == other.m_cols;
        }
        [[nodiscard]] bool isSquare() const noexcept { return m_rows == m_cols; }

        // --- Arithmetic ---

        [[nodiscard]] DynamicMatrix operator+(const DynamicMatrix &other) const
        {
            CFEM_ASSERT(m_rows == other.m_rows && m_cols == other.m_cols);
            DynamicMatrix res(m_rows, m_cols);
            for (size_t i = 0; i < m_data.size(); ++i)
                res.m_data[i] = m_data[i] + other.m_data[i];
            return res;
        }

        [[nodiscard]] DynamicMatrix operator-(const DynamicMatrix &other) const
        {
            CFEM_ASSERT(m_rows == other.m_rows && m_cols == other.m_cols);
            DynamicMatrix res(m_rows, m_cols);
            for (size_t i = 0; i < m_data.size(); ++i)
                res.m_data[i] = m_data[i] - other.m_data[i];
            return res;
        }

        [[nodiscard]] DynamicMatrix operator*(T scalar) const
        {
            DynamicMatrix res(m_rows, m_cols);
            for (size_t i = 0; i < m_data.size(); ++i)
                res.m_data[i] = m_data[i] * scalar;
            return res;
        }

        [[nodiscard]] DynamicMatrix operator/(T scalar) const
        {
            CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
            DynamicMatrix res(m_rows, m_cols);
            T invScalar = static_cast<T>(1) / scalar;
            for (size_t i = 0; i < m_data.size(); ++i)
                res.m_data[i] = m_data[i] * invScalar;
            return res;
        }

        // Compound Operators

        DynamicMatrix &operator+=(const DynamicMatrix &other)
        {
            CFEM_ASSERT(m_rows == other.m_rows && m_cols == other.m_cols);
            for (size_t i = 0; i < m_data.size(); ++i)
                m_data[i] += other.m_data[i];
            return *this;
        }

        DynamicMatrix &operator-=(const DynamicMatrix &other)
        {
            CFEM_ASSERT(m_rows == other.m_rows && m_cols == other.m_cols);
            for (size_t i = 0; i < m_data.size(); ++i)
                m_data[i] -= other.m_data[i];
            return *this;
        }

        DynamicMatrix &operator*=(T scalar)
        {
            for (auto &val : m_data)
                val *= scalar;
            return *this;
        }

        DynamicMatrix &operator/=(T scalar)
        {
            CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
            T invScalar = static_cast<T>(1) / scalar;
            for (auto &val : m_data)
                val *= invScalar;
            return *this;
        }

        // Linear Algebra

        inline void transpose(DynamicMatrix &out) const
        {
            //Case M.transpose(M);
            if (this == &out) 
            {
                if (isSquare()) {
                    for (CfemInt i = 0; i < m_rows; ++i) {
                        for (CfemInt j = i + 1; j < m_cols; ++j) {
                            std::swap(out(i, j), out(j, i));
                        }
                    }
                } else {
                    // Fallback for rectangular matrices
                    DynamicMatrix temp = this->transpose();
                    out = std::move(temp);
                }
                return;
            }

            // --- Comportement normal out-of-place (this != &out) ---
            if (out.n_rows() != m_cols || out.n_cols() != m_rows)
                out.resize(m_cols, m_rows);

            for (CfemInt i = 0; i < m_rows; ++i) {
                for (CfemInt j = 0; j < m_cols; ++j) {
                    out(j, i) = (*this)(i, j);
                }
            }
        }

        /**
         * @brief Returns a new transposed matrix.
         */
        inline DynamicMatrix transpose() const
        {
            DynamicMatrix res;
            transpose(res);
            return res;
        }
    };

    // Non-member Operations

    /**
     * @brief Matrix-Vector multiplication: y = A * x
     */
    template <typename T>
    [[nodiscard]] inline DynamicVector<T> operator*(const DynamicMatrix<T> &A, const DynamicVector<T> &x)
    {
        CFEM_ASSERT(A.n_cols() == x.size());
        DynamicVector<T> res(A.n_rows());
        for (CfemInt i = 0; i < A.n_rows(); ++i) {
            T sum = 0;
            for (CfemInt j = 0; j < A.n_cols(); ++j) {
                sum += A(i, j) * x[j];
            }
            res[i] = sum;
        }
        return res;
    }

    /**
     * @brief Matrix-Matrix multiplication: C = A * B (ikj optimized)
     *
     * @tparam T The matrix element type
     * @param A The left matrix
     * @param B The right matrix
     * @return DynamicMatrix<T> the result of the multiplication
     */
    template <typename T>
    inline DynamicMatrix<T> operator*(const DynamicMatrix<T> &A, const DynamicMatrix<T> &B)
    {
        CFEM_ASSERT(A.n_cols() == B.n_rows());

        DynamicMatrix<T> res(A.n_rows(), B.n_cols(), static_cast<T>(0));
        for (CfemInt i = 0; i < A.n_rows(); ++i) {
            for (CfemInt k = 0; k < A.n_cols(); ++k) {
                T val = A(i, k);
                for (CfemInt j = 0; j < B.n_cols(); ++j) {
                    res(i, j) += val * B(k, j);
                }
            }
        }
        return res;
    }

    /**
     * @brief Left-side scalar multiplication.
     */
    template <typename T>
    inline DynamicMatrix<T> operator*(T scalar, const DynamicMatrix<T> &m)
    {
        return m * scalar;
    }

    /**
     * @brief Output stream formatting.
     */
    template <typename T>
    inline std::ostream &operator<<(std::ostream &os, const DynamicMatrix<T> &m)
    {
        os << "[size=" << m.n_rows() << "x" << m.n_cols() << "]\n";
        for (CfemInt i = 0; i < m.n_rows(); ++i) {
            os << "| ";
            for (CfemInt j = 0; j < m.n_cols(); ++j) {
                os << m(i, j) << (j == m.n_cols() - 1 ? "" : ", ");
            }
            os << " |\n";
        }
        return os;
    }

} // namespace

#endif // CFEM_DYNAMIC_MATRIX