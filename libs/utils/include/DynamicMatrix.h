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

namespace cfem {

/**
 * @class DynamicMatrix
 * @brief A heap-allocated, dynamically-sized dense matrix.
 * * Storage is row-major: index = r * cols + c.
 * @tparam T The numeric type (e.g., CfemReal, double).
 */
template <typename T>
class DynamicMatrix {
private:
    std::vector<T> m_data;
    int m_rows = 0;
    int m_cols = 0;

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
    DynamicMatrix(int r, int c, T init = static_cast<T>(0)) 
        : m_rows(r), m_cols(c), m_data(static_cast<size_t>(r * c), init) {}

    /** 
     * @brief Initializer list constructor (row-major). 
     */
    DynamicMatrix(int r, int c, std::initializer_list<T> list) 
        : m_rows(r), m_cols(c), m_data(list) {
        CFEM_ASSERT(m_data.size() == static_cast<size_t>(r * c));
    }

    /**
     * @brief Conversion constructor from FixedMatrix.
     */
    template <int R, int C>
    DynamicMatrix(const FixedMatrix<T, R, C>& fixed) 
        : m_rows(R), m_cols(C), m_data(std::begin(fixed.data), std::end(fixed.data)) {}

    // Move constructor
    DynamicMatrix(DynamicMatrix&& other) noexcept 
        : m_data(std::move(other.m_data)), m_rows(other.m_rows), m_cols(other.m_cols) {
        other.m_rows = 0;
        other.m_cols = 0;
    }

    // move assignment
    DynamicMatrix& operator=(DynamicMatrix&& other) noexcept {
        if (this != &other) {
            m_data = std::move(other.m_data);
            m_rows = other.m_rows;
            m_cols = other.m_cols;
            other.m_rows = 0;
            other.m_cols = 0;
        }
        return *this;
    }

    DynamicMatrix(const DynamicMatrix& other) = default;
    DynamicMatrix& operator=(const DynamicMatrix& other) = default;

    inline int rows() const noexcept { return m_rows; }
    inline int cols() const noexcept { return m_cols; }

    void clear() {
        m_data.clear();
        m_rows = 0;
        m_cols = 0;
    }

    /** 
     * @brief Resizes the matrix. Note: data might be lost during resize.
     */
    void resize(int r, int c) {
        m_rows = r;
        m_cols = c;
        m_data.resize(static_cast<size_t>(r * c));
    }

    void reserve(int r, int c) { 
        m_data.reserve(static_cast<size_t>(r * c)); 
    }

    /**
     * @brief Redimensionne et initialise la matrice.
     * @param rows Nombre de lignes
     * @param cols Nombre de colonnes
     * @param value Valeur d'initialisation (par défaut 0)
     */
    void assign(int rows, int cols, const T& value = T()) {
        m_rows = rows;
        m_cols = cols;
        m_data.assign(rows * cols, value);
    }


    /** 
     * @brief Resets all entries to zero.
     */
    void setZero() { std::fill(m_data.begin(), m_data.end(), static_cast<T>(0)); }

    /** 
     * @brief Sets a square matrix to identity.
     */
    void setIdentity() {
        CFEM_ASSERT(m_rows == m_cols);
        setZero();
        for (int i = 0; i < m_rows; ++i) (*this)(i, i) = static_cast<T>(1);
    }

    T* data() noexcept { return m_data.data(); }
    const T* data() const noexcept { return m_data.data(); }


    /** @brief Itérateur sur le début des données. */
    auto begin() noexcept { return m_data.begin(); }
    auto end() noexcept { return m_data.end(); }

    /** @brief Itérateur constant sur le début des données. */
    auto begin() const noexcept { return m_data.begin(); }
    auto end() const noexcept { return m_data.end(); }

    /** @brief Itérateur constant explicite (cbegin/cend). */
    auto cbegin() const noexcept { return m_data.cbegin(); }
    auto cend() const noexcept { return m_data.cend(); }

    /**
     * @brief Element access (read-write). 
     */
    inline T& operator()(int r, int c) {
        CFEM_ASSERT(r >= 0 && r < m_rows && c >= 0 && c < m_cols);
        return m_data[static_cast<size_t>(r * m_cols + c)];
    }

    /** 
     * @brief Element access (read-only). 
     */
    inline const T& operator()(int r, int c) const {
        CFEM_ASSERT(r >= 0 && r < m_rows && c >= 0 && c < m_cols);
        return m_data[static_cast<size_t>(r * m_cols + c)];
    }

    // --- Arithmetic ---

    DynamicMatrix operator+(const DynamicMatrix& other) const {
        CFEM_ASSERT(m_rows == other.m_rows && m_cols == other.m_cols);
        DynamicMatrix res(m_rows, m_cols);
        for (size_t i = 0; i < m_data.size(); ++i) res.m_data[i] = m_data[i] + other.m_data[i];
        return res;
    }

    DynamicMatrix operator-(const DynamicMatrix& other) const {
        CFEM_ASSERT(m_rows == other.m_rows && m_cols == other.m_cols);
        DynamicMatrix res(m_rows, m_cols);
        for (size_t i = 0; i < m_data.size(); ++i) res.m_data[i] = m_data[i] - other.m_data[i];
        return res;
    }

    DynamicMatrix operator*(T scalar) const {
        DynamicMatrix res(m_rows, m_cols);
        for (size_t i = 0; i < m_data.size(); ++i) res.m_data[i] = m_data[i] * scalar;
        return res;
    }

    DynamicMatrix operator/(T scalar) const {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        DynamicMatrix res(m_rows, m_cols);
        T invScalar = static_cast<T>(1) / scalar;
        for (size_t i = 0; i < m_data.size(); ++i) res.m_data[i] = m_data[i] * invScalar;
        return res;
    }

    

    // --- Compound Operators ---

    DynamicMatrix& operator+=(const DynamicMatrix& other) {
        CFEM_ASSERT(m_rows == other.m_rows && m_cols == other.m_cols);
        for (size_t i = 0; i < m_data.size(); ++i) m_data[i] += other.m_data[i];
        return *this;
    }

    DynamicMatrix& operator-=(const DynamicMatrix& other) {
        CFEM_ASSERT(m_rows == other.m_rows && m_cols == other.m_cols);
        for (size_t i = 0; i < m_data.size(); ++i) m_data[i] -= other.m_data[i];
        return *this;
    }

    DynamicMatrix& operator*=(T scalar) {
        for (auto& val : m_data) val *= scalar;
        return *this;
    }

    DynamicMatrix& operator/=(T scalar) {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        T invScalar = static_cast<T>(1) / scalar;
        for (auto& val : m_data) val *= invScalar;
        return *this;
    }

    // --- Linear Algebra ---

    /** 
     * @brief Returns a new transposed matrix. 
     */
    DynamicMatrix transpose() const {
        DynamicMatrix res(m_cols, m_rows);
        for (int i = 0; i < m_rows; ++i)
            for (int j = 0; j < m_cols; ++j)
                res(j, i) = (*this)(i, j);
        return res;
    }
};

// --- Non-member Operations ---

/** 
 * @brief Matrix-Vector multiplication: y = A * x 
 */
template <typename T>
inline DynamicVector<T> operator*(const DynamicMatrix<T>& A, const DynamicVector<T>& x) {
    CFEM_ASSERT(A.cols() == x.size());
    DynamicVector<T> res(A.rows());
    for (int i = 0; i < A.rows(); ++i) {
        T sum = 0;
        for (int j = 0; j < A.cols(); ++j) {
            sum += A(i, j) * x[j];
        }
        res[i] = sum;
    }
    return res;
}

/** @brief Matrix-Matrix multiplication: C = A * B (ikj optimized) */
template <typename T>
inline DynamicMatrix<T> operator*(const DynamicMatrix<T>& A, const DynamicMatrix<T>& B) {
    CFEM_ASSERT(A.cols() == B.rows());
    DynamicMatrix<T> res(A.rows(), B.cols());
    res.setZero();
    for (int i = 0; i < A.rows(); ++i) {
        for (int k = 0; k < A.cols(); ++k) {
            T val = A(i, k);
            for (int j = 0; j < B.cols(); ++j) {
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
inline DynamicMatrix<T> operator*(T scalar, const DynamicMatrix<T>& m) {
    return m * scalar;
}

/** 
 * @brief Output stream formatting. 
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const DynamicMatrix<T>& m) {
    os << "[size=" << m.rows() << "x" << m.cols() << "]\n";
    for (int i = 0; i < m.rows(); ++i) {
        os << "| ";
        for (int j = 0; j < m.cols(); ++j) {
            os << m(i, j) << (j == m.cols() - 1 ? "" : ", ");
        }
        os << " |\n";
    }
    return os;
}

} //namespace

#endif // CFEM_DYNAMIC_MATRIX