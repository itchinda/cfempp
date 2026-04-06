/*!
 * @file DynamicVector.h
 * @author Itchinda
 * @brief Dynamic-size heap-allocated vector for global FE operations.
 * @version 0.1
 * @date 2026-01-27
 */

#ifndef CFEM_DYNAMIC_VECTOR_
#define CFEM_DYNAMIC_VECTOR_

#include "cfemutils.h"
#include "FixedVector.h"
#include <vector>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <iterator>

namespace cfem {

/*!
 * @class DynamicVector
 * @brief A dynamically-sized vector for general linear algebra.
 * @tparam T The numeric type (e.g., CfemReal, double).
 * @note Use this for operations where size is only known at runtime.
 */
template <typename T>
class DynamicVector {
private:
    std::vector<T> m_data;

public:

    /*! 
    * @brief Default constructor (empty vector). 
    */
    DynamicVector() = default;

    /*! 
     * @brief Size constructor. 
     * @param size Number of elements.
     * @param init Initial value for all elements (default 0).
     */
    explicit DynamicVector(int size, T init = T{}) 
        : m_data(static_cast<size_t>(size), init) {}

    /*!
     * @brief Initializer list constructor. 
     */
    DynamicVector(std::initializer_list<T> list) : m_data(list) {}

    /*! 
     * @brief Conversion constructor from FixedVector.
     *        Allows: DynamicVector<double> dv = someFixedVec;
     */
    template <int N>
    DynamicVector(const FixedVector<T, N>& fixed) {
        m_data.resize(N);
        for(int i = 0; i < N; ++i) m_data[i] = fixed[i];
    }

    // Move Constructor
    DynamicVector(DynamicVector&& other) noexcept 
        : m_data(std::move(other.m_data)) {}

    // Move Assignment
    DynamicVector& operator=(DynamicVector&& other) noexcept {
        if (this != &other) {
            m_data = std::move(other.m_data);
        }
        return *this;
    }

    DynamicVector(const DynamicVector& other) = default;
    DynamicVector& operator=(const DynamicVector& other) = default;

    /*! 
     * @brief Returns the number of elements. 
     */
    inline int size() const noexcept { return static_cast<int>(m_data.size()); }

    /*! 
     * @brief Resizes the vector. 
     */
    void resize(int newSize) { m_data.resize(static_cast<size_t>(newSize)); }

    /**
     * @brief Redimensionne le vecteur et initialise tous les éléments.
     * @param size Nouvelle taille du vecteur.
     * @param value Valeur à affecter à chaque élément.
     */
    void assign(size_t size, const T& value = T()) {
        m_data.assign(size, value);
    }

    /*! 
     * @brief Resets all elements to zero. 
     */
    void setZero() { std::fill(m_data.begin(), m_data.end(), static_cast<T>(0)); }

    T* data() noexcept { return m_data.data(); }
    const T* data() const noexcept { return m_data.data(); }
    void reserve(int n) { m_data.reserve(static_cast<size_t>(n)); }
    void clear() { m_data.clear(); }

    // Extraire une portion du vecteur
    DynamicVector getBlock(int start, int size) const {
        CFEM_ASSERT(start + size <= this->size());
        DynamicVector res(size);
        std::copy(m_data.begin() + start, m_data.begin() + start + size, res.m_data.begin());
        return res;
    }

    void push_back(T value) { m_data.push_back(value); }

    /** @brief Itérateur sur le début des données. */
    auto begin() noexcept { return m_data.begin(); }
    auto end() noexcept { return m_data.end(); }

    /** @brief Itérateur constant sur le début des données. */
    auto begin() const noexcept { return m_data.begin(); }
    auto end() const noexcept { return m_data.end(); }

    /** @brief Itérateur constant explicite (cbegin/cend). */
    auto cbegin() const noexcept { return m_data.cbegin(); }
    auto cend() const noexcept { return m_data.cend(); }

    /*! 
     * @brief Access element (read-write).
     */
    inline T& operator[](int i) {
        CFEM_ASSERT(i >= 0 && i < size());
        return m_data[static_cast<size_t>(i)];
    }

    /*! 
     * @brief Access element (read-only). 
     */
    inline const T& operator[](int i) const {
        CFEM_ASSERT(
            i >= 0 && i < size(),
            (("DynamicVector::operator[]: Index out of bounds. The vector size is "
            + std::to_string(size()) + ", but you are trying to access an element at index "
            + std::to_string(i)).c_str())
        );
        return m_data[static_cast<size_t>(i)];
    }

    // --- Arithmetic Operations ---

    /*!
     * @brief Component-wise addition. 
     */
    DynamicVector operator+(const DynamicVector& other) const {
        CFEM_ASSERT(size() == other.size());
        DynamicVector res(size());
        for (int i = 0; i < size(); ++i) res[i] = m_data[i] + other[i];
        return res;
    }

    /*!
     * @brief Component-wise subtraction. 
     */
    DynamicVector operator-(const DynamicVector& other) const {
        CFEM_ASSERT(size() == other.size());
        DynamicVector res(size());
        for (int i = 0; i < size(); ++i) res[i] = m_data[i] - other[i];
        return res;
    }

    /*! 
     * @brief Scalar multiplication. 
     */
    DynamicVector operator*(T scalar) const {
        DynamicVector res(size());
        for (int i = 0; i < size(); ++i) res[i] = m_data[i] * scalar;
        return res;
    }

    /*! 
     * @brief Scalar division. 
     * @throw CFEM_ASSERT if scalar is near zero.
     */
    DynamicVector operator/(T scalar) const {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        DynamicVector res(size());
        T invScalar = static_cast<T>(1) / scalar;
        for (int i = 0; i < size(); ++i) res[i] = m_data[i] * invScalar;
        return res;
    }

    /*! 
     * @brief In-place addition. 
     */
    DynamicVector& operator+=(const DynamicVector& other) {
        CFEM_ASSERT(size() == other.size());
        for (int i = 0; i < size(); ++i) m_data[i] += other[i];
        return *this;
    }

    /*! 
     * @brief In-place subtraction. 
     */
    DynamicVector& operator-=(const DynamicVector& other) {
        CFEM_ASSERT(size() == other.size());
        for (int i = 0; i < size(); ++i) m_data[i] -= other[i];
        return *this;
    }

    /*! 
     * @brief In-place scalar multiplication. 
     */
    DynamicVector& operator*=(T scalar) {
        for (auto& val : m_data) val *= scalar;
        return *this;
    }

    /*! 
     * @brief In-place scalar division. 
     */
    DynamicVector& operator/=(T scalar) {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        T invScalar = static_cast<T>(1) / scalar;
        for (auto& val : m_data) val *= invScalar;
        return *this;
    }

    // --- Linear Algebra Tools ---

    /*! 
     * @brief Dot product. 
     */
    T dot(const DynamicVector& other) const {
        CFEM_ASSERT(size() == other.size());
        T result = 0;
        for (int i = 0; i < size(); ++i) result += m_data[i] * other[i];
        return result;
    }

    /*! 
     * @brief Returns squared Euclidean norm. 
     * @note Faster than norm() as it avoids std::sqrt.
     */
    T normSquare() const { return dot(*this); }

    /*! 
     * @brief Euclidean norm. 
     */
    T norm() const { return std::sqrt(normSquare()); }

    /*! 
     * @brief Normalizes the vector in-place.
     * @throw CFEM_ASSERT if norm is near zero.
     */
    void normalize() {
        T n = norm();
        CFEM_ASSERT(n > cfem::ZERO_TOLERANCE);
        (*this) *= (static_cast<T>(1) / n);
    }
};

/*! 
 * @brief Left-side scalar multiplication. 
 * @relates FixedVector
 */
template <typename T>
inline DynamicVector<T> operator*(T scalar, const DynamicVector<T>& v) {
    return v * scalar;
}

/*! 
 * @brief Formats the vector for standard output. 
 * @relates FixedVector
 */
template <typename T>
inline std::ostream& operator<<(std::ostream& os, const DynamicVector<T>& v) {
    os << "[size=" << v.size() << "](";
    for (int i = 0; i < v.size(); ++i) {
        os << v[i] << (i == v.size() - 1 ? "" : ", ");
    }
    os << ")";
    return os;
}

}

#endif // CFEM_DYNAMIC_VECTOR