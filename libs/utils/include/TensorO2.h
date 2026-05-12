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
 * @file TensorO2.h
 * @brief Physical tensor definitions for CFEM.
 */

#ifndef CFEM_TENSORO2_H
#define CFEM_TENSORO2_H

#include "FixedMatrix.h"

namespace cfem
{

    /**
     * @brief Storage indices for full second-order tensors.
     */
    namespace idTO2
    {
        constexpr CfemInt xx = 0;
        constexpr CfemInt xy = 1;
        constexpr CfemInt xz = 2;
        constexpr CfemInt yx = 3;
        constexpr CfemInt yy = 4;
        constexpr CfemInt yz = 5;
        constexpr CfemInt zx = 6;
        constexpr CfemInt zy = 7;
        constexpr CfemInt zz = 8;
    }

    /**
     * @brief Storage indices for symmetric second-order tensors.
     */
    namespace idSTO2
    {
        constexpr CfemInt xx = 0;
        constexpr CfemInt yy = 1;
        constexpr CfemInt zz = 2;
        constexpr CfemInt yz = 3;
        constexpr CfemInt xz = 4;
        constexpr CfemInt xy = 5;
    }

    /**
     * @brief Storage indices for skew-symmetric second-order tensors.
     */
    namespace idSSTO2
    {
        constexpr CfemInt yz = 0;
        constexpr CfemInt xz = 1;
        constexpr CfemInt xy = 2;
    }

    /**
     * @brief General second-order tensor (9 components).
     * Used for non-symmetric quantities like the Deformation Gradient (F)
     * or the Velocity Gradient (L).
     */
    using TensorO2 = FixedMatrix<CfemReal, 3, 3>;

    /**
     * @brief Identity tensor constant.
     */
    [[nodiscard]] constexpr inline TensorO2 IdentityTensor() noexcept
    {
        TensorO2 I;
        I.setIdentity();
        return I;
    }

    /**
     * @brief Stores a symmetric 3x3 tensor using only 6 values.
     * Storage order (Voigt): [0]=xx, [1]=yy, [2]=zz, [3]=yz, [4]=xz, [5]=xy
     */
    struct SymTensorO2
    {
        CfemReal data[6]{};

    private:
        static constexpr CfemInt map[3][3] = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}};

    public:
        constexpr SymTensorO2() noexcept = default;

        constexpr SymTensorO2(std::initializer_list<CfemReal> list)
        {
            CFEM_ASSERT(list.size() == 6);
            CfemInt i = 0;
            for (auto val : list)
                data[i++] = val;
        }

        /** @brief Access by component index (0 to 5). */
        [[nodiscard]] constexpr inline CfemReal &operator[](CfemInt i)
        {
            CFEM_ASSERT(i >= 0 && i < 6);
            return data[i];
        }

        [[nodiscard]] constexpr inline const CfemReal &operator[](CfemInt i) const
        {
            CFEM_ASSERT(i >= 0 && i < 6);
            return data[i];
        }

        /**
         * @brief Map 3x3 indices to the symmetric storage.
         *
         * @param i Row index
         * @param j Column index
         * @return CfemReal
         */
        [[nodiscard]] constexpr inline CfemReal operator()(CfemInt i, CfemInt j) const
        {
            CFEM_ASSERT(i >= 0 && i < 3 && j >= 0 && j < 3);
            return data[map[i][j]];
        }

        /**
         * @brief Negate operator
         *
         * @return SymTensorO2
         */
        [[nodiscard]] constexpr SymTensorO2 operator-() const noexcept
        {
            return {-data[idSTO2::xx], -data[idSTO2::yy], -data[idSTO2::zz],
                    -data[idSTO2::yz], -data[idSTO2::xz], -data[idSTO2::xy]};
        }

        [[nodiscard]] constexpr inline CfemReal xx() const noexcept { return data[idSTO2::xx]; }
        [[nodiscard]] constexpr inline CfemReal yy() const noexcept { return data[idSTO2::yy]; }
        [[nodiscard]] constexpr inline CfemReal zz() const noexcept { return data[idSTO2::zz]; }
        [[nodiscard]] constexpr inline CfemReal yz() const noexcept { return data[idSTO2::yz]; }
        [[nodiscard]] constexpr inline CfemReal xz() const noexcept { return data[idSTO2::xz]; }
        [[nodiscard]] constexpr inline CfemReal xy() const noexcept { return data[idSTO2::xy]; }

        [[nodiscard]] constexpr inline CfemReal &xx() noexcept { return data[idSTO2::xx]; }
        [[nodiscard]] constexpr inline CfemReal &yy() noexcept { return data[idSTO2::yy]; }
        [[nodiscard]] constexpr inline CfemReal &zz() noexcept { return data[idSTO2::zz]; }
        [[nodiscard]] constexpr inline CfemReal &yz() noexcept { return data[idSTO2::yz]; }
        [[nodiscard]] constexpr inline CfemReal &xz() noexcept { return data[idSTO2::xz]; }
        [[nodiscard]] constexpr inline CfemReal &xy() noexcept { return data[idSTO2::xy]; }

        // Basic Math Ops

        [[nodiscard]] static constexpr SymTensorO2 Identity() noexcept
        {
            return {1, 1, 1, 0, 0, 0};
        }

        [[nodiscard]] static constexpr SymTensorO2 Zero() noexcept
        {
            return {0, 0, 0, 0, 0, 0};
        }

        constexpr void setZero() noexcept
        {
            for (CfemInt i = 0; i < 6; ++i)
                data[i] = 0;
        }

        constexpr void setTo(CfemReal value) noexcept
        {
            for (CfemInt i = 0; i < 6; ++i)
                data[i] = value;
        }

        constexpr inline void toMatrix(CfemReal M[3][3]) const noexcept
        {
            M[0][0] = data[idSTO2::xx];
            M[0][1] = data[idSTO2::xy];
            M[0][2] = data[idSTO2::xz];
            M[1][0] = data[idSTO2::xy];
            M[1][1] = data[idSTO2::yy];
            M[1][2] = data[idSTO2::yz];
            M[2][0] = data[idSTO2::xz];
            M[2][1] = data[idSTO2::yz];
            M[2][2] = data[idSTO2::zz];
        }

        [[nodiscard]] constexpr SymTensorO2 operator+(const SymTensorO2 &other) const noexcept
        {
            return {data[idSTO2::xx] + other.data[idSTO2::xx], data[idSTO2::yy] + other.data[idSTO2::yy], data[idSTO2::zz] + other.data[idSTO2::zz],
                    data[idSTO2::yz] + other.data[idSTO2::yz], data[idSTO2::xz] + other.data[idSTO2::xz], data[idSTO2::xy] + other.data[idSTO2::xy]};
        }

        [[nodiscard]] constexpr SymTensorO2 operator-(const SymTensorO2 &other) const noexcept
        {
            return {data[idSTO2::xx] - other.data[idSTO2::xx], data[idSTO2::yy] - other.data[idSTO2::yy], data[idSTO2::zz] - other.data[idSTO2::zz],
                    data[idSTO2::yz] - other.data[idSTO2::yz], data[idSTO2::xz] - other.data[idSTO2::xz], data[idSTO2::xy] - other.data[idSTO2::xy]};
        }

        [[nodiscard]] constexpr SymTensorO2 operator*(CfemReal scalar) const noexcept
        {
            return {data[idSTO2::xx] * scalar, data[idSTO2::yy] * scalar, data[idSTO2::zz] * scalar,
                    data[idSTO2::yz] * scalar, data[idSTO2::xz] * scalar, data[idSTO2::xy] * scalar};
        }

        [[nodiscard]] constexpr SymTensorO2 operator/(CfemReal scalar) const noexcept
        {
            if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
            {
                CFEM_WARN << "Division by zero in SymTensorO2::operator/";
            }
            return {data[idSTO2::xx] / scalar, data[idSTO2::yy] / scalar, data[idSTO2::zz] / scalar,
                    data[idSTO2::yz] / scalar, data[idSTO2::xz] / scalar, data[idSTO2::xy] / scalar};
        }

        [[nodiscard]] constexpr inline Vector3D operator*(const Vector3D &v) const noexcept
        {
            return {
                data[idSTO2::xx] * v.x + data[idSTO2::xy] * v.y + data[idSTO2::xz] * v.z, // Row x: [xx, xy, xz]
                data[idSTO2::xy] * v.x + data[idSTO2::yy] * v.y + data[idSTO2::yz] * v.z, // Row y: [yx, yy, yz]
                data[idSTO2::xz] * v.x + data[idSTO2::yz] * v.y + data[idSTO2::zz] * v.z  // Row z: [zx, zy, zz]
            };
        }

        constexpr SymTensorO2 &operator+=(const SymTensorO2 &other) noexcept
        {
            data[idSTO2::xx] += other.data[idSTO2::xx];
            data[idSTO2::yy] += other.data[idSTO2::yy];
            data[idSTO2::zz] += other.data[idSTO2::zz];
            data[idSTO2::yz] += other.data[idSTO2::yz];
            data[idSTO2::xz] += other.data[idSTO2::xz];
            data[idSTO2::xy] += other.data[idSTO2::xy];
            return *this;
        }

        constexpr SymTensorO2 &operator-=(const SymTensorO2 &other) noexcept
        {
            data[idSTO2::xx] -= other.data[idSTO2::xx];
            data[idSTO2::yy] -= other.data[idSTO2::yy];
            data[idSTO2::zz] -= other.data[idSTO2::zz];
            data[idSTO2::yz] -= other.data[idSTO2::yz];
            data[idSTO2::xz] -= other.data[idSTO2::xz];
            data[idSTO2::xy] -= other.data[idSTO2::xy];
            return *this;
        }

        constexpr SymTensorO2 &operator*=(CfemReal scalar) noexcept
        {
            data[idSTO2::xx] *= scalar;
            data[idSTO2::yy] *= scalar;
            data[idSTO2::zz] *= scalar;
            data[idSTO2::yz] *= scalar;
            data[idSTO2::xz] *= scalar;
            data[idSTO2::xy] *= scalar;
            return *this;
        }

        constexpr SymTensorO2 &operator/=(CfemReal scalar) noexcept
        {
            if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
            {
                CFEM_WARN << "Division by zero in SymTensorO2::operator/=";
            }
            data[idSTO2::xx] /= scalar;
            data[idSTO2::yy] /= scalar;
            data[idSTO2::zz] /= scalar;
            data[idSTO2::yz] /= scalar;
            data[idSTO2::xz] /= scalar;
            data[idSTO2::xy] /= scalar;
            return *this;
        }

        [[nodiscard]] constexpr bool operator==(const SymTensorO2 &other) const noexcept
        {
            for (CfemInt i = 0; i < 6; ++i)
                if (data[i] != other[i])
                    return false;
            return true;
        }

        [[nodiscard]] constexpr bool operator!=(const SymTensorO2 &other) const noexcept
        {
            return !(*this == other);
        }

        /**
         * @brief Computes the quadratic form M^T * S * M.
         *
         * @param M Transformation matrix.
         * @return Transformed symmetric tensor.
         *
         * @details
         * This operation preserves symmetry and is commonly used for second-order tensor transformations such as:
         *
         * H_phys = J^{-T} * H_{ref} * J^{-1}
         *
         * where H denotes a Hessian tensor.
         */
        [[nodiscard]] constexpr inline SymTensorO2 quadform(const TensorO2 &M) const noexcept
        {
            SymTensorO2 res;

            // On calcule d'abord le produit intermédiaire T = (*this) * M
            // Rappel stockage : [0]=xx, [1]=yy, [2]=zz, [3]=yz, [4]=xz, [5]=xy
            TensorO2 T;

            // Ligne 0 de (*this) : [xx, xy, xz]
            const CfemReal m00 = M.data[idTO2::xx], m01 = M.data[idTO2::xy], m02 = M.data[idTO2::xz];
            const CfemReal m10 = M.data[idTO2::yx], m11 = M.data[idTO2::yy], m12 = M.data[idTO2::yz];
            const CfemReal m20 = M.data[idTO2::zx], m21 = M.data[idTO2::zy], m22 = M.data[idTO2::zz];

            T.data[idTO2::xx] = data[idSTO2::xx] * m00 + data[idSTO2::xy] * m10 + data[idSTO2::xz] * m20;
            T.data[idTO2::xy] = data[idSTO2::xx] * m01 + data[idSTO2::xy] * m11 + data[idSTO2::xz] * m21;
            T.data[idTO2::xz] = data[idSTO2::xx] * m02 + data[idSTO2::xy] * m12 + data[idSTO2::xz] * m22;

            // Ligne 1 de (*this) : [xy, yy, yz]
            T.data[idTO2::yx] = data[idSTO2::xy] * m00 + data[idSTO2::yy] * m10 + data[idSTO2::yz] * m20;
            T.data[idTO2::yy] = data[idSTO2::xy] * m01 + data[idSTO2::yy] * m11 + data[idSTO2::yz] * m21;
            T.data[idTO2::yz] = data[idSTO2::xy] * m02 + data[idSTO2::yy] * m12 + data[idSTO2::yz] * m22;

            // Ligne 2 de (*this) : [xz, yz, zz]
            T.data[idTO2::zx] = data[idSTO2::xz] * m00 + data[idSTO2::yz] * m10 + data[idSTO2::zz] * m20;
            T.data[idTO2::zy] = data[idSTO2::xz] * m01 + data[idSTO2::yz] * m11 + data[idSTO2::zz] * m21;
            T.data[idTO2::zz] = data[idSTO2::xz] * m02 + data[idSTO2::yz] * m12 + data[idSTO2::zz] * m22;

            // On calcule res = M^T * T
            // On ne calcule que les 6 composantes symétriques
            res[0] = m00 * T.data[idTO2::xx] + m10 * T.data[idTO2::yx] + m20 * T.data[idTO2::zx]; // xx
            res[1] = m01 * T.data[idTO2::xy] + m11 * T.data[idTO2::yy] + m21 * T.data[idTO2::zy]; // yy
            res[2] = m02 * T.data[idTO2::xz] + m12 * T.data[idTO2::yz] + m22 * T.data[idTO2::zz]; // zz

            res[3] = m01 * T.data[idTO2::xz] + m11 * T.data[idTO2::yz] + m21 * T.data[idTO2::zz]; // yz
            res[4] = m00 * T.data[idTO2::xz] + m10 * T.data[idTO2::yz] + m20 * T.data[idTO2::zz]; // xz
            res[5] = m00 * T.data[idTO2::xy] + m10 * T.data[idTO2::yy] + m20 * T.data[idTO2::zy]; // xy

            return res;
        }

        // Physical Invariants

        /** @brief Trace of the tensor (sum of diagonals). */
        [[nodiscard]] constexpr inline CfemReal trace() const noexcept
        {
            return data[idSTO2::xx] + data[idSTO2::yy] + data[idSTO2::zz];
        }

        /** @brief Mean (hydrostatic) component: p = tr(sigma)/3 */
        [[nodiscard]] constexpr inline CfemReal mean() const noexcept
        {
            return trace() / 3.0;
        }

        /** @brief Returns the deviatoric part: s = sigma - p*I */
        [[nodiscard]] constexpr SymTensorO2 deviator() const noexcept
        {
            CfemReal p = mean();
            return {data[idSTO2::xx] - p, data[idSTO2::yy] - p, data[idSTO2::zz] - p, data[idSTO2::yz], data[idSTO2::xz], data[idSTO2::xy]};
        }

        /** @brief Double contraction with another symmetric tensor (A : B). */
        [[nodiscard]] constexpr inline CfemReal dotDot(const SymTensorO2 &other) const noexcept
        {
            // xx*xx + yy*yy + zz*zz + 2*(yz*yz + xz*xz + xy*xy)
            return data[idSTO2::xx] * other.data[idSTO2::xx] + data[idSTO2::yy] * other.data[idSTO2::yy] + data[idSTO2::zz] * other.data[idSTO2::zz] +
                   2.0 * (data[idSTO2::yz] * other.data[idSTO2::yz] + data[idSTO2::xz] * other.data[idSTO2::xz] + data[idSTO2::xy] * other.data[idSTO2::xy]);
        }

        /** @brief Second invariant of deviatoric tensor: J2 = 0.5 * (s : s) */
        [[nodiscard]] constexpr inline CfemReal J2() const noexcept
        {
            SymTensorO2 s = deviator();
            return 0.5 * s.dotDot(s);
        }

        /** @brief Computes Von Mises equivalent stress (if this is a stress tensor). */
        [[nodiscard]] inline CfemReal vonMises() const noexcept
        {
            CfemReal t = trace() / 3.0;
            CfemReal s_xx = data[idSTO2::xx] - t;
            CfemReal s_yy = data[idSTO2::yy] - t;
            CfemReal s_zz = data[idSTO2::zz] - t;
            return std::sqrt(1.5 * (s_xx * s_xx + s_yy * s_yy + s_zz * s_zz +
                                    2.0 * (data[idSTO2::yz] * data[idSTO2::yz] + data[idSTO2::xz] * data[idSTO2::xz] + data[idSTO2::xy] * data[idSTO2::xy])));
        }

        /** 
         * @brief Norm of the tensor (Frobenius norm) 
         */
        [[nodiscard]] inline CfemReal norm() const noexcept
        {
            return std::sqrt(this->dotDot(*this));
        }

        /**
         * @brief Compute the transpose of the matrix
         * 
         * @return constexpr SymTensorO2 The transposed matrix
         * 
         * @note the transpose of a symmetric matrix is exactely equal to the matrix itself
         */
        [[nodiscard]] constexpr inline SymTensorO2 transpose() const noexcept
        {
            return *this;
        }

        /**
         * @brief Computes the double contraction between this tensor and a general tensor.
         *
         * @param T General second-order tensor.
         * @return Scalar value defined as:
         *
         * S : T = \sum_{ij} S_{ij} T_{ij}
         */
        [[nodiscard]] constexpr inline CfemReal doubleDot(const TensorO2& T) const noexcept
        {
            // S : T = S_ij T_ij
            return data[idSTO2::xx] * T.data[idTO2::xx] + 
                data[idSTO2::yy] * T.data[idTO2::yy] + 
                data[idSTO2::zz] * T.data[idTO2::zz] + 
                data[idSTO2::xy] * (T.data[idTO2::xy] + T.data[idTO2::yx]) + 
                data[idSTO2::xz] * (T.data[idTO2::xz] + T.data[idTO2::zx]) + 
                data[idSTO2::yz] * (T.data[idTO2::yz] + T.data[idTO2::zy]);
        }

        [[nodiscard]] constexpr inline SymTensorO2 deviatoric() const noexcept
        {
            const CfemReal p = (data[idSTO2::xx] + data[idSTO2::yy] + data[idSTO2::zz]) / 3.0;
            SymTensorO2 res = *this;
            res.data[idSTO2::xx] -= p;
            res.data[idSTO2::yy] -= p;
            res.data[idSTO2::zz] -= p;
            return res;
        }

        /**
         * @brief Computes the outer product in Voigt notation.
         *
         * @param other Second symmetric tensor.
         * @return 6x6 matrix such that:
         *
         * C(i,j) = this[i] * other[j]
         *
         * @details
         * Commonly used to build tangent constitutive operators.
         */
        [[nodiscard]] constexpr inline FixedMatrix<CfemReal, 6, 6> outer(const SymTensorO2& other) const noexcept
        {
            FixedMatrix<CfemReal, 6, 6> res;
            for(CfemInt i = 0; i < 6; ++i) {
                for(CfemInt j = 0; j < 6; ++j) {
                    res(i, j) = this->data[i] * other.data[j];
                }
            }
            return res;
        }

        
    };



    // ------------------------------------------------------------------------------- //
    // ------------------------------ SKEW-SYMMETRIC --------------------------------- //
    // ------------------------------------------------------------------------------- //

    /**
     * @brief Skew-symmetric second-order tensor.
     * 
     * @details
     * Stores only the 3 independent off-diagonal components.
     * By definition, the diagonal is exactly zero and W_ij = -W_ji.
     * Storage order matches the off-diagonal parts of SymTensorO2:
     * [0]=yz, [1]=xz, [2]=xy
     */
    struct SkewSymTensorO2
    {
        CfemReal data[3]{};

    public:
        constexpr SkewSymTensorO2() noexcept = default;

        constexpr SkewSymTensorO2(std::initializer_list<CfemReal> list)
        {
            CFEM_ASSERT(list.size() == 3);
            CfemInt i = 0;
            for (auto val : list)
                data[i++] = val;
        }

        /** @brief Access by internal storage index (0 to 2). */
        [[nodiscard]] constexpr inline CfemReal &operator[](CfemInt i)
        {
            CFEM_ASSERT(i >= 0 && i < 3);
            return data[i];
        }

        [[nodiscard]] constexpr inline const CfemReal &operator[](CfemInt i) const
        {
            CFEM_ASSERT(i >= 0 && i < 3);
            return data[i];
        }

        /**
         * @brief Map 3x3 indices to the skew-symmetric storage.
         * Returns 0 for diagonal terms, and handles signs automatically.
         */
        [[nodiscard]] constexpr inline CfemReal operator()(CfemInt i, CfemInt j) const
        {
            CFEM_ASSERT(i >= 0 && i < 3 && j >= 0 && j < 3, "Index out of bounds in SkewSymTensorO2::operator()");

            switch (3*i + j)
            {
                case 1: return  data[idSSTO2::xy]; // (0,1)
                case 2: return  data[idSSTO2::xz]; // (0,2)

                case 3: return -data[idSSTO2::xy]; // (1,0)
                case 5: return  data[idSSTO2::yz]; // (1,2)

                case 6: return -data[idSSTO2::xz]; // (2,0)
                case 7: return -data[idSSTO2::yz]; // (2,1)

                default: return 0.0; // diagonal
            }
        }

        /** @brief Negate operator */
        [[nodiscard]] constexpr SkewSymTensorO2 operator-() const noexcept
        {
            return {-data[idSSTO2::yz], -data[idSSTO2::xz], -data[idSSTO2::xy]};
        }

        [[nodiscard]] constexpr inline CfemReal yz() const noexcept { return data[idSSTO2::yz]; }
        [[nodiscard]] constexpr inline CfemReal xz() const noexcept { return data[idSSTO2::xz]; }
        [[nodiscard]] constexpr inline CfemReal xy() const noexcept { return data[idSSTO2::xy]; }

        [[nodiscard]] constexpr inline CfemReal &yz() noexcept { return data[idSSTO2::yz]; }
        [[nodiscard]] constexpr inline CfemReal &xz() noexcept { return data[idSSTO2::xz]; }
        [[nodiscard]] constexpr inline CfemReal &xy() noexcept { return data[idSSTO2::xy]; }

        // Basic Math Ops

        [[nodiscard]] static constexpr SkewSymTensorO2 Zero() noexcept
        {
            return {0, 0, 0};
        }

        constexpr void setZero() noexcept
        {
            data[idSSTO2::yz] = 0;
            data[idSSTO2::xz] = 0;
            data[idSSTO2::xy] = 0;
        }

        constexpr void setTo(CfemReal value) noexcept
        {
            data[idSSTO2::yz] = value;
            data[idSSTO2::xz] = value;
            data[idSSTO2::xy] = value;
        }

        constexpr inline void fromMatrix(CfemReal M[3][3]) noexcept
        {
            // We assume M is already skew-symmetric, we extract the upper triangle
            data[idSSTO2::yz] = M[1][2]; // yz
            data[idSSTO2::xz] = M[0][2]; // xz
            data[idSSTO2::xy] = M[0][1]; // xy
        }

        [[nodiscard]] constexpr SkewSymTensorO2 operator+(const SkewSymTensorO2 &other) const noexcept
        {
            return {data[idSSTO2::yz] + other.data[idSSTO2::yz], data[idSSTO2::xz] + other.data[idSSTO2::xz], data[idSSTO2::xy] + other.data[idSSTO2::xy]};
        }

        [[nodiscard]] constexpr SkewSymTensorO2 operator-(const SkewSymTensorO2 &other) const noexcept
        {
            return {data[idSSTO2::yz] - other.data[idSSTO2::yz], data[idSSTO2::xz] - other.data[idSSTO2::xz], data[idSSTO2::xy] - other.data[idSSTO2::xy]};
        }

        [[nodiscard]] constexpr SkewSymTensorO2 operator*(CfemReal scalar) const noexcept
        {
            return {data[idSSTO2::yz] * scalar, data[idSSTO2::xz] * scalar, data[idSSTO2::xy] * scalar};
        }

        [[nodiscard]] constexpr SkewSymTensorO2 operator/(CfemReal scalar) const noexcept
        {
            if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
            {
                CFEM_WARN << "Division by zero in SkewSymTensorO2::operator/";
            }
            return {data[idSSTO2::yz] / scalar, data[idSSTO2::xz] / scalar, data[idSSTO2::xy] / scalar};
        }

        /**
         * @brief Computes W * v (Equivalent to the cross product: axial_vector x v)
         */
        [[nodiscard]] constexpr inline Vector3D operator*(const Vector3D &v) const noexcept
        {
            return {
                data[idSSTO2::xy] * v.y + data[idSSTO2::xz] * v.z,  // Row x: [0, xy, xz]
                -data[idSSTO2::xy] * v.x + data[idSSTO2::yz] * v.z, // Row y: [-xy, 0, yz]
                -data[idSSTO2::xz] * v.x - data[idSSTO2::yz] * v.y  // Row z: [-xz, -yz, 0]
            };
        }

        constexpr SkewSymTensorO2 &operator+=(const SkewSymTensorO2 &other) noexcept
        {
            data[idSSTO2::yz] += other.data[idSSTO2::yz];
            data[idSSTO2::xz] += other.data[idSSTO2::xz];
            data[idSSTO2::xy] += other.data[idSSTO2::xy];
            return *this;
        }

        constexpr SkewSymTensorO2 &operator-=(const SkewSymTensorO2 &other) noexcept
        {
            data[idSSTO2::yz] -= other.data[idSSTO2::yz];
            data[idSSTO2::xz] -= other.data[idSSTO2::xz];
            data[idSSTO2::xy] -= other.data[idSSTO2::xy];
            return *this;
        }

        constexpr SkewSymTensorO2 &operator*=(CfemReal scalar) noexcept
        {
            data[idSSTO2::yz] *= scalar;
            data[idSSTO2::xz] *= scalar;
            data[idSSTO2::xy] *= scalar;
            return *this;
        }

        constexpr SkewSymTensorO2 &operator/=(CfemReal scalar) noexcept
        {
            if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
            {
                CFEM_WARN << "Division by zero in SkewSymTensorO2::operator/=";
            }
            data[idSSTO2::yz] /= scalar;
            data[idSSTO2::xz] /= scalar;
            data[idSSTO2::xy] /= scalar;
            return *this;
        }

        [[nodiscard]] constexpr bool operator==(const SkewSymTensorO2 &other) const noexcept
        {
            for (CfemInt i = 0; i < 3; ++i)
                if (data[i] != other.data[i])
                    return false;
            return true;
        }

        [[nodiscard]] constexpr bool operator!=(const SkewSymTensorO2 &other) const noexcept
        {
            return !(*this == other);
        }

        /**
         * @brief Computes the quadratic form Mᵀ * S * M.
         *
         * @param M Transformation matrix.
         * @return Transformed skew-symmetric tensor.
         *
         * @details
         * This operation preserves skew-symmetry.
         */
        [[nodiscard]] constexpr inline SkewSymTensorO2 quadform(const TensorO2 &M) const noexcept
        {
            SkewSymTensorO2 res;

            // 1. T = (*this) * M
            TensorO2 T;
            const CfemReal m00 = M.data[idTO2::xx], m01 = M.data[idTO2::xy], m02 = M.data[idTO2::xz];
            const CfemReal m10 = M.data[idTO2::yx], m11 = M.data[idTO2::yy], m12 = M.data[idTO2::yz];
            const CfemReal m20 = M.data[idTO2::zx], m21 = M.data[idTO2::zy], m22 = M.data[idTO2::zz];

            // Row 0 of W: [0, xy, xz]
            T.data[idTO2::xx] = data[idSSTO2::xy] * m10 + data[idSSTO2::xz] * m20;
            T.data[idTO2::xy] = data[idSSTO2::xy] * m11 + data[idSSTO2::xz] * m21;
            T.data[idTO2::xz] = data[idSSTO2::xy] * m12 + data[idSSTO2::xz] * m22;

            // Row 1 of W: [-xy, 0, yz]
            T.data[idTO2::yx] = -data[idSSTO2::xy] * m00 + data[idSSTO2::yz] * m20;
            T.data[idTO2::yy] = -data[idSSTO2::xy] * m01 + data[idSSTO2::yz] * m21;
            T.data[idTO2::yz] = -data[idSSTO2::xy] * m02 + data[idSSTO2::yz] * m22;

            // Row 2 of W: [-xz, -yz, 0]
            T.data[idTO2::zx] = -data[idSSTO2::xz] * m00 - data[idSSTO2::yz] * m10;
            T.data[idTO2::zy] = -data[idSSTO2::xz] * m01 - data[idSSTO2::yz] * m11;
            T.data[idTO2::zz] = -data[idSSTO2::xz] * m02 - data[idSSTO2::yz] * m12;

            // 2. res = M^T * T (Only extract the 3 skew components)
            res.data[idSSTO2::yz] = m01 * T.data[idTO2::xz] + m11 * T.data[idTO2::yz] + m21 * T.data[idTO2::zz]; // yz
            res.data[idSSTO2::xz] = m00 * T.data[idTO2::xz] + m10 * T.data[idTO2::yz] + m20 * T.data[idTO2::zz]; // xz
            res.data[idSSTO2::xy] = m00 * T.data[idTO2::xy] + m10 * T.data[idTO2::yy] + m20 * T.data[idTO2::zy]; // xy

            return res;
        }

        // Physical Invariants

        /** @brief Trace is always exactly zero for a skew-symmetric tensor. */
        [[nodiscard]] constexpr inline CfemReal trace() const noexcept
        {
            return 0.0;
        }

        /** @brief Double contraction with another skew-symmetric tensor. */
        [[nodiscard]] constexpr inline CfemReal dotDot(const SkewSymTensorO2 &other) const noexcept
        {
            // A:B = 2 * (A_xy*B_xy + A_xz*B_xz + A_yz*B_yz)
            return 2.0 * (data[idSSTO2::yz] * other.data[idSSTO2::yz] + data[idSSTO2::xz] * other.data[idSSTO2::xz] + data[idSSTO2::xy] * other.data[idSSTO2::xy]);
        }

        /** @brief Frobenius norm */
        [[nodiscard]] inline CfemReal norm() const noexcept
        {
            return std::sqrt(this->dotDot(*this));
        }


        /**
         * @brief Returns the axial vector associated with the skew-symmetric tensor.
         *
         * @details
         * For a tensor W:
         *
         *      [  0  -wz   wy ]
         * W =  [ wz    0  -wx ]
         *      [-wy   wx    0 ]
         *
         * this function returns the vector (wx, wy, wz).
         */
        [[nodiscard]] constexpr inline FixedVector<CfemReal, 3> axialVector() const noexcept
        {
            FixedVector<CfemReal, 3> w;
            w.data[0] = -data[idSSTO2::yz];
            w.data[1] =  data[idSSTO2::xz];
            w.data[2] = -data[idSSTO2::xy];
            return w;
        }

        /**
         * @brief Builds a skew-symmetric tensor from its axial vector.
         *
         * @param w Axial vector.
         * @return Associated skew-symmetric tensor.
         */
        [[nodiscard]] static constexpr inline SkewSymTensorO2 fromAxialVector(const FixedVector<CfemReal, 3>& w) noexcept
        {
            SkewSymTensorO2 W;
            W.data[idSSTO2::xy] = -w.data[2];
            W.data[idSSTO2::xz] =  w.data[1];
            W.data[idSSTO2::yz] = -w.data[0];
            return W;
        }

        [[nodiscard]] constexpr inline SkewSymTensorO2 transpose() const noexcept
        {
            SkewSymTensorO2 W_T;
            W_T.data[idSSTO2::xy] = -data[idSSTO2::xy];
            W_T.data[idSSTO2::xz] = -data[idSSTO2::xz];
            W_T.data[idSSTO2::yz] = -data[idSSTO2::yz];
            return W_T;
        }

        /**
         * @brief Computes the matrix exponential using Rodrigues' formula.
         *
         * @return Proper orthogonal rotation tensor.
         *
         * @details
         * For a skew-symmetric tensor W associated with an axial vector \omega, this computes:
         *
         * R = exp(W)
         *
         * which belongs to SO(3).
         */
        [[nodiscard]] inline TensorO2 exp() const noexcept
        {
            FixedVector<CfemReal, 3> w = axialVector();
            const CfemReal theta = std::sqrt(w.data[0]*w.data[0] + 
                                             w.data[1]*w.data[1] + 
                                             w.data[2]*w.data[2] );
            
            TensorO2 res = IdentityTensor();
            
            if (theta < ZERO_TOLERANCE) {
                // Taylor expansion for small angles
                res(0,1) += data[idSSTO2::xy]; res(0,2) += data[idSSTO2::xz];
                res(1,0) -= data[idSSTO2::xy]; res(1,2) += data[idSSTO2::yz];
                res(2,0) -= data[idSSTO2::xz]; res(2,1) -= data[idSSTO2::yz];
                return res;
            }

            const CfemReal coef1 = std::sin(theta) / theta;
            const CfemReal coef2 = (1.0 - std::cos(theta)) / (theta * theta);

            // Optimisation : compute R = I + coef1 * W + coef2 * W^2 component by component
            res(0,0) += coef2 * (-w.data[1]*w.data[1] - w.data[2]*w.data[2]);
            res(0,1) += coef1 * data[idSSTO2::xy] + coef2 * (w.data[0]*w.data[1]);
            res(0,2) += coef1 * data[idSSTO2::xz] + coef2 * (w.data[0]*w.data[2]);

            res(1,0) -= coef1 * data[idSSTO2::xy] - coef2 * (w.data[0]*w.data[1]);
            res(1,1) += coef2 * (-w.data[0]*w.data[0] - w.data[2]*w.data[2]);
            res(1,2) += coef1 * data[idSSTO2::yz] + coef2 * (w.data[1]*w.data[2]);

            res(2,0) -= coef1 * data[idSSTO2::xz] - coef2 * (w.data[0]*w.data[2]);
            res(2,1) -= coef1 * data[idSSTO2::yz] - coef2 * (w.data[1]*w.data[2]);
            res(2,2) += coef2 * (-w.data[0]*w.data[0] - w.data[1]*w.data[1]);

            return res;
        }

        /**
         * @brief Returns the magnitude of the imaginary eigen values (||\omega||).
         * 
         * @note The eigen values of a 2nd order skew-symmetric tensor are either {0, i||w||, -i||w||} , where ||w|| is the magnitude of the axial vector.
         * This means ||w|| also corresponds to the magnitude of the pure imaginary eigen values...
         */
        [[nodiscard]] inline CfemReal imaginaryEigenvalueMagnitude() const noexcept
        {
            return std::sqrt( xy()*xy() + xz()*xz() + yz()*yz() );
        }
        
    };

    // Global Ops

    [[nodiscard]] constexpr inline SymTensorO2 operator*(CfemReal scalar, const SymTensorO2 &S) noexcept
    {
        return S * scalar;
    }

    [[nodiscard]] constexpr inline TensorO2 operator*(const TensorO2 &T, const SymTensorO2 &S) noexcept
    {
        TensorO2 res;
        const CfemReal a00 = T.data[idTO2::xx], a01 = T.data[idTO2::xy], a02 = T.data[idTO2::xz];
        const CfemReal a10 = T.data[idTO2::yx], a11 = T.data[idTO2::yy], a12 = T.data[idTO2::yz];
        const CfemReal a20 = T.data[idTO2::zx], a21 = T.data[idTO2::zy], a22 = T.data[idTO2::zz];

        res.data[idTO2::xx] = a00 * S.xx() + a01 * S.xy() + a02 * S.xz();
        res.data[idTO2::xy] = a00 * S.xy() + a01 * S.yy() + a02 * S.yz();
        res.data[idTO2::xz] = a00 * S.xz() + a01 * S.yz() + a02 * S.zz();

        res.data[idTO2::yx] = a10 * S.xx() + a11 * S.xy() + a12 * S.xz();
        res.data[idTO2::yy] = a10 * S.xy() + a11 * S.yy() + a12 * S.yz();
        res.data[idTO2::yz] = a10 * S.xz() + a11 * S.yz() + a12 * S.zz();

        res.data[idTO2::zx] = a20 * S.xx() + a21 * S.xy() + a22 * S.xz();
        res.data[idTO2::zy] = a20 * S.xy() + a21 * S.yy() + a22 * S.yz();
        res.data[idTO2::zz] = a20 * S.xz() + a21 * S.yz() + a22 * S.zz();

        return res;
    }

    [[nodiscard]] constexpr inline TensorO2 operator*( const SymTensorO2 &S,
                                                       const TensorO2 &T ) noexcept
    {
        TensorO2 res;
        const CfemReal a00 = T.data[idTO2::xx], a01 = T.data[idTO2::xy], a02 = T.data[idTO2::xz];
        const CfemReal a10 = T.data[idTO2::yx], a11 = T.data[idTO2::yy], a12 = T.data[idTO2::yz];
        const CfemReal a20 = T.data[idTO2::zx], a21 = T.data[idTO2::zy], a22 = T.data[idTO2::zz];

        res.data[idTO2::xx] = S.xx() * a00 + S.xy() * a10 + S.xz() * a20;
        res.data[idTO2::xy] = S.xx() * a01 + S.xy() * a11 + S.xz() * a21;
        res.data[idTO2::xz] = S.xx() * a02 + S.xy() * a12 + S.xz() * a22;

        res.data[idTO2::yx] = S.xy() * a00 + S.yy() * a10 + S.yz() * a20;
        res.data[idTO2::yy] = S.xy() * a01 + S.yy() * a11 + S.yz() * a21;
        res.data[idTO2::yz] = S.xy() * a02 + S.yy() * a12 + S.yz() * a22;

        res.data[idTO2::zx] = S.xz() * a00 + S.yz() * a10 + S.zz() * a20;
        res.data[idTO2::zy] = S.xz() * a01 + S.yz() * a11 + S.zz() * a21;
        res.data[idTO2::zz] = S.xz() * a02 + S.yz() * a12 + S.zz() * a22;

        return res;
    }

    [[nodiscard]] constexpr inline TensorO2 operator*( const SymTensorO2 &A,
                                                       const SymTensorO2 &B ) noexcept
    {
        TensorO2 res;

        // Ligne 0 de A : [xx, xy, xz] par colonnes de B
        res.data[idTO2::xx] = A.xx() * B.xx() + A.xy() * B.xy() + A.xz() * B.xz();
        res.data[idTO2::xy] = A.xx() * B.xy() + A.xy() * B.yy() + A.xz() * B.yz();
        res.data[idTO2::xz] = A.xx() * B.xz() + A.xy() * B.yz() + A.xz() * B.zz();

        // Ligne 1 de A : [xy, yy, yz] par colonnes de B
        res.data[idTO2::yx] = A.xy() * B.xx() + A.yy() * B.xy() + A.yz() * B.xz();
        res.data[idTO2::yy] = A.xy() * B.xy() + A.yy() * B.yy() + A.yz() * B.yz();
        res.data[idTO2::yz] = A.xy() * B.xz() + A.yy() * B.yz() + A.yz() * B.zz();

        // Ligne 2 de A : [xz, yz, zz] par colonnes de B
        res.data[idTO2::zx] = A.xz() * B.xx() + A.yz() * B.xy() + A.zz() * B.xz();
        res.data[idTO2::zy] = A.xz() * B.xy() + A.yz() * B.yy() + A.zz() * B.yz();
        res.data[idTO2::zz] = A.xz() * B.xz() + A.yz() * B.yz() + A.zz() * B.zz();

        return res;
    }

    // Global Ops SkewSymTensorO2

    [[nodiscard]] constexpr inline SkewSymTensorO2 operator*(CfemReal scalar, const SkewSymTensorO2 &W) noexcept
    {
        return W * scalar;
    }

    /**
     * @brief TensorO2 * SkewSymTensorO2
     * Result is a general TensorO2
     */
    [[nodiscard]] constexpr inline TensorO2 operator*(const TensorO2 &T, const SkewSymTensorO2 &W) noexcept
    {
        TensorO2 res;
        const CfemReal a00 = T.data[idTO2::xx], a01 = T.data[idTO2::xy], a02 = T.data[idTO2::xz];
        const CfemReal a10 = T.data[idTO2::yx], a11 = T.data[idTO2::yy], a12 = T.data[idTO2::yz];
        const CfemReal a20 = T.data[idTO2::zx], a21 = T.data[idTO2::zy], a22 = T.data[idTO2::zz];

        res.data[idTO2::xx] = -a01 * W.xy() - a02 * W.xz();
        res.data[idTO2::yx] = -a11 * W.xy() - a12 * W.xz();
        res.data[idTO2::zx] = -a21 * W.xy() - a22 * W.xz();

        res.data[idTO2::xy] = a00 * W.xy() - a02 * W.yz();
        res.data[idTO2::yy] = a10 * W.xy() - a12 * W.yz();
        res.data[idTO2::zy] = a20 * W.xy() - a22 * W.yz();

        res.data[idTO2::xz] = a00 * W.xz() + a01 * W.yz();
        res.data[idTO2::yz] = a10 * W.xz() + a11 * W.yz();
        res.data[idTO2::zz] = a20 * W.xz() + a21 * W.yz();

        return res;
    }

    /**
     * @brief SkewSymTensorO2 * TensorO2
     * Result is a general TensorO2
     */
    [[nodiscard]] constexpr inline TensorO2 operator*(const SkewSymTensorO2 &W, const TensorO2 &T) noexcept
    {
        TensorO2 res;
        const CfemReal a00 = T.data[idTO2::xx], a01 = T.data[idTO2::xy], a02 = T.data[idTO2::xz];
        const CfemReal a10 = T.data[idTO2::yx], a11 = T.data[idTO2::yy], a12 = T.data[idTO2::yz];
        const CfemReal a20 = T.data[idTO2::zx], a21 = T.data[idTO2::zy], a22 = T.data[idTO2::zz];

        res.data[idTO2::xx] = W.xy() * a10 + W.xz() * a20;
        res.data[idTO2::xy] = W.xy() * a11 + W.xz() * a21;
        res.data[idTO2::xz] = W.xy() * a12 + W.xz() * a22;

        res.data[idTO2::yx] = -W.xy() * a00 + W.yz() * a20;
        res.data[idTO2::yy] = -W.xy() * a01 + W.yz() * a21;
        res.data[idTO2::yz] = -W.xy() * a02 + W.yz() * a22;

        res.data[idTO2::zx] = -W.xz() * a00 - W.yz() * a10;
        res.data[idTO2::zy] = -W.xz() * a01 - W.yz() * a11;
        res.data[idTO2::zz] = -W.xz() * a02 - W.yz() * a12;

        return res;
    }

    /**
     * @brief SkewSymTensorO2 * SkewSymTensorO2
     * Result is a general TensorO2 (symmetric part + trace part)
     */
    [[nodiscard]] constexpr inline TensorO2 operator*(const SkewSymTensorO2 &A, const SkewSymTensorO2 &B) noexcept
    {
        TensorO2 res;

        res.data[idTO2::xx] = -A.xy() * B.xy() - A.xz() * B.xz();
        res.data[idTO2::xy] = -A.xz() * B.yz();
        res.data[idTO2::xz] =  A.xy() * B.yz();

        res.data[idTO2::yx] = -A.yz() * B.xz();
        res.data[idTO2::yy] = -A.xy() * B.xy() - A.yz() * B.yz();
        res.data[idTO2::yz] = -A.xy() * B.xz();

        res.data[idTO2::zx] =  A.yz() * B.xy();
        res.data[idTO2::zy] = -A.xz() * B.xy();
        res.data[idTO2::zz] = -A.xz() * B.xz() - A.yz() * B.yz();

        return res;
    }

    /**
     * @brief SymTensorO2 * SkewSymTensorO2
     * Result is a general TensorO2
     */
    [[nodiscard]] constexpr inline TensorO2 operator*(const SymTensorO2 &S, const SkewSymTensorO2 &W) noexcept
    {
        TensorO2 res;

        res.data[idTO2::xx] =-S.xy() * W.xy() - S.xz() * W.xz();
        res.data[idTO2::xy] = S.xx() * W.xy() - S.xz() * W.yz();
        res.data[idTO2::xz] = S.xx() * W.xz() + S.xy() * W.yz();

        res.data[idTO2::yx] =-S.yy() * W.xy() - S.yz() * W.xz();
        res.data[idTO2::yy] = S.xy() * W.xy() - S.yz() * W.yz();
        res.data[idTO2::yz] = S.xy() * W.xz() + S.yy() * W.yz();

        res.data[idTO2::zx] =-S.yz() * W.xy() - S.zz() * W.xz();
        res.data[idTO2::zy] = S.xz() * W.xy() - S.zz() * W.yz();
        res.data[idTO2::zz] = S.xz() * W.xz() + S.yz() * W.yz();

        return res;
    }

    /**
     * @brief SkewSymTensorO2 * SymTensorO2
     * Result is a general TensorO2
     */
    [[nodiscard]] constexpr inline TensorO2 operator*(const SkewSymTensorO2 &W, const SymTensorO2 &S) noexcept
    {
        TensorO2 res;

        res.data[idTO2::xx] = W.xy() * S.xy() + W.xz() * S.xz();
        res.data[idTO2::xy] = W.xy() * S.yy() + W.xz() * S.yz();
        res.data[idTO2::xz] = W.xy() * S.yz() + W.xz() * S.zz();

        res.data[idTO2::yx] = -W.xy() * S.xx() + W.yz() * S.xz();
        res.data[idTO2::yy] = -W.xy() * S.xy() + W.yz() * S.yz();
        res.data[idTO2::yz] = -W.xy() * S.xz() + W.yz() * S.zz();

        res.data[idTO2::zx] = -W.xz() * S.xx() - W.yz() * S.xy();
        res.data[idTO2::zy] = -W.xz() * S.xy() - W.yz() * S.yy();
        res.data[idTO2::zz] = -W.xz() * S.xz() - W.yz() * S.yz();

        return res;
    }

    // ========================================================================
    // Operators : TensorO2 vs SymTensorO2
    // ========================================================================

    [[nodiscard]] constexpr inline TensorO2 operator+(const TensorO2 &T, const SymTensorO2 &S) noexcept
    {
        TensorO2 res = T;
        res.data[idTO2::xx] += S.xx(); res.data[idTO2::xy] += S.xy(); res.data[idTO2::xz] += S.xz();
        res.data[idTO2::yx] += S.xy(); res.data[idTO2::yy] += S.yy(); res.data[idTO2::yz] += S.yz();
        res.data[idTO2::zx] += S.xz(); res.data[idTO2::zy] += S.yz(); res.data[idTO2::zz] += S.zz();
        return res;
    }

    [[nodiscard]] constexpr inline TensorO2 operator+(const SymTensorO2 &S, const TensorO2 &T) noexcept
    {
        return T + S;
    }

    [[nodiscard]] constexpr inline TensorO2 operator-(const TensorO2 &T, const SymTensorO2 &S) noexcept
    {
        TensorO2 res = T;
        res.data[idTO2::xx] -= S.xx(); res.data[idTO2::xy] -= S.xy(); res.data[idTO2::xz] -= S.xz();
        res.data[idTO2::yx] -= S.xy(); res.data[idTO2::yy] -= S.yy(); res.data[idTO2::yz] -= S.yz();
        res.data[idTO2::zx] -= S.xz(); res.data[idTO2::zy] -= S.yz(); res.data[idTO2::zz] -= S.zz();
        return res;
    }

    [[nodiscard]] constexpr inline TensorO2 operator-(const SymTensorO2 &S, const TensorO2 &T) noexcept
    {
        TensorO2 res;
        res.data[idTO2::xx] = S.xx() - T.data[idTO2::xx]; res.data[idTO2::xy] = S.xy() - T.data[idTO2::xy]; res.data[idTO2::xz] = S.xz() - T.data[idTO2::xz];
        res.data[idTO2::yx] = S.xy() - T.data[idTO2::yx]; res.data[idTO2::yy] = S.yy() - T.data[idTO2::yy]; res.data[idTO2::yz] = S.yz() - T.data[idTO2::yz];
        res.data[idTO2::zx] = S.xz() - T.data[idTO2::zx]; res.data[idTO2::zy] = S.yz() - T.data[idTO2::zy]; res.data[idTO2::zz] = S.zz() - T.data[idTO2::zz];
        return res;
    }

    // ========================================================================
    // Operators : TensorO2 vs SkewSymTensorO2
    // ========================================================================

    [[nodiscard]] constexpr inline TensorO2 operator+(const TensorO2 &T, const SkewSymTensorO2 &W) noexcept
    {
        TensorO2 res = T;
        res.data[idTO2::xy] += W.xy();  res.data[idTO2::xz] += W.xz();
        res.data[idTO2::yx] -= W.xy();  res.data[idTO2::yz] += W.yz();
        res.data[idTO2::zx] -= W.xz();  res.data[idTO2::zy] -= W.yz();
        return res;
    }

    [[nodiscard]] constexpr inline TensorO2 operator+(const SkewSymTensorO2 &W, const TensorO2 &T) noexcept
    {
        return T + W;
    }

    [[nodiscard]] constexpr inline TensorO2 operator-(const TensorO2 &T, const SkewSymTensorO2 &W) noexcept
    {
        TensorO2 res = T;
        res.data[idTO2::xy] -= W.xy();  res.data[idTO2::xz] -= W.xz();
        res.data[idTO2::yx] += W.xy();  res.data[idTO2::yz] -= W.yz();
        res.data[idTO2::zx] += W.xz();  res.data[idTO2::zy] += W.yz();
        return res;
    }

    [[nodiscard]] constexpr inline TensorO2 operator-(const SkewSymTensorO2 &W, const TensorO2 &T) noexcept
    {
        TensorO2 res;
        res.data[idTO2::xx] =         - T.data[idTO2::xx]; res.data[idTO2::xy] =  W.xy() - T.data[idTO2::xy]; res.data[idTO2::xz] =  W.xz() - T.data[idTO2::xz];
        res.data[idTO2::yx] = -W.xy() - T.data[idTO2::yx]; res.data[idTO2::yy] =         - T.data[idTO2::yy]; res.data[idTO2::yz] =  W.yz() - T.data[idTO2::yz];
        res.data[idTO2::zx] = -W.xz() - T.data[idTO2::zx]; res.data[idTO2::zy] = -W.yz() - T.data[idTO2::zy]; res.data[idTO2::zz] =         - T.data[idTO2::zz];
        return res;
    }

    [[nodiscard]] constexpr inline TensorO2 operator+(const SymTensorO2 &S, const SkewSymTensorO2 &W) noexcept
    {
        TensorO2 res;
        res.data[idTO2::xx] = S.xx();    res.data[idTO2::xy] = S.xy() + W.xy(); res.data[idTO2::xz] = S.xz() + W.xz();
        res.data[idTO2::yx] = S.xy() - W.xy(); res.data[idTO2::yy] = S.yy();          res.data[idTO2::yz] = S.yz() + W.yz();
        res.data[idTO2::zx] = S.xz() - W.xz(); res.data[idTO2::zy] = S.yz() - W.yz(); res.data[idTO2::zz] = S.zz();
        return res;
    }

    [[nodiscard]] constexpr inline TensorO2 operator+(const SkewSymTensorO2 &W, const SymTensorO2 &S) noexcept
    {
        return S + W; 
    }

    [[nodiscard]] constexpr inline TensorO2 operator-(const SymTensorO2 &S, const SkewSymTensorO2 &W) noexcept
    {
        TensorO2 res;
        res.data[idTO2::xx] = S.xx();          res.data[idTO2::xy] = S.xy() - W.xy(); res.data[idTO2::xz] = S.xz() - W.xz();
        res.data[idTO2::yx] = S.xy() + W.xy(); res.data[idTO2::yy] = S.yy();          res.data[idTO2::yz] = S.yz() - W.yz();
        res.data[idTO2::zx] = S.xz() + W.xz(); res.data[idTO2::zy] = S.yz() + W.yz(); res.data[idTO2::zz] = S.zz();
        return res;
    }

    [[nodiscard]] constexpr inline TensorO2 operator-(const SkewSymTensorO2 &W, const SymTensorO2 &S) noexcept
    {
        TensorO2 res;
        res.data[idTO2::xx] = -S.xx();          res.data[idTO2::xy] =  W.xy() - S.xy(); res.data[idTO2::xz] =  W.xz() - S.xz();
        res.data[idTO2::yx] = -W.xy() - S.xy(); res.data[idTO2::yy] = -S.yy();          res.data[idTO2::yz] =  W.yz() - S.yz();
        res.data[idTO2::zx] = -W.xz() - S.xz(); res.data[idTO2::zy] = -W.yz() - S.yz(); res.data[idTO2::zz] = -S.zz();
        return res;
    }

    [[nodiscard]] constexpr inline CfemReal det(const SymTensorO2 &S) noexcept
    {
        return S.data[idSTO2::xx] * (S.data[idSTO2::yy] * S.data[idSTO2::zz] - S.data[idSTO2::yz] * S.data[idSTO2::yz])
             - S.data[idSTO2::xy] * (S.data[idSTO2::xy] * S.data[idSTO2::zz] - S.data[idSTO2::yz] * S.data[idSTO2::xz])
             + S.data[idSTO2::xz] * (S.data[idSTO2::xy] * S.data[idSTO2::yz] - S.data[idSTO2::yy] * S.data[idSTO2::xz]);
    }

    [[nodiscard]] constexpr inline CfemReal det(const SkewSymTensorO2 &S) noexcept
    {
        return 0.0;
    }

    [[nodiscard]] constexpr inline SymTensorO2 inverse(const SymTensorO2 &S, CfemReal determinant) noexcept
    {
        if (std::abs(determinant) < cfem::ZERO_TOLERANCE)
            CFEM_WARN << "Singular 3x3 matrix/Jacobian detected.";

        // En pratique, tu peux ajouter ici un assert(std::abs(det) > 1e-14);
        const CfemReal inv_det = 1.0 / determinant;
        
        SymTensorO2 res;
        res.data[idSTO2::xx] = (S.data[idSTO2::yy] * S.data[idSTO2::zz] - S.data[idSTO2::yz] * S.data[idSTO2::yz]) * inv_det;
        res.data[idSTO2::yy] = (S.data[idSTO2::xx] * S.data[idSTO2::zz] - S.data[idSTO2::xz] * S.data[idSTO2::xz]) * inv_det;
        res.data[idSTO2::zz] = (S.data[idSTO2::xx] * S.data[idSTO2::yy] - S.data[idSTO2::xy] * S.data[idSTO2::xy]) * inv_det;
        res.data[idSTO2::xy] = (S.data[idSTO2::xz] * S.data[idSTO2::yz] - S.data[idSTO2::xy] * S.data[idSTO2::zz]) * inv_det;
        res.data[idSTO2::xz] = (S.data[idSTO2::xy] * S.data[idSTO2::yz] - S.data[idSTO2::xx] * S.data[idSTO2::xz]) * inv_det;
        res.data[idSTO2::yz] = (S.data[idSTO2::xy] * S.data[idSTO2::xz] - S.data[idSTO2::yy] * S.data[idSTO2::xz]) * inv_det;
        return res;
    }

    [[nodiscard]] constexpr inline SymTensorO2 inverse(const SymTensorO2 &S) noexcept
    {
        return inverse(S, det(S) );
    }

    /**
     * @brief Tries to inverse a skew-symmetric second order tensor (This operation is impossible).
     * @throw std::runtime_error.
     */
    [[nodiscard]] inline SkewSymTensorO2 inverse(const SkewSymTensorO2 &S) 
    {
        CFEM_THROW("3x3 skew-symmetric tensor is not invertible since det = 0.");
        return SkewSymTensorO2{}; 
    }

    /**
     * @brief Computes the eigenvalues of a symmetric 3x3 tensor analytically.
     *
     * @param S Symmetric tensor.
     * @return Eigenvalues sorted in descending order.
     *
     * @details
     * Uses the closed-form Cardano method for symmetric tensors.
     */
    [[nodiscard]] inline FixedVector<CfemReal, 3> eigenvalues(const SymTensorO2& S) noexcept
    {
        const CfemReal I1 = S.data[idSTO2::xx] + S.data[idSTO2::yy] + S.data[idSTO2::zz];
        const CfemReal I2 = S.data[idSTO2::xx]*S.data[idSTO2::yy] + S.data[idSTO2::yy]*S.data[idSTO2::zz] + S.data[idSTO2::xx]*S.data[idSTO2::zz] 
                          - S.data[idSTO2::xy]*S.data[idSTO2::xy] - S.data[idSTO2::yz]*S.data[idSTO2::yz] - S.data[idSTO2::xz]*S.data[idSTO2::xz];
        const CfemReal I3 = det(S);

        const CfemReal p = I2 - I1 * I1 / 3.0;
        const CfemReal q = I3 - I1 * I2 / 3.0 + 2.0 * I1 * I1 * I1 / 27.0;

        FixedVector<CfemReal, 3> eig;
            
        if (std::abs(p) < ZERO_TOLERANCE) {
            eig[0] = I1 / 3.0; eig[1] = I1 / 3.0; eig[2] = I1 / 3.0;
        } else {
            const CfemReal R = std::sqrt(-p / 3.0);
            CfemReal arg = 3.0 * q / (2.0 * p * R);

            arg = std::clamp(arg, -1.0, 1.0);
                
            const CfemReal theta = std::acos(arg) / 3.0;
            const CfemReal offset = I1 / 3.0;

            eig[0] = 2.0 * R * std::cos(theta) + offset;
            eig[1] = 2.0 * R * std::cos(theta - CFEM_PI_2_3) + offset;
            eig[2] = 2.0 * R * std::cos(theta + CFEM_PI_2_3) + offset;
                
            // Sort to have eig[0] >= eig[1] >= eig[2]
            if (eig[0] < eig[1]) std::swap(eig[0], eig[1]);
            if (eig[1] < eig[2]) std::swap(eig[1], eig[2]);
            if (eig[0] < eig[1]) std::swap(eig[0], eig[1]);
        }
        return eig;
    }

    inline std::ostream &operator<<(std::ostream &os, const SkewSymTensorO2 &W)
    {
        os << "SkewSymTensor2[yz:" << W[0] << ", xz:" << W[1] << ", xy:" << W[2] << "]";
        return os;
    }

        inline std::ostream &operator<<(std::ostream &os, const SymTensorO2 &T)
    {
        os << "SymTensor2[xx:" << T[0] << ", yy:" << T[1] << ", zz:" << T[2]
           << ", yz:" << T[3] << ", xz:" << T[4] << ", xy:" << T[5] << "]";
        return os;
    }


}

#endif