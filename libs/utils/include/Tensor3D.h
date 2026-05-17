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
 * @file Tensor3D.h
 * @brief Physical tensor definitions for CFEM.
 */

#pragma once 

#include "FixedMatrix.h"

namespace cfem
{

    using Mat3x3 = std::array<std::array<CfemReal, 3>, 3>;

    /**
     * @brief General second-order tensor (9 components).
     * Used for non-symmetric quantities like the Deformation Gradient (F)
     * or the Velocity Gradient (L).
     * @note
     * Inherits from FixedMatrix to add specific physical indices.
     */
    struct Tensor3D : public FixedMatrix<CfemReal, 3, 3>
    {
        using FixedMatrix<CfemReal, 3, 3>::FixedMatrix;

        constexpr Tensor3D(const FixedMatrix<CfemReal, 3, 3>& base) noexcept 
            : FixedMatrix<CfemReal, 3, 3>(base) {}

        /**
         * @brief Storage indices for general 3x3 tensors.
         * 
         * @{
         */
        static constexpr CfemInt idx_xx = 0;
        static constexpr CfemInt idx_xy = 1;
        static constexpr CfemInt idx_xz = 2;
        static constexpr CfemInt idx_yx = 3;
        static constexpr CfemInt idx_yy = 4;
        static constexpr CfemInt idx_yz = 5;
        static constexpr CfemInt idx_zx = 6;
        static constexpr CfemInt idx_zy = 7;
        static constexpr CfemInt idx_zz = 8;
        /** @} */

        [[nodiscard]] constexpr CfemReal xx() const noexcept { return data[idx_xx]; }
        [[nodiscard]] constexpr CfemReal xz() const noexcept { return data[idx_xz]; }
        [[nodiscard]] constexpr CfemReal xy() const noexcept { return data[idx_xy]; }
        [[nodiscard]] constexpr CfemReal yx() const noexcept { return data[idx_yx]; }
        [[nodiscard]] constexpr CfemReal yy() const noexcept { return data[idx_yy]; }
        [[nodiscard]] constexpr CfemReal yz() const noexcept { return data[idx_yz]; }
        [[nodiscard]] constexpr CfemReal zx() const noexcept { return data[idx_zx]; }
        [[nodiscard]] constexpr CfemReal zy() const noexcept { return data[idx_zy]; }
        [[nodiscard]] constexpr CfemReal zz() const noexcept { return data[idx_zz]; }

        constexpr CfemReal &xx() noexcept { return data[idx_xx]; }
        constexpr CfemReal &xz() noexcept { return data[idx_xz]; }
        constexpr CfemReal &xy() noexcept { return data[idx_xy]; }
        constexpr CfemReal &yx() noexcept { return data[idx_yx]; }
        constexpr CfemReal &yy() noexcept { return data[idx_yy]; }
        constexpr CfemReal &yz() noexcept { return data[idx_yz]; }
        constexpr CfemReal &zx() noexcept { return data[idx_zx]; }
        constexpr CfemReal &zy() noexcept { return data[idx_zy]; }
        constexpr CfemReal &zz() noexcept { return data[idx_zz]; }

        [[nodiscard]] constexpr Mat3x3 toMat3x3() const noexcept {
            return {{{data[idx_xx], data[idx_xy], data[idx_xz]},
                     {data[idx_yx], data[idx_yy], data[idx_yz]},
                     {data[idx_zx], data[idx_zy], data[idx_zz]}}};
        }
    };

    /**
     * @brief Identity tensor constant.
     */
    [[nodiscard]] constexpr Tensor3D IdentityTensor3D() noexcept
    {
        Tensor3D I;
        I.setIdentity();
        return I;
    }

    /**
     * @brief Stores a symmetric 3x3 tensor using only 6 values.
     * Storage order (Voigt): [0]=xx, [1]=yy, [2]=zz, [3]=yz, [4]=xz, [5]=xy
     */
    struct SymTensor3D
    {
        CfemReal data[6]{};

        /**
         * @brief Storage indices for symmetric 3x3 tensors.
         * 
         * @{
         */
        static constexpr CfemInt idx_xx = 0;
        static constexpr CfemInt idx_yy = 1;
        static constexpr CfemInt idx_zz = 2;
        static constexpr CfemInt idx_yz = 3;
        static constexpr CfemInt idx_xz = 4;
        static constexpr CfemInt idx_xy = 5;

        /** @} */

    private:
        static constexpr CfemInt map[3][3] = {{0, 5, 4}, {5, 1, 3}, {4, 3, 2}};

    public:
        constexpr SymTensor3D() noexcept = default;
        constexpr SymTensor3D(CfemReal xx, CfemReal yy, CfemReal zz,
                              CfemReal yz, CfemReal xz, CfemReal xy) noexcept 
        : data{xx, yy, zz, yz, xz, xy}
        {}

        /** @brief Access by component index (0 to 5). */
        constexpr CfemReal &operator[](CfemInt i)
        {
            CFEM_ASSERT(i >= 0 && i < 6);
            return data[i];
        }

        [[nodiscard]] constexpr const CfemReal &operator[](CfemInt i) const
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
        [[nodiscard]] constexpr CfemReal operator()(CfemInt i, CfemInt j) const
        {
            CFEM_ASSERT(i >= 0 && i < 3 && j >= 0 && j < 3);
            return data[map[i][j]];
        }

        [[nodiscard]] constexpr Mat3x3 toMat3x3() const noexcept{
            return {{{data[idx_xx], data[idx_xy], data[idx_xz]},
                     {data[idx_xy], data[idx_yy], data[idx_yz]},
                     {data[idx_xz], data[idx_yz], data[idx_zz]}}};
        }

        /**j;
         * @brief Negate operator
         *
         * @return SymTensor3D
         */
        [[nodiscard]] constexpr SymTensor3D operator-() const noexcept
        {
            return {-data[idx_xx], -data[idx_yy], -data[idx_zz],
                    -data[idx_yz], -data[idx_xz], -data[idx_xy]};
        }

        [[nodiscard]] constexpr CfemReal xx() const noexcept { return data[idx_xx]; }
        [[nodiscard]] constexpr CfemReal yy() const noexcept { return data[idx_yy]; }
        [[nodiscard]] constexpr CfemReal zz() const noexcept { return data[idx_zz]; }
        [[nodiscard]] constexpr CfemReal yz() const noexcept { return data[idx_yz]; }
        [[nodiscard]] constexpr CfemReal xz() const noexcept { return data[idx_xz]; }
        [[nodiscard]] constexpr CfemReal xy() const noexcept { return data[idx_xy]; }

        constexpr CfemReal &xx() noexcept { return data[idx_xx]; }
        constexpr CfemReal &yy() noexcept { return data[idx_yy]; }
        constexpr CfemReal &zz() noexcept { return data[idx_zz]; }
        constexpr CfemReal &yz() noexcept { return data[idx_yz]; }
        constexpr CfemReal &xz() noexcept { return data[idx_xz]; }
        constexpr CfemReal &xy() noexcept { return data[idx_xy]; }

        inline operator Tensor3D() const noexcept{
            return Tensor3D{ data[SymTensor3D::idx_xx], data[SymTensor3D::idx_xy], data[SymTensor3D::idx_xz],
                             data[SymTensor3D::idx_xy], data[SymTensor3D::idx_yy], data[SymTensor3D::idx_yz],
                             data[SymTensor3D::idx_xz], data[SymTensor3D::idx_yz], data[SymTensor3D::idx_zz] };
        }
        // Basic Math Ops

        [[nodiscard]] static constexpr SymTensor3D Identity() noexcept
        {
            return {1, 1, 1, 0, 0, 0};
        }

        [[nodiscard]] static constexpr SymTensor3D Zero() noexcept
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

        [[nodiscard]] constexpr SymTensor3D operator+(const SymTensor3D &other) const noexcept
        {
            return {data[idx_xx] + other.data[idx_xx], data[idx_yy] + other.data[idx_yy], data[idx_zz] + other.data[idx_zz],
                    data[idx_yz] + other.data[idx_yz], data[idx_xz] + other.data[idx_xz], data[idx_xy] + other.data[idx_xy]};
        }

        [[nodiscard]] constexpr SymTensor3D operator-(const SymTensor3D &other) const noexcept
        {
            return {data[idx_xx] - other.data[idx_xx], data[idx_yy] - other.data[idx_yy], data[idx_zz] - other.data[idx_zz],
                    data[idx_yz] - other.data[idx_yz], data[idx_xz] - other.data[idx_xz], data[idx_xy] - other.data[idx_xy]};
        }

        [[nodiscard]] constexpr SymTensor3D operator*(CfemReal scalar) const noexcept
        {
            return {data[idx_xx] * scalar, data[idx_yy] * scalar, data[idx_zz] * scalar,
                    data[idx_yz] * scalar, data[idx_xz] * scalar, data[idx_xy] * scalar};
        }

        [[nodiscard]] constexpr SymTensor3D operator/(CfemReal scalar) const
        {
            CFEM_ASSERT(std::abs(scalar) > ZERO_TOLERANCE, "Division by zero in SymTensor3D::operator/");

            CfemReal inv = 1.0 / scalar;
            return {data[idx_xx] * inv, data[idx_yy] * inv, data[idx_zz] * inv,
                    data[idx_yz] * inv, data[idx_xz] * inv, data[idx_xy] * inv};
        }

        [[nodiscard]] constexpr Vector3D operator*(const Vector3D &v) const noexcept
        {
            return {
                data[idx_xx] * v.x + data[idx_xy] * v.y + data[idx_xz] * v.z, // Row x: [xx, xy, xz]
                data[idx_xy] * v.x + data[idx_yy] * v.y + data[idx_yz] * v.z, // Row y: [yx, yy, yz]
                data[idx_xz] * v.x + data[idx_yz] * v.y + data[idx_zz] * v.z  // Row z: [zx, zy, zz]
            };
        }

        [[nodiscard]] constexpr Tensor3D operator*(const SymTensor3D &B) const noexcept
        {
            Tensor3D res;

            const CfemReal Axx = data[idx_xx], Ayy = data[idx_yy], Azz = data[idx_zz], Axy = data[idx_xy], Axz = data[idx_xz], Ayz = data[idx_yz];
            const CfemReal Bxx = B.data[SymTensor3D::idx_xx], Byy = B.data[SymTensor3D::idx_yy], Bzz = B.data[SymTensor3D::idx_zz];
            const CfemReal Bxy = B.data[SymTensor3D::idx_xy], Bxz = B.data[SymTensor3D::idx_xz], Byz = B.data[SymTensor3D::idx_yz];

            res.data[Tensor3D::idx_xx] = Axx * Bxx + Axy * Bxy + Axz * Bxz;
            res.data[Tensor3D::idx_xy] = Axx * Bxy + Axy * Byy + Axz * Byz;
            res.data[Tensor3D::idx_xz] = Axx * Bxz + Axy * Byz + Axz * Bzz;

            res.data[Tensor3D::idx_yx] = Axy * Bxx + Ayy * Bxy + Ayz * Bxz;
            res.data[Tensor3D::idx_yy] = Axy * Bxy + Ayy * Byy + Ayz * Byz;
            res.data[Tensor3D::idx_yz] = Axy * Bxz + Ayy * Byz + Ayz * Bzz;

            res.data[Tensor3D::idx_zx] = Axz * Bxx + Ayz * Bxy + Azz * Bxz;
            res.data[Tensor3D::idx_zy] = Axz * Bxy + Ayz * Byy + Azz * Byz;
            res.data[Tensor3D::idx_zz] = Axz * Bxz + Ayz * Byz + Azz * Bzz;

            return res;
        }

        constexpr SymTensor3D &operator+=(const SymTensor3D &other) noexcept
        {
            data[idx_xx] += other.data[idx_xx];
            data[idx_yy] += other.data[idx_yy];
            data[idx_zz] += other.data[idx_zz];
            data[idx_yz] += other.data[idx_yz];
            data[idx_xz] += other.data[idx_xz];
            data[idx_xy] += other.data[idx_xy];
            return *this;
        }

        constexpr SymTensor3D &operator-=(const SymTensor3D &other) noexcept
        {
            data[idx_xx] -= other.data[idx_xx];
            data[idx_yy] -= other.data[idx_yy];
            data[idx_zz] -= other.data[idx_zz];
            data[idx_yz] -= other.data[idx_yz];
            data[idx_xz] -= other.data[idx_xz];
            data[idx_xy] -= other.data[idx_xy];
            return *this;
        }

        constexpr SymTensor3D &operator*=(CfemReal scalar) noexcept
        {
            data[idx_xx] *= scalar;
            data[idx_yy] *= scalar;
            data[idx_zz] *= scalar;
            data[idx_yz] *= scalar;
            data[idx_xz] *= scalar;
            data[idx_xy] *= scalar;
            return *this;
        }

        constexpr SymTensor3D &operator/=(CfemReal scalar)
        {
            CFEM_ASSERT(std::abs(scalar) > ZERO_TOLERANCE, "Division by zero in SymTensor3D::operator/=");

            CfemReal inv = 1.0 / scalar;

            data[idx_yy] *= inv;
            data[idx_zz] *= inv;
            data[idx_yz] *= inv;
            data[idx_xz] *= inv;
            data[idx_xy] *= inv;
            data[idx_xx] *= inv;
            return *this;
        }

        [[nodiscard]] constexpr bool operator==(const SymTensor3D &S) const noexcept
        {
            for (CfemInt i = 0; i < 6; ++i)
                if (data[i] != S.data[i])
                    return false;
            return true;
        }

        [[nodiscard]] constexpr bool operator!=(const SymTensor3D &S) const noexcept
        {
            return !(*this == S);
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
        [[nodiscard]] constexpr SymTensor3D quadform(const Tensor3D &M) const noexcept
        {
            SymTensor3D res;

            // On calcule d'abord le produit intermédiaire T = (*this) * M
            // Rappel stockage : [0]=xx, [1]=yy, [2]=zz, [3]=yz, [4]=xz, [5]=xy
            Tensor3D T;

            const CfemReal m00 = M.data[Tensor3D::idx_xx], m01 = M.data[Tensor3D::idx_xy], m02 = M.data[Tensor3D::idx_xz];
            const CfemReal m10 = M.data[Tensor3D::idx_yx], m11 = M.data[Tensor3D::idx_yy], m12 = M.data[Tensor3D::idx_yz];
            const CfemReal m20 = M.data[Tensor3D::idx_zx], m21 = M.data[Tensor3D::idx_zy], m22 = M.data[Tensor3D::idx_zz];

            // Ligne 0 de (*this) : [xx, xy, xz]
            T.data[Tensor3D::idx_xx] = data[idx_xx] * m00 + data[idx_xy] * m10 + data[idx_xz] * m20;
            T.data[Tensor3D::idx_xy] = data[idx_xx] * m01 + data[idx_xy] * m11 + data[idx_xz] * m21;
            T.data[Tensor3D::idx_xz] = data[idx_xx] * m02 + data[idx_xy] * m12 + data[idx_xz] * m22;

            // Ligne 1 de (*this) : [xy, yy, yz]
            T.data[Tensor3D::idx_yx] = data[idx_xy] * m00 + data[idx_yy] * m10 + data[idx_yz] * m20;
            T.data[Tensor3D::idx_yy] = data[idx_xy] * m01 + data[idx_yy] * m11 + data[idx_yz] * m21;
            T.data[Tensor3D::idx_yz] = data[idx_xy] * m02 + data[idx_yy] * m12 + data[idx_yz] * m22;

            // Ligne 2 de (*this) : [xz, yz, zz]
            T.data[Tensor3D::idx_zx] = data[idx_xz] * m00 + data[idx_yz] * m10 + data[idx_zz] * m20;
            T.data[Tensor3D::idx_zy] = data[idx_xz] * m01 + data[idx_yz] * m11 + data[idx_zz] * m21;
            T.data[Tensor3D::idx_zz] = data[idx_xz] * m02 + data[idx_yz] * m12 + data[idx_zz] * m22;

            // On calcule res = M^T * T
            // On ne calcule que les 6 composantes symétriques
            res[0] = m00 * T.data[Tensor3D::idx_xx] + m10 * T.data[Tensor3D::idx_yx] + m20 * T.data[Tensor3D::idx_zx]; // xx
            res[1] = m01 * T.data[Tensor3D::idx_xy] + m11 * T.data[Tensor3D::idx_yy] + m21 * T.data[Tensor3D::idx_zy]; // yy
            res[2] = m02 * T.data[Tensor3D::idx_xz] + m12 * T.data[Tensor3D::idx_yz] + m22 * T.data[Tensor3D::idx_zz]; // zz

            res[3] = m01 * T.data[Tensor3D::idx_xz] + m11 * T.data[Tensor3D::idx_yz] + m21 * T.data[Tensor3D::idx_zz]; // yz
            res[4] = m00 * T.data[Tensor3D::idx_xz] + m10 * T.data[Tensor3D::idx_yz] + m20 * T.data[Tensor3D::idx_zz]; // xz
            res[5] = m00 * T.data[Tensor3D::idx_xy] + m10 * T.data[Tensor3D::idx_yy] + m20 * T.data[Tensor3D::idx_zy]; // xy

            return res;
        }

        // Physical Invariants

        /** @brief Trace of the tensor (sum of diagonals). */
        [[nodiscard]] constexpr CfemReal trace() const noexcept
        {
            return data[idx_xx] + data[idx_yy] + data[idx_zz];
        }

        /** @brief Mean (hydrostatic) component: p = tr(sigma)/3 */
        [[nodiscard]] constexpr CfemReal mean() const noexcept
        {
            return trace() / 3.0;
        }

        /** @brief Returns the deviatoric part: s = sigma - p*I */
        [[nodiscard]] constexpr SymTensor3D deviator() const noexcept
        {
            CfemReal p = mean();
            return {data[idx_xx] - p, data[idx_yy] - p, data[idx_zz] - p, data[idx_yz], data[idx_xz], data[idx_xy]};
        }

        /** @brief Double contraction with another symmetric tensor (A : B). */
        [[nodiscard]] constexpr CfemReal doubleDot(const SymTensor3D &other) const noexcept
        {
            // xx*xx + yy*yy + zz*zz + 2*(yz*yz + xz*xz + xy*xy)
            return data[idx_xx] * other.data[idx_xx] + data[idx_yy] * other.data[idx_yy] + data[idx_zz] * other.data[idx_zz] +
                   2.0 * (data[idx_yz] * other.data[idx_yz] + data[idx_xz] * other.data[idx_xz] + data[idx_xy] * other.data[idx_xy]);
        }

        /** @brief Second invariant of deviatoric tensor: J2 = 0.5 * (s : s) */
        [[nodiscard]] constexpr CfemReal J2() const noexcept
        {
            SymTensor3D s = deviator();
            return 0.5 * s.doubleDot(s);
        }

        /** @brief Computes Von Mises equivalent stress (if this is a stress tensor). */
        [[nodiscard]] inline CfemReal vonMises() const noexcept
        {
            CfemReal t = trace() / 3.0;
            CfemReal s_xx = data[idx_xx] - t;
            CfemReal s_yy = data[idx_yy] - t;
            CfemReal s_zz = data[idx_zz] - t;
            return std::sqrt(1.5 * (s_xx * s_xx + s_yy * s_yy + s_zz * s_zz +
                                    2.0 * (data[idx_yz] * data[idx_yz] + data[idx_xz] * data[idx_xz] + data[idx_xy] * data[idx_xy])));
        }

        /** 
         * @brief Norm of the tensor (Frobenius norm) 
         */
        [[nodiscard]] inline CfemReal norm() const noexcept
        {
            return std::sqrt(this->doubleDot(*this));
        }

        /**
         * @brief Compute the transpose of the matrix
         * 
         * @return constexpr SymTensor3D The transposed matrix
         * 
         * @note the transpose of a symmetric matrix is exactely equal to the matrix itself
         */
        [[nodiscard]] constexpr SymTensor3D transpose() const noexcept
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
        [[nodiscard]] constexpr CfemReal doubleDot(const Tensor3D& T) const noexcept
        {
            // S : T = S_ij T_ij
            return data[idx_xx] * T.data[Tensor3D::idx_xx] + 
                data[idx_yy] * T.data[Tensor3D::idx_yy] + 
                data[idx_zz] * T.data[Tensor3D::idx_zz] + 
                data[idx_xy] *(T.data[Tensor3D::idx_xy] + T.data[Tensor3D::idx_yx]) + 
                data[idx_xz] *(T.data[Tensor3D::idx_xz] + T.data[Tensor3D::idx_zx]) + 
                data[idx_yz] *(T.data[Tensor3D::idx_yz] + T.data[Tensor3D::idx_zy]);
        }

        [[nodiscard]] constexpr SymTensor3D deviatoric() const noexcept
        {
            const CfemReal p = (data[idx_xx] + data[idx_yy] + data[idx_zz]) / 3.0;
            SymTensor3D res = *this;
            res.data[idx_xx] -= p;
            res.data[idx_yy] -= p;
            res.data[idx_zz] -= p;
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
        [[nodiscard]] constexpr FixedMatrix<CfemReal, 6, 6> outer(const SymTensor3D &other) const noexcept
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


    /**
     * @brief Symmetric Axisymmetric Tensor (4 components).
     * Components stored as: [xx, yy, zz, xy] (e.g., rr, zz, theta_theta, rz).
     */
    struct AxiSymTensor3D
    {
        static constexpr CfemInt idx_xx = 0; 
        static constexpr CfemInt idx_yy = 1;
        static constexpr CfemInt idx_zz = 2;
        static constexpr CfemInt idx_xy = 3;

        CfemReal data[4]{0.0};

    private:
        static constexpr CfemInt map2d[2][2] = {{idx_xx, idx_xy}, {idx_xy, idx_yy}};

    public:

        // Constructors
        constexpr AxiSymTensor3D() noexcept = default;
        constexpr AxiSymTensor3D(CfemReal xx, CfemReal yy, CfemReal zz, CfemReal xy) noexcept 
            : data{xx, yy, zz, xy} {}

        // Accessors
        [[nodiscard]] constexpr CfemReal xx() const noexcept { return data[idx_xx]; }
        [[nodiscard]] constexpr CfemReal yy() const noexcept { return data[idx_yy]; }
        [[nodiscard]] constexpr CfemReal zz() const noexcept { return data[idx_zz]; }
        [[nodiscard]] constexpr CfemReal xy() const noexcept { return data[idx_xy]; }
        
        constexpr CfemReal &xx() noexcept { return data[idx_xx]; }
        constexpr CfemReal &yy() noexcept { return data[idx_yy]; }
        constexpr CfemReal &zz() noexcept { return data[idx_zz]; }
        constexpr CfemReal &xy() noexcept { return data[idx_xy]; }

        inline operator Tensor3D() const noexcept{
            return Tensor3D{ data[idx_xx], data[idx_xy], 0.0, 
                             data[idx_xy], data[idx_yy], 0.0,
                             0.0,  0.0,  data[idx_zz] };
        }

        inline operator SymTensor3D() const noexcept{
            return SymTensor3D{ data[idx_xx],  data[idx_yy], data[idx_zz], 0.0, 0.0, data[idx_xy] };
        }


        [[nodiscard]] constexpr Mat3x3 toMat3x3() const noexcept{
            return {{{data[idx_xx], data[idx_xy], 0.0}, 
                     {data[idx_xy], data[idx_yy], 0.0},
                     {0.0,  0.0,  data[idx_zz]}}};
        }

        [[nodiscard]] constexpr AxiSymTensor3D operator-() const noexcept
        {
            return {-data[idx_xx], -data[idx_yy], -data[idx_zz], -data[idx_xy]};
        }

        [[nodiscard]] constexpr AxiSymTensor3D operator*(CfemReal scalar) const noexcept
        {
            return { data[idx_xx]*scalar, data[idx_yy]*scalar, data[idx_zz]*scalar, data[idx_xy]*scalar };
        }

        [[nodiscard]] constexpr AxiSymTensor3D operator/(CfemReal scalar) const
        {
            CFEM_ASSERT(std::abs(scalar) > ZERO_TOLERANCE, "Division by zero in SymTensor3D::operator/");

            const CfemReal inv = 1.0 / scalar;
            return { data[idx_xx]*inv, data[idx_yy]*inv, data[idx_zz]*inv, data[idx_xy]*inv };
        }

        [[nodiscard]] constexpr AxiSymTensor3D operator+(const AxiSymTensor3D &A) const noexcept
        {
            return { data[idx_xx] + A.data[AxiSymTensor3D::idx_xx], data[idx_yy] + A.data[AxiSymTensor3D::idx_yy], data[idx_zz] + A.data[AxiSymTensor3D::idx_zz], data[idx_xy] + A.data[AxiSymTensor3D::idx_xy] };
        }

        [[nodiscard]] constexpr AxiSymTensor3D operator-(const AxiSymTensor3D &A) const noexcept
        {
            return { data[idx_xx] - A.data[AxiSymTensor3D::idx_xx], data[idx_yy] - A.data[AxiSymTensor3D::idx_yy], data[idx_zz] - A.data[AxiSymTensor3D::idx_zz], data[idx_xy] - A.data[AxiSymTensor3D::idx_xy] };
        }

        [[nodiscard]] constexpr Vector3D operator*(const Vector3D &v) const noexcept
        {
            return {
                data[idx_xx] * v.x + data[idx_xy] * v.y, // Row x: [xx, xy, xz]
                data[idx_xy] * v.x + data[idx_yy] * v.y, // Row y: [yx, yy, yz]
                data[idx_zz] * v.z                       // Row z: [zx, zy, zz]
            };
        }

        [[nodiscard]] constexpr Tensor3D operator*(const AxiSymTensor3D &A) const noexcept
        {
            Tensor3D res;
            res.data[Tensor3D::idx_xx] = data[idx_xx]*A.data[AxiSymTensor3D::idx_xx] + data[idx_xy]*A.data[AxiSymTensor3D::idx_xy];
            res.data[Tensor3D::idx_xy] = data[idx_xx]*A.data[AxiSymTensor3D::idx_xy] + data[idx_xy]*A.data[AxiSymTensor3D::idx_yy];
            res.data[Tensor3D::idx_xz] = 0.0;
            
            res.data[Tensor3D::idx_yx] = data[idx_xy]*A.data[AxiSymTensor3D::idx_xx] + data[idx_yy]*A.data[AxiSymTensor3D::idx_xy];
            res.data[Tensor3D::idx_yy] = data[idx_xy]*A.data[AxiSymTensor3D::idx_xy] + data[idx_yy]*A.data[AxiSymTensor3D::idx_yy];
            res.data[Tensor3D::idx_yz] = 0.0;
            
            res.data[Tensor3D::idx_zx] = 0.0;
            res.data[Tensor3D::idx_zy] = 0.0;
            res.data[Tensor3D::idx_zz] = data[idx_zz]*A.data[AxiSymTensor3D::idx_zz];
            return res;
        }
        
        constexpr AxiSymTensor3D& operator+=(const AxiSymTensor3D &A) noexcept
        {
            data[idx_xx] += A.data[idx_xx]; data[idx_yy] += A.data[idx_yy];
            data[idx_zz] += A.data[idx_zz]; data[idx_xy] += A.data[idx_xy];
            return *this;
        }

        constexpr AxiSymTensor3D& operator-=(const AxiSymTensor3D &A) noexcept
        {
            data[idx_xx] -= A.data[idx_xx]; data[idx_yy] -= A.data[idx_yy];
            data[idx_zz] -= A.data[idx_zz]; data[idx_xy] -= A.data[idx_xy];
            return *this;
        }

        constexpr AxiSymTensor3D& operator*=(CfemReal scalar) noexcept
        {
            data[idx_xx] *= scalar; data[idx_yy] *= scalar;
            data[idx_zz] *= scalar; data[idx_xy] *= scalar;
            return *this;
        }

        constexpr AxiSymTensor3D& operator/=(CfemReal scalar)
        {
            CFEM_ASSERT(std::abs(scalar) > ZERO_TOLERANCE, "Division by zero in AxiSymTensor3D::operator/=");

            const CfemReal inv = 1.0 / scalar;
            data[idx_xx] *= inv; data[idx_yy] *= inv;
            data[idx_zz] *= inv; data[idx_xy] *= inv;
            return *this;
        }
        // Math Operations
        
        [[nodiscard]] constexpr CfemReal trace() const noexcept
        {
            return data[idx_xx] + data[idx_yy] + data[idx_zz];
        }

        [[nodiscard]] constexpr CfemReal determinant() const noexcept
        {
            return data[idx_zz] * (data[idx_xx] * data[idx_yy] - data[idx_xy] * data[idx_xy]);
        }

        [[nodiscard]] constexpr AxiSymTensor3D transpose() const noexcept
        {
            return *this; // Symétrique !
        }

        [[nodiscard]] constexpr AxiSymTensor3D deviatoric() const noexcept
        {
            const CfemReal mean_tr = trace() / 3.0;
            return { data[idx_xx] - mean_tr, 
                     data[idx_yy] - mean_tr, 
                     data[idx_zz] - mean_tr, 
                     data[idx_xy] };
        }

        [[nodiscard]] constexpr AxiSymTensor3D inverse() const noexcept
        {
            const CfemReal det_2d = data[idx_xx] * data[idx_yy] - data[idx_xy] * data[idx_xy];
            const CfemReal inv_det_2d = 1.0 / det_2d;
            const CfemReal inv_zz = 1.0 / data[idx_zz];

            return {  data[idx_yy] * inv_det_2d,
                      data[idx_xx] * inv_det_2d,
                      inv_zz,
                     -data[idx_xy] * inv_det_2d };
        }

        // Produit doublement contracté (S : S2) = S_ij S2_ij
        [[nodiscard]] constexpr CfemReal doubleDot(const AxiSymTensor3D &A) const noexcept
        {
            return data[idx_xx] * A.data[idx_xx] + 
                   data[idx_yy] * A.data[idx_yy] + 
                   data[idx_zz] * A.data[idx_zz] + 
                   2.0 * data[idx_xy] * A.data[idx_xy];
        }

        // Produit doublement contracté (S : T) 
        [[nodiscard]] constexpr CfemReal doubleDot(const Tensor3D &T) const noexcept
        {
            return data[idx_xx]*T.data[Tensor3D::idx_xx] + 
                   data[idx_yy]*T.data[Tensor3D::idx_yy] + 
                   data[idx_zz]*T.data[Tensor3D::idx_zz] + 
                   data[idx_xy]*(T.data[Tensor3D::idx_xy] + T.data[Tensor3D::idx_yx]);
        }

        // Produit tensoriel externe S (x) S2 -> Matrice de Voigt 4x4 (souvent K_tangent local)
        [[nodiscard]] constexpr FixedMatrix<CfemReal, 4, 4> outer(const AxiSymTensor3D& other) const noexcept
        {
            FixedMatrix<CfemReal, 4, 4> res;
            for(CfemInt i = 0; i < 4; ++i) {
                for(CfemInt j = 0; j < 4; ++j) {
                    res(i, j) = data[i] * other.data[j];
                }
            }
            return res;
        }

        [[nodiscard]] inline FixedVector<CfemReal, 3> eigenvalues() const noexcept
        {
            const CfemReal tr_2d = data[idx_xx] + data[idx_yy];
            const CfemReal gap = data[idx_xx] - data[idx_yy];
            const CfemReal desc = std::sqrt(gap * gap + 4.0 * data[idx_xy] * data[idx_xy]);

            FixedVector<CfemReal, 3> eig;
            eig[0] = 0.5 * (tr_2d + desc);
            eig[1] = 0.5 * (tr_2d - desc);
            eig[2] = data[idx_zz];

            // Tri décroissant
            if (eig[0] < eig[1]) std::swap(eig[0], eig[1]);
            if (eig[1] < eig[2]) std::swap(eig[1], eig[2]);
            if (eig[0] < eig[1]) std::swap(eig[0], eig[1]);
            
            return eig;
        }

        // Mapping pour operator(i,j) si nécessaire pour compatibilité
        [[nodiscard]] constexpr CfemReal operator()(CfemInt i, CfemInt j) const noexcept
        {
            if (i == 2 && j == 2) return data[idx_zz];
            if (i == 2 || j == 2) return 0.0; // xz = yz = 0
            
            return data[map2d[i][j]];
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
     * Storage order matches the off-diagonal parts of SymTensor3D:
     * [0]=yz, [1]=xz, [2]=xy
     */
    struct SkewSymTensor3D
    {
        CfemReal data[3]{};

        /**
         * @brief Storage indices for skew-symmetric 3x3 tensors.
         * 
         * @{
         */
        static constexpr CfemInt idx_yz = 0;
        static constexpr CfemInt idx_xz = 1;
        static constexpr CfemInt idx_xy = 2;
        /** @} */

    public:
        constexpr SkewSymTensor3D() noexcept = default;
        constexpr SkewSymTensor3D(CfemReal yz, CfemReal xz, CfemReal xy ):
        data{yz, xz, xy}
        { }

        /** @brief Access by internal storage index (0 to 2). */
        constexpr CfemReal &operator[](CfemInt i)
        {
            CFEM_ASSERT(i >= 0 && i < 3);
            return data[i];
        }

        [[nodiscard]] constexpr CfemReal operator[](CfemInt i) const
        {
            CFEM_ASSERT(i >= 0 && i < 3);
            return data[i];
        }

        /**
         * @brief Map 3x3 indices to the skew-symmetric storage.
         * Returns 0 for diagonal terms, and handles signs automatically.
         */
        [[nodiscard]] constexpr CfemReal operator()(CfemInt i, CfemInt j) const
        {
            CFEM_ASSERT(i >= 0 && i < 3 && j >= 0 && j < 3, "Index out of bounds in SkewSymTensor3D::operator()");

            switch (3*i + j)
            {
                case 1: return  data[idx_xy]; // (0,1)
                case 2: return  data[idx_xz]; // (0,2)

                case 3: return -data[idx_xy]; // (1,0)
                case 5: return  data[idx_yz]; // (1,2)

                case 6: return -data[idx_xz]; // (2,0)
                case 7: return -data[idx_yz]; // (2,1)

                default: return 0.0; // diagonal
            }
        }

        [[nodiscard]] constexpr Mat3x3 toMat3x3() const noexcept {
            return {{{0.0,   data[idx_xy], data[idx_xz]}, 
                     {-data[idx_xy], 0.0,  data[idx_yz]},
                     {-data[idx_xz], -data[idx_yz], 0.0}}};
        }
        

        inline operator Tensor3D() const noexcept{
            return Tensor3D{ 0.0,   data[idx_xy], data[idx_xz], 
                            -data[idx_xy],  0.0,  data[idx_yz],
                            -data[idx_xz], -data[idx_yz], 0.0 };
        }

        /** @brief Negate operator */
        [[nodiscard]] constexpr SkewSymTensor3D operator-() const noexcept
        {
            return {-data[idx_yz], -data[idx_xz], -data[idx_xy]};
        }

        [[nodiscard]] constexpr CfemReal yz() const noexcept { return data[idx_yz]; }
        [[nodiscard]] constexpr CfemReal xz() const noexcept { return data[idx_xz]; }
        [[nodiscard]] constexpr CfemReal xy() const noexcept { return data[idx_xy]; }

        constexpr CfemReal &yz() noexcept { return data[idx_yz]; }
        constexpr CfemReal &xz() noexcept { return data[idx_xz]; }
        constexpr CfemReal &xy() noexcept { return data[idx_xy]; }

        // Basic Math Ops

        [[nodiscard]] static constexpr SkewSymTensor3D Zero() noexcept
        {
            return {0, 0, 0};
        }

        constexpr void setZero() noexcept
        {
            data[idx_yz] = 0;
            data[idx_xz] = 0;
            data[idx_xy] = 0;
        }

        constexpr void setTo(CfemReal value) noexcept
        {
            data[idx_yz] = value;
            data[idx_xz] = value;
            data[idx_xy] = value;
        }

        constexpr void fromMatrix(CfemReal M[3][3]) noexcept
        {
            // We assume M is already skew-symmetric, we extract the upper triangle
            data[idx_yz] = M[1][2]; // yz
            data[idx_xz] = M[0][2]; // xz
            data[idx_xy] = M[0][1]; // xy
        }

        [[nodiscard]] constexpr SkewSymTensor3D operator+(const SkewSymTensor3D &other) const noexcept
        {
            return {data[idx_yz] + other.data[idx_yz], data[idx_xz] + other.data[idx_xz], data[idx_xy] + other.data[idx_xy]};
        }

        [[nodiscard]] constexpr SkewSymTensor3D operator-(const SkewSymTensor3D &other) const noexcept
        {
            return {data[idx_yz] - other.data[idx_yz], data[idx_xz] - other.data[idx_xz], data[idx_xy] - other.data[idx_xy]};
        }

        [[nodiscard]] constexpr SkewSymTensor3D operator*(CfemReal scalar) const noexcept
        {
            return {data[idx_yz] * scalar, data[idx_xz] * scalar, data[idx_xy] * scalar};
        }

        [[nodiscard]] constexpr SkewSymTensor3D operator/(CfemReal scalar) const
        {
            CFEM_ASSERT(std::abs(scalar) > ZERO_TOLERANCE, "Division by zero in SkewSymTensor3D::operator/");

            CfemReal inv = 1.0 / scalar;
            return {data[idx_yz] * inv, data[idx_xz] * inv, data[idx_xy] * inv};
        }

        /**
         * @brief Computes W * v (Equivalent to the cross product: axial_vector x v)
         */
        [[nodiscard]] constexpr Vector3D operator*(const Vector3D &v) const noexcept
        {
            return {
                 data[idx_xy] * v.y + data[idx_xz] * v.z, // Row x: [0, xy, xz]
                -data[idx_xy] * v.x + data[idx_yz] * v.z, // Row y: [-xy, 0, yz]
                -data[idx_xz] * v.x - data[idx_yz] * v.y  // Row z: [-xz, -yz, 0]
            };
        }

        [[nodiscard]] constexpr Tensor3D operator*(const SkewSymTensor3D &W) const noexcept
        {
            Tensor3D res;
            const CfemReal Axy = data[idx_xy], Axz = data[idx_xz], Ayz = data[idx_yz];
            const CfemReal Bxy = W.data[SkewSymTensor3D::idx_xy], Bxz = W.data[SkewSymTensor3D::idx_xz], Byz = W.data[SkewSymTensor3D::idx_yz];

            res.data[Tensor3D::idx_xx] = -Axy * Bxy - Axz * Bxz;
            res.data[Tensor3D::idx_xy] = -Axz * Byz;
            res.data[Tensor3D::idx_xz] =  Axy * Byz;

            res.data[Tensor3D::idx_yx] = -Ayz * Bxz;
            res.data[Tensor3D::idx_yy] = -Axy * Bxy - Ayz * Byz;
            res.data[Tensor3D::idx_yz] = -Axy * Bxz;

            res.data[Tensor3D::idx_zx] =  Ayz * Bxy;
            res.data[Tensor3D::idx_zy] = -Axz * Bxy;
            res.data[Tensor3D::idx_zz] = -Axz * Bxz - Ayz * Byz;

            return res;
        }

        constexpr SkewSymTensor3D &operator+=(const SkewSymTensor3D &other) noexcept
        {
            data[idx_yz] += other.data[idx_yz];
            data[idx_xz] += other.data[idx_xz];
            data[idx_xy] += other.data[idx_xy];
            return *this;
        }

        constexpr SkewSymTensor3D &operator-=(const SkewSymTensor3D &other) noexcept
        {
            data[idx_yz] -= other.data[idx_yz];
            data[idx_xz] -= other.data[idx_xz];
            data[idx_xy] -= other.data[idx_xy];
            return *this;
        }

        constexpr SkewSymTensor3D &operator*=(CfemReal scalar) noexcept
        {
            data[idx_yz] *= scalar;
            data[idx_xz] *= scalar;
            data[idx_xy] *= scalar;
            return *this;
        }

        constexpr SkewSymTensor3D &operator/=(CfemReal scalar)
        {
            CFEM_ASSERT(std::abs(scalar) > ZERO_TOLERANCE, "Division by zero in SkewSymTensor3D::operator/=");

            CfemReal inv = 1.0/ scalar;
            data[idx_yz] *= inv;
            data[idx_xz] *= inv;
            data[idx_xy] *= inv;
            return *this;
        }

        [[nodiscard]] constexpr bool operator==(const SkewSymTensor3D &other) const noexcept
        {
            for (CfemInt i = 0; i < 3; ++i)
                if (data[i] != other.data[i])
                    return false;
            return true;
        }

        [[nodiscard]] constexpr bool operator!=(const SkewSymTensor3D &other) const noexcept
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
        [[nodiscard]] constexpr SkewSymTensor3D quadform(const Tensor3D &M) const noexcept
        {
            SkewSymTensor3D res;

            // 1. T = (*this) * M
            Tensor3D T;
            const CfemReal m00 = M.data[Tensor3D::idx_xx], m01 = M.data[Tensor3D::idx_xy], m02 = M.data[Tensor3D::idx_xz];
            const CfemReal m10 = M.data[Tensor3D::idx_yx], m11 = M.data[Tensor3D::idx_yy], m12 = M.data[Tensor3D::idx_yz];
            const CfemReal m20 = M.data[Tensor3D::idx_zx], m21 = M.data[Tensor3D::idx_zy], m22 = M.data[Tensor3D::idx_zz];

            // Row 0 of W: [0, xy, xz]
            T.data[Tensor3D::idx_xx] = data[idx_xy] * m10 + data[idx_xz] * m20;
            T.data[Tensor3D::idx_xy] = data[idx_xy] * m11 + data[idx_xz] * m21;
            T.data[Tensor3D::idx_xz] = data[idx_xy] * m12 + data[idx_xz] * m22;

            // Row 1 of W: [-xy, 0, yz]
            T.data[Tensor3D::idx_yx] = -data[idx_xy] * m00 + data[idx_yz] * m20;
            T.data[Tensor3D::idx_yy] = -data[idx_xy] * m01 + data[idx_yz] * m21;
            T.data[Tensor3D::idx_yz] = -data[idx_xy] * m02 + data[idx_yz] * m22;

            // Row 2 of W: [-xz, -yz, 0]
            T.data[Tensor3D::idx_zx] = -data[idx_xz] * m00 - data[idx_yz] * m10;
            T.data[Tensor3D::idx_zy] = -data[idx_xz] * m01 - data[idx_yz] * m11;
            T.data[Tensor3D::idx_zz] = -data[idx_xz] * m02 - data[idx_yz] * m12;

            // 2. res = M^T * T (Only extract the 3 skew components)
            res.data[idx_yz] = m01 * T.data[Tensor3D::idx_xz] + m11 * T.data[Tensor3D::idx_yz] + m21 * T.data[Tensor3D::idx_zz]; // yz
            res.data[idx_xz] = m00 * T.data[Tensor3D::idx_xz] + m10 * T.data[Tensor3D::idx_yz] + m20 * T.data[Tensor3D::idx_zz]; // xz
            res.data[idx_xy] = m00 * T.data[Tensor3D::idx_xy] + m10 * T.data[Tensor3D::idx_yy] + m20 * T.data[Tensor3D::idx_zy]; // xy

            return res;
        }

        // Physical Invariants

        /** @brief Trace is always exactly zero for a skew-symmetric tensor. */
        [[nodiscard]] constexpr CfemReal trace() const noexcept
        {
            return 0.0;
        }

        /** @brief Double contraction with another skew-symmetric tensor. */
        [[nodiscard]] constexpr CfemReal doubleDot(const SkewSymTensor3D &other) const noexcept
        {
            // A:B = 2 * (A_xy*B_xy + A_xz*B_xz + A_yz*B_yz)
            return 2.0 * (data[idx_yz] * other.data[idx_yz] + data[idx_xz] * other.data[idx_xz] + data[idx_xy] * other.data[idx_xy]);
        }

        /** @brief Frobenius norm */
        [[nodiscard]] inline CfemReal norm() const noexcept
        {
            return std::sqrt(this->doubleDot(*this));
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
        [[nodiscard]] constexpr FixedVector<CfemReal, 3> axialVector() const noexcept
        {
            FixedVector<CfemReal, 3> w;
            w.data[0] = -data[idx_yz];
            w.data[1] =  data[idx_xz];
            w.data[2] = -data[idx_xy];
            return w;
        }

        /**
         * @brief Builds a skew-symmetric tensor from its axial vector.
         *
         * @param w Axial vector.
         * @return Associated skew-symmetric tensor.
         */
        [[nodiscard]] static constexpr SkewSymTensor3D fromAxialVector(const FixedVector<CfemReal, 3>& w) noexcept
        {
            SkewSymTensor3D W;
            W.data[idx_xy] = -w.data[2];
            W.data[idx_xz] =  w.data[1];
            W.data[idx_yz] = -w.data[0];
            return W;
        }

        [[nodiscard]] constexpr SkewSymTensor3D transpose() const noexcept
        {
            SkewSymTensor3D W_T;
            W_T.data[idx_xy] = -data[idx_xy];
            W_T.data[idx_xz] = -data[idx_xz];
            W_T.data[idx_yz] = -data[idx_yz];
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
        [[nodiscard]] inline Tensor3D exp() const noexcept
        {
            FixedVector<CfemReal, 3> w = axialVector();
            const CfemReal theta = std::sqrt(w.data[0]*w.data[0] + 
                                             w.data[1]*w.data[1] + 
                                             w.data[2]*w.data[2] );
            
            Tensor3D res = IdentityTensor3D();
            
            if (theta < ZERO_TOLERANCE) {
                // Taylor expansion for small angles
                res(0,1) += data[idx_xy]; res(0,2) += data[idx_xz];
                res(1,0) -= data[idx_xy]; res(1,2) += data[idx_yz];
                res(2,0) -= data[idx_xz]; res(2,1) -= data[idx_yz];
                return res;
            }

            const CfemReal coef1 = std::sin(theta) / theta;
            const CfemReal coef2 = (1.0 - std::cos(theta)) / (theta * theta);

            // Optimisation : compute R = I + coef1 * W + coef2 * W^2 component by component
            res(0,0) += coef2 * (-w.data[1]*w.data[1] - w.data[2]*w.data[2]);
            res(0,1) += coef1 * data[idx_xy] + coef2 * (w.data[0]*w.data[1]);
            res(0,2) += coef1 * data[idx_xz] + coef2 * (w.data[0]*w.data[2]);

            res(1,0) -= coef1 * data[idx_xy] - coef2 * (w.data[0]*w.data[1]);
            res(1,1) += coef2 * (-w.data[0]*w.data[0] - w.data[2]*w.data[2]);
            res(1,2) += coef1 * data[idx_yz] + coef2 * (w.data[1]*w.data[2]);

            res(2,0) -= coef1 * data[idx_xz] - coef2 * (w.data[0]*w.data[2]);
            res(2,1) -= coef1 * data[idx_yz] - coef2 * (w.data[1]*w.data[2]);
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
            return std::sqrt( data[idx_xy]*data[idx_xy] + data[idx_xz]*data[idx_xz] + data[idx_yz]*data[idx_yz] );
        }
        
    };

    // Global Ops
    
    // operator*
    [[nodiscard]] constexpr SymTensor3D operator*(CfemReal scalar, const SymTensor3D &S) noexcept
    {
        return S * scalar;
    }

    [[nodiscard]] constexpr SkewSymTensor3D operator*(CfemReal scalar, const SkewSymTensor3D &W) noexcept
    {
        return W * scalar;
    }

    [[nodiscard]] constexpr AxiSymTensor3D operator*(CfemReal scalar, const AxiSymTensor3D &A) noexcept
    {
        return A * scalar;
    }

    /**
     * @ingroup tensor3d_operators
     * @name Matrix-matrix multiplication operators
     * @brief Matrix-matrix tensor products and algebraic operations.
     * @{
     */

    /** Tensor3D x SymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const Tensor3D &T, const SymTensor3D &S) noexcept
    {
        Tensor3D res;
        const CfemReal Txx = T.data[Tensor3D::idx_xx], Txy = T.data[Tensor3D::idx_xy], Txz = T.data[Tensor3D::idx_xz];
        const CfemReal Tyx = T.data[Tensor3D::idx_yx], Tyy = T.data[Tensor3D::idx_yy], Tyz = T.data[Tensor3D::idx_yz];
        const CfemReal Tzx = T.data[Tensor3D::idx_zx], Tzy = T.data[Tensor3D::idx_zy], Tzz = T.data[Tensor3D::idx_zz];

        const CfemReal Sxx = S.data[SymTensor3D::idx_xx], Syy = S.data[SymTensor3D::idx_yy], Szz = S.data[SymTensor3D::idx_zz], Sxy = S.data[SymTensor3D::idx_xy], Sxz = S.data[SymTensor3D::idx_xz], Syz = S.data[SymTensor3D::idx_yz];

        res.data[Tensor3D::idx_xx] = Txx * Sxx + Txy * Sxy + Txz * Sxz;
        res.data[Tensor3D::idx_xy] = Txx * Sxy + Txy * Syy + Txz * Syz;
        res.data[Tensor3D::idx_xz] = Txx * Sxz + Txy * Syz + Txz * Szz;

        res.data[Tensor3D::idx_yx] = Tyx * Sxx + Tyy * Sxy + Tyz * Sxz;
        res.data[Tensor3D::idx_yy] = Tyx * Sxy + Tyy * Syy + Tyz * Syz;
        res.data[Tensor3D::idx_yz] = Tyx * Sxz + Tyy * Syz + Tyz * Szz;

        res.data[Tensor3D::idx_zx] = Tzx * Sxx + Tzy * Sxy + Tzz * Sxz;
        res.data[Tensor3D::idx_zy] = Tzx * Sxy + Tzy * Syy + Tzz * Syz;
        res.data[Tensor3D::idx_zz] = Tzx * Sxz + Tzy * Syz + Tzz * Szz;

        return res;
    }

    /** Tensor3D x AxiSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const Tensor3D &T, const AxiSymTensor3D &A) noexcept
    {
        Tensor3D res;
        const CfemReal Txx = T.data[Tensor3D::idx_xx], Txy = T.data[Tensor3D::idx_xy], Txz = T.data[Tensor3D::idx_xz];
        const CfemReal Tyx = T.data[Tensor3D::idx_yx], Tyy = T.data[Tensor3D::idx_yy], Tyz = T.data[Tensor3D::idx_yz];
        const CfemReal Tzx = T.data[Tensor3D::idx_zx], Tzy = T.data[Tensor3D::idx_zy], Tzz = T.data[Tensor3D::idx_zz];

        const CfemReal Sxx = A.data[AxiSymTensor3D::idx_xx], Syy = A.data[AxiSymTensor3D::idx_yy], Szz = A.data[AxiSymTensor3D::idx_zz], Sxy = A.data[AxiSymTensor3D::idx_xy];

        res.data[Tensor3D::idx_xx] = Txx * Sxx + Txy * Sxy;
        res.data[Tensor3D::idx_xy] = Txx * Sxy + Txy * Syy;
        res.data[Tensor3D::idx_xz] = Txz * Szz;
        
        res.data[Tensor3D::idx_yx] = Tyx * Sxx + Tyy * Sxy;
        res.data[Tensor3D::idx_yy] = Tyx * Sxy + Tyy * Syy;
        res.data[Tensor3D::idx_yz] = Tyz * Szz;
        
        res.data[Tensor3D::idx_zx] = Tzx * Sxx + Tzy * Sxy;
        res.data[Tensor3D::idx_zy] = Tzx * Sxy + Tzy * Syy;
        res.data[Tensor3D::idx_zz] = Tzz * Szz;
        return res;
    }

    /** Tensor3D x SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const Tensor3D &T, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res;
        const CfemReal Txx = T.data[Tensor3D::idx_xx], Txy = T.data[Tensor3D::idx_xy], Txz = T.data[Tensor3D::idx_xz];
        const CfemReal Tyx = T.data[Tensor3D::idx_yx], Tyy = T.data[Tensor3D::idx_yy], Tyz = T.data[Tensor3D::idx_yz];
        const CfemReal Tzx = T.data[Tensor3D::idx_zx], Tzy = T.data[Tensor3D::idx_zy], Tzz = T.data[Tensor3D::idx_zz];

        const CfemReal Wxy = W.data[SkewSymTensor3D::idx_xy], Wxz = W.data[SkewSymTensor3D::idx_xz], Wyz = W.data[SkewSymTensor3D::idx_yz];

        res.data[Tensor3D::idx_xx] = -Txy * Wxy - Txz * Wxz;
        res.data[Tensor3D::idx_yx] = -Tyy * Wxy - Tyz * Wxz;
        res.data[Tensor3D::idx_zx] = -Tzy * Wxy - Tzz * Wxz;

        res.data[Tensor3D::idx_xy] = Txx * Wxy - Txz * Wyz;
        res.data[Tensor3D::idx_yy] = Tyx * Wxy - Tyz * Wyz;
        res.data[Tensor3D::idx_zy] = Tzx * Wxy - Tzz * Wyz;

        res.data[Tensor3D::idx_xz] = Txx * Wxz + Txy * Wyz;
        res.data[Tensor3D::idx_yz] = Tyx * Wxz + Tyy * Wyz;
        res.data[Tensor3D::idx_zz] = Tzx * Wxz + Tzy * Wyz;

        return res;
    }

    /** AxiSymTensor3D x Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const AxiSymTensor3D &A, const Tensor3D &T) noexcept
    {
        Tensor3D res;
        const CfemReal Txx = T.data[Tensor3D::idx_xx], Txy = T.data[Tensor3D::idx_xy], Txz = T.data[Tensor3D::idx_xz];
        const CfemReal Tyx = T.data[Tensor3D::idx_yx], Tyy = T.data[Tensor3D::idx_yy], Tyz = T.data[Tensor3D::idx_yz];
        const CfemReal Tzx = T.data[Tensor3D::idx_zx], Tzy = T.data[Tensor3D::idx_zy], Tzz = T.data[Tensor3D::idx_zz];

        const CfemReal Sxx = A.data[AxiSymTensor3D::idx_xx], Syy = A.data[AxiSymTensor3D::idx_yy], Szz = A.data[AxiSymTensor3D::idx_zz], Sxy = A.data[AxiSymTensor3D::idx_xy];

        res.data[Tensor3D::idx_xx] = Sxx * Txx + Sxy * Tyx;
        res.data[Tensor3D::idx_xy] = Sxx * Txy + Sxy * Tyy;
        res.data[Tensor3D::idx_xz] = Sxx * Txz + Sxy * Tyz;
        
        res.data[Tensor3D::idx_yx] = Sxy * Txx + Syy * Tyx;
        res.data[Tensor3D::idx_yy] = Sxy * Txy + Syy * Tyy;
        res.data[Tensor3D::idx_yz] = Sxy * Txz + Syy * Tyz;
        
        res.data[Tensor3D::idx_zx] = Szz * Tzx;
        res.data[Tensor3D::idx_zy] = Szz * Tzy;
        res.data[Tensor3D::idx_zz] = Szz * Tzz;
        return res;
    }

    /** AxiSymTensor3D x SymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const AxiSymTensor3D &A, const SymTensor3D &S) noexcept
    {
        Tensor3D res;

        res.data[Tensor3D::idx_xx] = A.data[AxiSymTensor3D::idx_xx]*S.data[SymTensor3D::idx_xx] + A.data[AxiSymTensor3D::idx_xy]*S.data[SymTensor3D::idx_xy];
        res.data[Tensor3D::idx_xy] = A.data[AxiSymTensor3D::idx_xx]*S.data[SymTensor3D::idx_xy] + A.data[AxiSymTensor3D::idx_xy]*S.data[SymTensor3D::idx_yy];
        res.data[Tensor3D::idx_xz] = A.data[AxiSymTensor3D::idx_xx]*S.data[SymTensor3D::idx_xz] + A.data[AxiSymTensor3D::idx_xy]*S.data[SymTensor3D::idx_yz];

        res.data[Tensor3D::idx_yx] = A.data[AxiSymTensor3D::idx_xy]*S.data[SymTensor3D::idx_xx] + A.data[AxiSymTensor3D::idx_yy]*S.data[SymTensor3D::idx_xy];
        res.data[Tensor3D::idx_yy] = A.data[AxiSymTensor3D::idx_xy]*S.data[SymTensor3D::idx_xy] + A.data[AxiSymTensor3D::idx_yy]*S.data[SymTensor3D::idx_yy];
        res.data[Tensor3D::idx_yz] = A.data[AxiSymTensor3D::idx_xy]*S.data[SymTensor3D::idx_xz] + A.data[AxiSymTensor3D::idx_yy]*S.data[SymTensor3D::idx_yz];

        res.data[Tensor3D::idx_zx] = A.data[AxiSymTensor3D::idx_zz]*S.data[SymTensor3D::idx_xz];
        res.data[Tensor3D::idx_zy] = A.data[AxiSymTensor3D::idx_zz]*S.data[SymTensor3D::idx_yz];
        res.data[Tensor3D::idx_zz] = A.data[AxiSymTensor3D::idx_zz]*S.data[SymTensor3D::idx_zz];

        return res;
    }

    /** AxiSymTensor3D x SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const AxiSymTensor3D &A, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res;
        const CfemReal Wxy = W.data[SkewSymTensor3D::idx_xy], Wxz = W.data[SkewSymTensor3D::idx_xz], Wyz = W.data[SkewSymTensor3D::idx_yz];
        const CfemReal Sxx = A.data[AxiSymTensor3D::idx_xx], Syy = A.data[AxiSymTensor3D::idx_yy], Szz = A.data[AxiSymTensor3D::idx_zz], Sxy = A.data[AxiSymTensor3D::idx_xy];

        res.data[Tensor3D::idx_xx] =-Sxy * Wxy;
        res.data[Tensor3D::idx_xy] = Sxx * Wxy;
        res.data[Tensor3D::idx_xz] = Sxx * Wxz + Sxy * Wyz;

        res.data[Tensor3D::idx_yx] =-Syy * Wxy;
        res.data[Tensor3D::idx_yy] = Sxy * Wxy;
        res.data[Tensor3D::idx_yz] = Sxy * Wxz + Syy * Wyz;

        res.data[Tensor3D::idx_zx] =-Szz * Wxz;
        res.data[Tensor3D::idx_zy] =-Szz * Wyz;
        res.data[Tensor3D::idx_zz] = 0.0;

        return res;
    }

    /** SymTensor3D x Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const SymTensor3D &S, const Tensor3D &T) noexcept
    {
        Tensor3D res;
        const CfemReal Txx = T.data[Tensor3D::idx_xx], Txy = T.data[Tensor3D::idx_xy], Txz = T.data[Tensor3D::idx_xz];
        const CfemReal Tyx = T.data[Tensor3D::idx_yx], Tyy = T.data[Tensor3D::idx_yy], Tyz = T.data[Tensor3D::idx_yz];
        const CfemReal Tzx = T.data[Tensor3D::idx_zx], Tzy = T.data[Tensor3D::idx_zy], Tzz = T.data[Tensor3D::idx_zz];

        const CfemReal Sxx = S.data[SymTensor3D::idx_xx], Syy = S.data[SymTensor3D::idx_yy], Szz = S.data[SymTensor3D::idx_zz], Sxy = S.data[SymTensor3D::idx_xy], Sxz = S.data[SymTensor3D::idx_xz], Syz = S.data[SymTensor3D::idx_yz];

        res.data[Tensor3D::idx_xx] = Sxx * Txx + Sxy * Tyx + Sxz * Tzx;
        res.data[Tensor3D::idx_xy] = Sxx * Txy + Sxy * Tyy + Sxz * Tzy;
        res.data[Tensor3D::idx_xz] = Sxx * Txz + Sxy * Tyz + Sxz * Tzz;

        res.data[Tensor3D::idx_yx] = Sxy * Txx + Syy * Tyx + Syz * Tzx;
        res.data[Tensor3D::idx_yy] = Sxy * Txy + Syy * Tyy + Syz * Tzy;
        res.data[Tensor3D::idx_yz] = Sxy * Txz + Syy * Tyz + Syz * Tzz;

        res.data[Tensor3D::idx_zx] = Sxz * Txx + Syz * Tyx + Szz * Tzx;
        res.data[Tensor3D::idx_zy] = Sxz * Txy + Syz * Tyy + Szz * Tzy;
        res.data[Tensor3D::idx_zz] = Sxz * Txz + Syz * Tyz + Szz * Tzz;

        return res;
    }

    /** SymTensor3D x SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const SymTensor3D &S, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res;
        const CfemReal Sxx = S.data[SymTensor3D::idx_xx], Syy = S.data[SymTensor3D::idx_yy], Szz = S.data[SymTensor3D::idx_zz], Sxy = S.data[SymTensor3D::idx_xy], Sxz = S.data[SymTensor3D::idx_xz], Syz = S.data[SymTensor3D::idx_yz];
        const CfemReal Wxy = W.data[SkewSymTensor3D::idx_xy], Wxz = W.data[SkewSymTensor3D::idx_xz], Wyz = W.data[SkewSymTensor3D::idx_yz];

        res.data[Tensor3D::idx_xx] =-Sxy * Wxy - Sxz * Wxz;
        res.data[Tensor3D::idx_xy] = Sxx * Wxy - Sxz * Wyz;
        res.data[Tensor3D::idx_xz] = Sxx * Wxz + Sxy * Wyz;

        res.data[Tensor3D::idx_yx] =-Syy * Wxy - Syz * Wxz;
        res.data[Tensor3D::idx_yy] = Sxy * Wxy - Syz * Wyz;
        res.data[Tensor3D::idx_yz] = Sxy * Wxz + Syy * Wyz;

        res.data[Tensor3D::idx_zx] =-Syz * Wxy - Szz * Wxz;
        res.data[Tensor3D::idx_zy] = Sxz * Wxy - Szz * Wyz;
        res.data[Tensor3D::idx_zz] = Sxz * Wxz + Syz * Wyz;

        return res;
    }

    /** SymTensor3D x AxiSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const SymTensor3D &S, const AxiSymTensor3D &A) noexcept
    {
        Tensor3D res;
        const CfemReal Sxx = S.data[SymTensor3D::idx_xx], Syy = S.data[SymTensor3D::idx_yy], Szz = S.data[SymTensor3D::idx_zz], Sxy = S.data[SymTensor3D::idx_xy], Sxz = S.data[SymTensor3D::idx_xz], Syz = S.data[SymTensor3D::idx_yz];
        const CfemReal Axx = A.data[AxiSymTensor3D::idx_xx], Ayy = A.data[AxiSymTensor3D::idx_yy], Azz = A.data[AxiSymTensor3D::idx_zz], Axy = A.data[AxiSymTensor3D::idx_xy];

        res.data[Tensor3D::idx_xx] = Sxx * Axx + Sxy * Axy;
        res.data[Tensor3D::idx_xy] = Sxx * Axy + Sxy * Ayy;
        res.data[Tensor3D::idx_xz] = Sxz * Azz;

        res.data[Tensor3D::idx_yx] = Sxy * Axx + Syy * Axy;
        res.data[Tensor3D::idx_yy] = Sxy * Axy + Syy * Ayy;
        res.data[Tensor3D::idx_yz] = Syz * Azz;

        res.data[Tensor3D::idx_zx] = Sxz * Axx + Syz * Axy;
        res.data[Tensor3D::idx_zy] = Sxz * Axy + Syz * Ayy;
        res.data[Tensor3D::idx_zz] = Szz * Azz;

        return res;
    }

    /** SkewSymTensor3D x Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const SkewSymTensor3D &W, const Tensor3D &T) noexcept
    {
        Tensor3D res;
        const CfemReal Txx = T.data[Tensor3D::idx_xx], Txy = T.data[Tensor3D::idx_xy], Txz = T.data[Tensor3D::idx_xz];
        const CfemReal Tyx = T.data[Tensor3D::idx_yx], Tyy = T.data[Tensor3D::idx_yy], Tyz = T.data[Tensor3D::idx_yz];
        const CfemReal Tzx = T.data[Tensor3D::idx_zx], Tzy = T.data[Tensor3D::idx_zy], Tzz = T.data[Tensor3D::idx_zz];

        const CfemReal Wxy = W.data[SkewSymTensor3D::idx_xy], Wxz = W.data[SkewSymTensor3D::idx_xz], Wyz = W.data[SkewSymTensor3D::idx_yz];

        res.data[Tensor3D::idx_xx] = Wxy * Tyx + Wxz * Tzx;
        res.data[Tensor3D::idx_xy] = Wxy * Tyy + Wxz * Tzy;
        res.data[Tensor3D::idx_xz] = Wxy * Tyz + Wxz * Tzz;

        res.data[Tensor3D::idx_yx] = -Wxy * Txx + Wyz * Tzx;
        res.data[Tensor3D::idx_yy] = -Wxy * Txy + Wyz * Tzy;
        res.data[Tensor3D::idx_yz] = -Wxy * Txz + Wyz * Tzz;

        res.data[Tensor3D::idx_zx] = -Wxz * Txx - Wyz * Tyx;
        res.data[Tensor3D::idx_zy] = -Wxz * Txy - Wyz * Tyy;
        res.data[Tensor3D::idx_zz] = -Wxz * Txz - Wyz * Tyz;

        return res;
    }

    /** SkewSymTensor3D x SymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const SkewSymTensor3D &W, const SymTensor3D &S) noexcept
    {
        Tensor3D res;
        const CfemReal Sxx = S.data[SymTensor3D::idx_xx], Syy = S.data[SymTensor3D::idx_yy], Szz = S.data[SymTensor3D::idx_zz], Sxy = S.data[SymTensor3D::idx_xy], Sxz = S.data[SymTensor3D::idx_xz], Syz = S.data[SymTensor3D::idx_yz];
        const CfemReal Wxy = W.data[SkewSymTensor3D::idx_xy], Wxz = W.data[SkewSymTensor3D::idx_xz], Wyz = W.data[SkewSymTensor3D::idx_yz];

        res.data[Tensor3D::idx_xx] = Wxy * Sxy + Wxz * Sxz;
        res.data[Tensor3D::idx_xy] = Wxy * Syy + Wxz * Syz;
        res.data[Tensor3D::idx_xz] = Wxy * Syz + Wxz * Szz;

        res.data[Tensor3D::idx_yx] = -Wxy * Sxx + Wyz * Sxz;
        res.data[Tensor3D::idx_yy] = -Wxy * Sxy + Wyz * Syz;
        res.data[Tensor3D::idx_yz] = -Wxy * Sxz + Wyz * Szz;

        res.data[Tensor3D::idx_zx] = -Wxz * Sxx - Wyz * Sxy;
        res.data[Tensor3D::idx_zy] = -Wxz * Sxy - Wyz * Syy;
        res.data[Tensor3D::idx_zz] = -Wxz * Sxz - Wyz * Syz;

        return res;
    }

    /** SkewSymTensor3D x AxiSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const SkewSymTensor3D &W, const AxiSymTensor3D &A) noexcept
    {
        Tensor3D res;
        const CfemReal Sxx = A.data[AxiSymTensor3D::idx_xx], Syy = A.data[AxiSymTensor3D::idx_yy], Szz = A.data[AxiSymTensor3D::idx_zz], Sxy = A.data[AxiSymTensor3D::idx_xy];
        const CfemReal Wxy = W.data[SkewSymTensor3D::idx_xy], Wxz = W.data[SkewSymTensor3D::idx_xz], Wyz = W.data[SkewSymTensor3D::idx_yz];

        res.data[Tensor3D::idx_xx] =  Wxy * Sxy;             res.data[Tensor3D::idx_xy] = Wxy * Syy;             res.data[Tensor3D::idx_xz] = Wxz * Szz;
        res.data[Tensor3D::idx_yx] = -Wxy * Sxx;             res.data[Tensor3D::idx_yy] =-Wxy * Sxy;             res.data[Tensor3D::idx_yz] = Wyz * Szz;
        res.data[Tensor3D::idx_zx] = -Wxz * Sxx - Wyz * Sxy; res.data[Tensor3D::idx_zy] =-Wxz * Sxy - Wyz * Syy; res.data[Tensor3D::idx_zz] = 0.0;

        return res;
    }
    /** @} */

    /**
     * @ingroup tensor3d_operators
     * @name Matrix-matrix addition operators
     * @brief Matrix-matrix tensor addition
     * @{
     */

    /** Tensor3D + SymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const Tensor3D &T, const SymTensor3D &S) noexcept
    {
        Tensor3D res = T;
        res.data[Tensor3D::idx_xx] += S.data[SymTensor3D::idx_xx]; res.data[Tensor3D::idx_xy] += S.data[SymTensor3D::idx_xy]; res.data[Tensor3D::idx_xz] += S.data[SymTensor3D::idx_xz];
        res.data[Tensor3D::idx_yx] += S.data[SymTensor3D::idx_xy];                   res.data[Tensor3D::idx_yy] += S.data[SymTensor3D::idx_yy];                   res.data[Tensor3D::idx_yz] += S.data[SymTensor3D::idx_yz];
        res.data[Tensor3D::idx_zx] += S.data[SymTensor3D::idx_xz];                   res.data[Tensor3D::idx_zy] += S.data[SymTensor3D::idx_yz];                   res.data[Tensor3D::idx_zz] += S.data[SymTensor3D::idx_zz];
        return res;
    }

    /** Tensor3D + AxiSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const Tensor3D &T, const AxiSymTensor3D &A) noexcept
    {
        Tensor3D res = T;
        res.data[Tensor3D::idx_xx] += A.data[AxiSymTensor3D::idx_xx]; res.data[Tensor3D::idx_xy] += A.data[AxiSymTensor3D::idx_xy];
        res.data[Tensor3D::idx_yx] += A.data[AxiSymTensor3D::idx_xy]; res.data[Tensor3D::idx_yy] += A.data[AxiSymTensor3D::idx_yy]; 
        res.data[Tensor3D::idx_zz] += A.data[AxiSymTensor3D::idx_zz];
        return res;
    }

    /** Tensor3D + SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const Tensor3D &T, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res = T;
        res.data[Tensor3D::idx_xy] += W.data[SkewSymTensor3D::idx_xy];  res.data[Tensor3D::idx_xz] += W.data[SkewSymTensor3D::idx_xz];
        res.data[Tensor3D::idx_yx] -= W.data[SkewSymTensor3D::idx_xy];  res.data[Tensor3D::idx_yz] += W.data[SkewSymTensor3D::idx_yz];
        res.data[Tensor3D::idx_zx] -= W.data[SkewSymTensor3D::idx_xz];  res.data[Tensor3D::idx_zy] -= W.data[SkewSymTensor3D::idx_yz];
        return res;
    }

    /** SymTensor3D + Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const SymTensor3D &S, const Tensor3D &T) noexcept
    {
        return T + S;
    }

    /** SymTensor3D + AxiSymTensor3D */
    [[nodiscard]] constexpr SymTensor3D operator+(const SymTensor3D &S, const AxiSymTensor3D &A) noexcept
    {
        SymTensor3D res = S;
        res.data[SymTensor3D::idx_xx] += A.data[AxiSymTensor3D::idx_xx]; 
        res.data[SymTensor3D::idx_yy] += A.data[AxiSymTensor3D::idx_yy]; 
        res.data[SymTensor3D::idx_zz] += A.data[AxiSymTensor3D::idx_zz]; 
        res.data[SymTensor3D::idx_xy] += A.data[AxiSymTensor3D::idx_xy];
        return res;
    }

    /** SymTensor3D + SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const SymTensor3D &S, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res;
        res.data[Tensor3D::idx_xx] = S.data[SymTensor3D::idx_xx]; 
        res.data[Tensor3D::idx_xy] = S.data[SymTensor3D::idx_xy] + W.data[SkewSymTensor3D::idx_xy]; 
        res.data[Tensor3D::idx_xz] = S.data[SymTensor3D::idx_xz] + W.data[SkewSymTensor3D::idx_xz];

        res.data[Tensor3D::idx_yx] = S.data[SymTensor3D::idx_xy] - W.data[SkewSymTensor3D::idx_xy]; 
        res.data[Tensor3D::idx_yy] = S.data[SymTensor3D::idx_yy];          
        res.data[Tensor3D::idx_yz] = S.data[SymTensor3D::idx_yz] + W.data[SkewSymTensor3D::idx_yz];

        res.data[Tensor3D::idx_zx] = S.data[SymTensor3D::idx_xz] - W.data[SkewSymTensor3D::idx_xz]; 
        res.data[Tensor3D::idx_zy] = S.data[SymTensor3D::idx_yz] - W.data[SkewSymTensor3D::idx_yz]; 
        res.data[Tensor3D::idx_zz] = S.data[SymTensor3D::idx_zz];
        return res;
    }

    /** AxiSymTensor3D + Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const AxiSymTensor3D &A, const Tensor3D &T) noexcept
    {
        return T + A;
    }

    /** AxiSymTensor3D + SymTensor3D */
    [[nodiscard]] constexpr SymTensor3D operator+(const AxiSymTensor3D &A, const SymTensor3D &S) noexcept
    {
        return S + A;
    }

    /** AxiSymTensor3D + SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const AxiSymTensor3D &A, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res;
        res.data[Tensor3D::idx_xx] = A.data[AxiSymTensor3D::idx_xx];
        res.data[Tensor3D::idx_xy] = A.data[AxiSymTensor3D::idx_xy] + W.data[SkewSymTensor3D::idx_xy]; 
        res.data[Tensor3D::idx_xz] =                                 W.data[SkewSymTensor3D::idx_xz];
        res.data[Tensor3D::idx_yx] = A.data[AxiSymTensor3D::idx_xy] - W.data[SkewSymTensor3D::idx_xy]; 
        res.data[Tensor3D::idx_yy] = A.data[AxiSymTensor3D::idx_yy];            
        res.data[Tensor3D::idx_yz] =                                 W.data[SkewSymTensor3D::idx_yz];
        res.data[Tensor3D::idx_zx] =                               - W.data[SkewSymTensor3D::idx_xz]; 
        res.data[Tensor3D::idx_zy] =                               - W.data[SkewSymTensor3D::idx_yz];                        
        res.data[Tensor3D::idx_zz] = A.data[AxiSymTensor3D::idx_zz];
        return res;
    }

    /** SkewSymTensor3D + Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const SkewSymTensor3D &W, const Tensor3D &T) noexcept
    {
        return T + W;
    }

    /** SkewSymTensor3D + SymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const SkewSymTensor3D &W, const SymTensor3D &S) noexcept
    {
        return S + W; 
    }

    /** SkewSymTensor3D + AxiSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const SkewSymTensor3D &W, const AxiSymTensor3D &A) noexcept
    {
        return A + W;
    }
    /** @} */

    /**
     * @ingroup tensor3d_operators
     * @name Matrix-matrix substraction operators
     * @brief Matrix-matrix tensor subtraction
     * @{
     */

    /** Tensor3D - SymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const Tensor3D &T, const SymTensor3D &S) noexcept
    {
        Tensor3D res = T;
        res.data[Tensor3D::idx_xx] -= S.data[SymTensor3D::idx_xx]; 
        res.data[Tensor3D::idx_xy] -= S.data[SymTensor3D::idx_xy]; 
        res.data[Tensor3D::idx_xz] -= S.data[SymTensor3D::idx_xz];

        res.data[Tensor3D::idx_yx] -= S.data[SymTensor3D::idx_xy]; 
        res.data[Tensor3D::idx_yy] -= S.data[SymTensor3D::idx_yy]; 
        res.data[Tensor3D::idx_yz] -= S.data[SymTensor3D::idx_yz];

        res.data[Tensor3D::idx_zx] -= S.data[SymTensor3D::idx_xz]; 
        res.data[Tensor3D::idx_zy] -= S.data[SymTensor3D::idx_yz]; 
        res.data[Tensor3D::idx_zz] -= S.data[SymTensor3D::idx_zz];
        return res;
    }

    /** Tensor3D - AxiSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const Tensor3D &T, const AxiSymTensor3D &A) noexcept
    {
        Tensor3D res = T;
        res.data[Tensor3D::idx_xx] -= A.data[AxiSymTensor3D::idx_xx];
        res.data[Tensor3D::idx_xy] -= A.data[AxiSymTensor3D::idx_xy];
        res.data[Tensor3D::idx_yx] -= A.data[AxiSymTensor3D::idx_xy]; 
        res.data[Tensor3D::idx_yy] -= A.data[AxiSymTensor3D::idx_yy];
        res.data[Tensor3D::idx_zz] -= A.data[AxiSymTensor3D::idx_zz];
        return res;
    }

    /** Tensor3D - SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const Tensor3D &T, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res = T;
        res.data[Tensor3D::idx_xy] -= W.data[SkewSymTensor3D::idx_xy]; res.data[Tensor3D::idx_xz] -= W.data[SkewSymTensor3D::idx_xz];
        res.data[Tensor3D::idx_yx] += W.data[SkewSymTensor3D::idx_xy]; res.data[Tensor3D::idx_yz] -= W.data[SkewSymTensor3D::idx_yz];
        res.data[Tensor3D::idx_zx] += W.data[SkewSymTensor3D::idx_xz]; res.data[Tensor3D::idx_zy] += W.data[SkewSymTensor3D::idx_yz];

        return res;
    }

    /** SymTensor3D - Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const SymTensor3D &S, const Tensor3D &T) noexcept
    {
        Tensor3D res;
        res.data[Tensor3D::idx_xx] = S.data[SymTensor3D::idx_xx] - T.data[Tensor3D::idx_xx]; 
        res.data[Tensor3D::idx_xy] = S.data[SymTensor3D::idx_xy] - T.data[Tensor3D::idx_xy]; 
        res.data[Tensor3D::idx_xz] = S.data[SymTensor3D::idx_xz] - T.data[Tensor3D::idx_xz];

        res.data[Tensor3D::idx_yx] = S.data[SymTensor3D::idx_xy] - T.data[Tensor3D::idx_yx]; 
        res.data[Tensor3D::idx_yy] = S.data[SymTensor3D::idx_yy] - T.data[Tensor3D::idx_yy]; 
        res.data[Tensor3D::idx_yz] = S.data[SymTensor3D::idx_yz] - T.data[Tensor3D::idx_yz];

        res.data[Tensor3D::idx_zx] = S.data[SymTensor3D::idx_xz] - T.data[Tensor3D::idx_zx]; 
        res.data[Tensor3D::idx_zy] = S.data[SymTensor3D::idx_yz] - T.data[Tensor3D::idx_zy]; 
        res.data[Tensor3D::idx_zz] = S.data[SymTensor3D::idx_zz] - T.data[Tensor3D::idx_zz];
        return res;
    }

    /** SymTensor3D - SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const SymTensor3D &S, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res;
        res.data[Tensor3D::idx_xx] = S.data[SymTensor3D::idx_xx];
        res.data[Tensor3D::idx_xy] = S.data[SymTensor3D::idx_xy] - W.data[SkewSymTensor3D::idx_xy]; 
        res.data[Tensor3D::idx_xz] = S.data[SymTensor3D::idx_xz] - W.data[SkewSymTensor3D::idx_xz];

        res.data[Tensor3D::idx_yx] = S.data[SymTensor3D::idx_xy] + W.data[SkewSymTensor3D::idx_xy]; 
        res.data[Tensor3D::idx_yy] = S.data[SymTensor3D::idx_yy]; 
        res.data[Tensor3D::idx_yz] = S.data[SymTensor3D::idx_yz] - W.data[SkewSymTensor3D::idx_yz];

        res.data[Tensor3D::idx_zx] = S.data[SymTensor3D::idx_xz] + W.data[SkewSymTensor3D::idx_xz]; 
        res.data[Tensor3D::idx_zy] = S.data[SymTensor3D::idx_yz] + W.data[SkewSymTensor3D::idx_yz];
        res.data[Tensor3D::idx_zz] = S.data[SymTensor3D::idx_zz];
        return res;
    }

    /** SymTensor3D - AxiSymTensor3D */
    [[nodiscard]] constexpr SymTensor3D operator-(const SymTensor3D &S, const AxiSymTensor3D &A) noexcept
    {
        return { S.data[SymTensor3D::idx_xx] - A.data[AxiSymTensor3D::idx_xx], 
                 S.data[SymTensor3D::idx_yy] - A.data[AxiSymTensor3D::idx_yy], 
                 S.data[SymTensor3D::idx_zz] - A.data[AxiSymTensor3D::idx_zz], 
                 S.data[SymTensor3D::idx_yz], 
                 S.data[SymTensor3D::idx_xz], 
                 S.data[SymTensor3D::idx_xy] - A.data[AxiSymTensor3D::idx_xy] };
    }

    /** AxiSymTensor3D - Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const AxiSymTensor3D &A, const Tensor3D &T) noexcept
    {
        Tensor3D res;
        res.data[Tensor3D::idx_xx] = A.data[AxiSymTensor3D::idx_xx] - T.data[Tensor3D::idx_xx]; 
        res.data[Tensor3D::idx_xy] = A.data[AxiSymTensor3D::idx_xy] - T.data[Tensor3D::idx_xy]; 
        res.data[Tensor3D::idx_xz] =                               - T.data[Tensor3D::idx_xz];
        res.data[Tensor3D::idx_yx] = A.data[AxiSymTensor3D::idx_xy] - T.data[Tensor3D::idx_yx]; 
        res.data[Tensor3D::idx_yy] = A.data[AxiSymTensor3D::idx_yy] - T.data[Tensor3D::idx_yy]; 
        res.data[Tensor3D::idx_yz] =                               - T.data[Tensor3D::idx_yz];
        res.data[Tensor3D::idx_zx] =                               - T.data[Tensor3D::idx_zx];
        res.data[Tensor3D::idx_zy] =                               - T.data[Tensor3D::idx_zy]; 
        res.data[Tensor3D::idx_zz] = A.data[AxiSymTensor3D::idx_zz] - T.data[Tensor3D::idx_zz];
        return res;
    }

    /** AxiSymTensor3D - SymTensor3D */
    [[nodiscard]] constexpr SymTensor3D operator-(const AxiSymTensor3D &A, const SymTensor3D &S) noexcept
    {
        return { A.data[AxiSymTensor3D::idx_xx] - S.data[SymTensor3D::idx_xx], 
                 A.data[AxiSymTensor3D::idx_yy] - S.data[SymTensor3D::idx_yy], 
                 A.data[AxiSymTensor3D::idx_zz] - S.data[SymTensor3D::idx_zz], 
                                               - S.data[SymTensor3D::idx_yz], 
                                               - S.data[SymTensor3D::idx_xz], 
                 A.data[AxiSymTensor3D::idx_xy] - S.data[SymTensor3D::idx_xy] };
    }

    /** AxiSymTensor3D - SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const AxiSymTensor3D &A, const SkewSymTensor3D &W) noexcept
    {
        return { A.data[AxiSymTensor3D::idx_xx],
                 A.data[AxiSymTensor3D::idx_xy] - W.data[SkewSymTensor3D::idx_xy],
                                               - W.data[SkewSymTensor3D::idx_xz],
                 A.data[AxiSymTensor3D::idx_xy] + W.data[SkewSymTensor3D::idx_xy], 
                 A.data[AxiSymTensor3D::idx_yy],                  
                                               - W.data[SkewSymTensor3D::idx_yz],
                                                 W.data[SkewSymTensor3D::idx_xz],
                                                 W.data[SkewSymTensor3D::idx_yz],  
                 A.data[AxiSymTensor3D::idx_zz] };
    }

    /** SkewSymTensor3D - Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const SkewSymTensor3D &W, const Tensor3D &T) noexcept
    {
        Tensor3D res;
        res.data[Tensor3D::idx_xx] =                                 - T.data[Tensor3D::idx_xx]; 
        res.data[Tensor3D::idx_xy] = W.data[SkewSymTensor3D::idx_xy] - T.data[Tensor3D::idx_xy]; 
        res.data[Tensor3D::idx_xz] = W.data[SkewSymTensor3D::idx_xz] - T.data[Tensor3D::idx_xz];
        res.data[Tensor3D::idx_yx] =-W.data[SkewSymTensor3D::idx_xy] - T.data[Tensor3D::idx_yx]; 
        res.data[Tensor3D::idx_yy] =                                 - T.data[Tensor3D::idx_yy]; 
        res.data[Tensor3D::idx_yz] = W.data[SkewSymTensor3D::idx_yz] - T.data[Tensor3D::idx_yz];
        res.data[Tensor3D::idx_zx] =-W.data[SkewSymTensor3D::idx_xz] - T.data[Tensor3D::idx_zx]; 
        res.data[Tensor3D::idx_zy] =-W.data[SkewSymTensor3D::idx_yz] - T.data[Tensor3D::idx_zy]; 
        res.data[Tensor3D::idx_zz] =                                 - T.data[Tensor3D::idx_zz];
        return res;
    }

    /** SkewSymTensor3D - SymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const SkewSymTensor3D &W, const SymTensor3D &S) noexcept
    {
        Tensor3D res;
        res.data[Tensor3D::idx_xx] =                                  - S.data[SymTensor3D::idx_xx];
        res.data[Tensor3D::idx_xy] =  W.data[SkewSymTensor3D::idx_xy] - S.data[SymTensor3D::idx_xy]; 
        res.data[Tensor3D::idx_xz] =  W.data[SkewSymTensor3D::idx_xz] - S.data[SymTensor3D::idx_xz];
        res.data[Tensor3D::idx_yx] = -W.data[SkewSymTensor3D::idx_xy] - S.data[SymTensor3D::idx_xy]; 
        res.data[Tensor3D::idx_yy] =                                  - S.data[SymTensor3D::idx_yy];
        res.data[Tensor3D::idx_yz] =  W.data[SkewSymTensor3D::idx_yz] - S.data[SymTensor3D::idx_yz];
        res.data[Tensor3D::idx_zx] = -W.data[SkewSymTensor3D::idx_xz] - S.data[SymTensor3D::idx_xz]; 
        res.data[Tensor3D::idx_zy] = -W.data[SkewSymTensor3D::idx_yz] - S.data[SymTensor3D::idx_yz]; 
        res.data[Tensor3D::idx_zz] =                                  - S.data[SymTensor3D::idx_zz];
        return res;
    }

    /** SkewSymTensor3D - AxiSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const SkewSymTensor3D &W, const AxiSymTensor3D &A) noexcept
    {
        Tensor3D res;
        res.data[Tensor3D::idx_xx] =                                  - A.data[AxiSymTensor3D::idx_xx]; 
        res.data[Tensor3D::idx_xy] =  W.data[SkewSymTensor3D::idx_xy] - A.data[AxiSymTensor3D::idx_xy];
        res.data[Tensor3D::idx_xz] =  W.data[SkewSymTensor3D::idx_xz];
        res.data[Tensor3D::idx_yx] = -W.data[SkewSymTensor3D::idx_xy] - A.data[AxiSymTensor3D::idx_xy];
        res.data[Tensor3D::idx_yy] =                                  - A.data[AxiSymTensor3D::idx_yy];
        res.data[Tensor3D::idx_yz] =  W.data[SkewSymTensor3D::idx_yz];
        res.data[Tensor3D::idx_zx] = -W.data[SkewSymTensor3D::idx_xz];
        res.data[Tensor3D::idx_zy] = -W.data[SkewSymTensor3D::idx_yz];
        res.data[Tensor3D::idx_zz] =                                  - A.data[AxiSymTensor3D::idx_zz];
        return res;
    }
    /** @} */

    [[nodiscard]] constexpr CfemReal det(const SymTensor3D &S) noexcept
    {
        return S.data[SymTensor3D::idx_xx] * (S.data[SymTensor3D::idx_yy] * S.data[SymTensor3D::idx_zz] - S.data[SymTensor3D::idx_yz] * S.data[SymTensor3D::idx_yz])
             - S.data[SymTensor3D::idx_xy] * (S.data[SymTensor3D::idx_xy] * S.data[SymTensor3D::idx_zz] - S.data[SymTensor3D::idx_yz] * S.data[SymTensor3D::idx_xz])
             + S.data[SymTensor3D::idx_xz] * (S.data[SymTensor3D::idx_xy] * S.data[SymTensor3D::idx_yz] - S.data[SymTensor3D::idx_yy] * S.data[SymTensor3D::idx_xz]);
    }

    [[nodiscard]] constexpr CfemReal det(const AxiSymTensor3D &A) noexcept
    {
        return (  A.data[AxiSymTensor3D::idx_xx] * A.data[AxiSymTensor3D::idx_yy]
                - A.data[AxiSymTensor3D::idx_xy] * A.data[AxiSymTensor3D::idx_xy]) * A.data[AxiSymTensor3D::idx_zz];
    }

    [[nodiscard]] constexpr CfemReal det(const SkewSymTensor3D &W) noexcept
    {
        return 0.0;
    }

    [[nodiscard]] constexpr SymTensor3D inverse(const SymTensor3D &S, CfemReal determinant)
    {
        CFEM_ASSERT(std::abs(determinant) > cfem::ZERO_TOLERANCE, "Singular 3x3 matrix/Jacobian detected.");

        const CfemReal inv_det = 1.0 / determinant;
        
        SymTensor3D res;
        res.data[SymTensor3D::idx_xx] = ( S.data[SymTensor3D::idx_yy] * S.data[SymTensor3D::idx_zz] - S.data[SymTensor3D::idx_yz] * S.data[SymTensor3D::idx_yz] ) * inv_det;
        res.data[SymTensor3D::idx_yy] = ( S.data[SymTensor3D::idx_xx] * S.data[SymTensor3D::idx_zz] - S.data[SymTensor3D::idx_xz] * S.data[SymTensor3D::idx_xz] ) * inv_det;
        res.data[SymTensor3D::idx_zz] = ( S.data[SymTensor3D::idx_xx] * S.data[SymTensor3D::idx_yy] - S.data[SymTensor3D::idx_xy] * S.data[SymTensor3D::idx_xy] ) * inv_det;
        res.data[SymTensor3D::idx_xy] = ( S.data[SymTensor3D::idx_xz] * S.data[SymTensor3D::idx_yz] - S.data[SymTensor3D::idx_xy] * S.data[SymTensor3D::idx_zz] ) * inv_det;
        res.data[SymTensor3D::idx_xz] = ( S.data[SymTensor3D::idx_xy] * S.data[SymTensor3D::idx_yz] - S.data[SymTensor3D::idx_xx] * S.data[SymTensor3D::idx_xz] ) * inv_det;
        res.data[SymTensor3D::idx_yz] = ( S.data[SymTensor3D::idx_xy] * S.data[SymTensor3D::idx_xz] - S.data[SymTensor3D::idx_yy] * S.data[SymTensor3D::idx_xz] ) * inv_det;
        return res;
    }

    [[nodiscard]] inline SymTensor3D inverse(const SymTensor3D &S)
    {
        return inverse(S, det(S) );
    }

    [[nodiscard]] constexpr AxiSymTensor3D inverse(const AxiSymTensor3D &S, CfemReal determinant)
    {
        CFEM_ASSERT(std::abs(determinant) > cfem::ZERO_TOLERANCE, "Singular 3x3 matrix/Jacobian detected.");

        const CfemReal inv_det = 1.0 / determinant;
        
        AxiSymTensor3D res;
        res.data[AxiSymTensor3D::idx_xx] = ( S.data[AxiSymTensor3D::idx_yy] * S.data[AxiSymTensor3D::idx_zz] ) * inv_det;
        res.data[AxiSymTensor3D::idx_yy] = ( S.data[AxiSymTensor3D::idx_xx] * S.data[AxiSymTensor3D::idx_zz] ) * inv_det;
        res.data[AxiSymTensor3D::idx_zz] = ( S.data[AxiSymTensor3D::idx_xx] * S.data[AxiSymTensor3D::idx_yy] - S.data[AxiSymTensor3D::idx_xy] * S.data[AxiSymTensor3D::idx_xy] ) * inv_det;
        res.data[AxiSymTensor3D::idx_xy] =-( S.data[AxiSymTensor3D::idx_xy] * S.data[AxiSymTensor3D::idx_zz] ) * inv_det;
        return res;
    }

    [[nodiscard]] inline AxiSymTensor3D inverse(const AxiSymTensor3D &S)
    {
        return inverse(S, det(S) );
    }

    /**
     * @brief Tries to inverse a skew-symmetric second order tensor (This operation is impossible).
     * @throw std::runtime_error.
     */
    [[nodiscard]] inline SkewSymTensor3D inverse(const SkewSymTensor3D &W) 
    {
        CFEM_THROW("3x3 skew-symmetric tensor is not invertible since det = 0.");
        return SkewSymTensor3D{}; 
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
    [[nodiscard]] inline FixedVector<CfemReal, 3> eigenvalues(const SymTensor3D &S) noexcept
    {
        const CfemReal I1 = S.data[SymTensor3D::idx_xx] + S.data[SymTensor3D::idx_yy] + S.data[SymTensor3D::idx_zz];
        const CfemReal I2 = S.data[SymTensor3D::idx_xx]*S.data[SymTensor3D::idx_yy] + S.data[SymTensor3D::idx_yy]*S.data[SymTensor3D::idx_zz] + S.data[SymTensor3D::idx_xx]*S.data[SymTensor3D::idx_zz] 
                          - S.data[SymTensor3D::idx_xy]*S.data[SymTensor3D::idx_xy] - S.data[SymTensor3D::idx_yz]*S.data[SymTensor3D::idx_yz] - S.data[SymTensor3D::idx_xz]*S.data[SymTensor3D::idx_xz];
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

    inline std::ostream &operator<<(std::ostream &os, const SkewSymTensor3D &W)
    {
        os << "SkewSymTensor2[yz:" << W[0] << ", xz:" << W[1] << ", xy:" << W[2] << "]";
        return os;
    }

        inline std::ostream &operator<<(std::ostream &os, const SymTensor3D &T)
    {
        os << "SymTensor2[xx:" << T[0] << ", yy:" << T[1] << ", zz:" << T[2]
           << ", yz:" << T[3] << ", xz:" << T[4] << ", xy:" << T[5] << "]";
        return os;
    }


}