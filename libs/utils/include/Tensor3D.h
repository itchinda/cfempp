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
            return {{{xx(), xy(), xz()},
                     {yx(), yy(), yz()},
                     {zx(), zy(), zz()}}};
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
            return {{{xx(), xy(), xz()},
                     {xy(), yy(), yz()},
                     {xz(), yz(), zz()}}};
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

        [[nodiscard]] inline SymTensor3D operator/(CfemReal scalar) const noexcept
        {
            if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
            {
                CFEM_WARN << "Division by zero in SymTensor3D::operator/";
            }

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

            const CfemReal Axx = xx(), Ayy = yy(), Azz = zz(), Axy = xy(), Axz = xz(), Ayz = yz();
            const CfemReal Bxx = B.xx(), Byy = B.yy(), Bzz = B.zz(), Bxy = B.xy(), Bxz = B.xz(), Byz = B.yz();

            res.xx() = Axx * Bxx + Axy * Bxy + Axz * Bxz;
            res.xy() = Axx * Bxy + Axy * Byy + Axz * Byz;
            res.xz() = Axx * Bxz + Axy * Byz + Axz * Bzz;

            res.yx() = Axy * Bxx + Ayy * Bxy + Ayz * Bxz;
            res.yy() = Axy * Bxy + Ayy * Byy + Ayz * Byz;
            res.yz() = Axy * Bxz + Ayy * Byz + Ayz * Bzz;

            res.zx() = Axz * Bxx + Ayz * Bxy + Azz * Bxz;
            res.zy() = Axz * Bxy + Ayz * Byy + Azz * Byz;
            res.zz() = Axz * Bxz + Ayz * Byz + Azz * Bzz;

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

        inline SymTensor3D &operator/=(CfemReal scalar) noexcept
        {
            if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
            {
                CFEM_WARN << "Division by zero in SymTensor3D::operator/=";
            }

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

            // Ligne 0 de (*this) : [xx, xy, xz]
            const CfemReal m00 = M.xx(), m01 = M.xy(), m02 = M.xz();
            const CfemReal m10 = M.yx(), m11 = M.yy(), m12 = M.yz();
            const CfemReal m20 = M.zx(), m21 = M.zy(), m22 = M.zz();

            T.xx() = data[idx_xx] * m00 + data[idx_xy] * m10 + data[idx_xz] * m20;
            T.xy() = data[idx_xx] * m01 + data[idx_xy] * m11 + data[idx_xz] * m21;
            T.xz() = data[idx_xx] * m02 + data[idx_xy] * m12 + data[idx_xz] * m22;

            // Ligne 1 de (*this) : [xy, yy, yz]
            T.yx() = data[idx_xy] * m00 + data[idx_yy] * m10 + data[idx_yz] * m20;
            T.yy() = data[idx_xy] * m01 + data[idx_yy] * m11 + data[idx_yz] * m21;
            T.yz() = data[idx_xy] * m02 + data[idx_yy] * m12 + data[idx_yz] * m22;

            // Ligne 2 de (*this) : [xz, yz, zz]
            T.zx() = data[idx_xz] * m00 + data[idx_yz] * m10 + data[idx_zz] * m20;
            T.zy() = data[idx_xz] * m01 + data[idx_yz] * m11 + data[idx_zz] * m21;
            T.zz() = data[idx_xz] * m02 + data[idx_yz] * m12 + data[idx_zz] * m22;

            // On calcule res = M^T * T
            // On ne calcule que les 6 composantes symétriques
            res[0] = m00 * T.xx() + m10 * T.yx() + m20 * T.zx(); // xx
            res[1] = m01 * T.xy() + m11 * T.yy() + m21 * T.zy(); // yy
            res[2] = m02 * T.xz() + m12 * T.yz() + m22 * T.zz(); // zz

            res[3] = m01 * T.xz() + m11 * T.yz() + m21 * T.zz(); // yz
            res[4] = m00 * T.xz() + m10 * T.yz() + m20 * T.zz(); // xz
            res[5] = m00 * T.xy() + m10 * T.yy() + m20 * T.zy(); // xy

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
            return data[idx_xx] * T.xx() + 
                data[idx_yy] * T.yy() + 
                data[idx_zz] * T.zz() + 
                data[idx_xy] * (T.xy() + T.yx()) + 
                data[idx_xz] * (T.xz() + T.zx()) + 
                data[idx_yz] * (T.yz() + T.zy());
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
        [[nodiscard]] constexpr FixedMatrix<CfemReal, 6, 6> outer(const SymTensor3D& other) const noexcept
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
    struct AxisSymTensor
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
        constexpr AxisSymTensor() noexcept = default;
        constexpr AxisSymTensor(CfemReal xx, CfemReal yy, CfemReal zz, CfemReal xy) noexcept 
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

        [[nodiscard]] constexpr SymTensor3D toSymTensor3D() const noexcept
        {
            return { xx(), yy(), zz(), 0.0, 0.0, xy() };
        }

        inline operator Tensor3D() const noexcept{
            return Tensor3D{ xx(), xy(), 0.0, 
                             xy(), yy(), 0.0,
                             0.0,  0.0,  zz() };
        }

        inline operator SymTensor3D() const noexcept{
            CFEM_INFO << "Implicit converter called";
            return SymTensor3D{ xx(),  yy(), zz(), 0.0, 0.0, xy() };
        }


        [[nodiscard]] constexpr Mat3x3 toMat3x3() const noexcept{
            return {{{xx(), xy(), 0.0}, 
                     {xy(), yy(), 0.0},
                     {0.0,  0.0,  zz()}}};
        }

        [[nodiscard]] constexpr AxisSymTensor operator-() const noexcept
        {
            return {-data[idx_xx], -data[idx_yy], -data[idx_zz], -data[idx_xy]};
        }

        [[nodiscard]] constexpr AxisSymTensor operator*(CfemReal scalar) const noexcept
        {
            return { xx()*scalar, yy()*scalar, zz()*scalar, xy()*scalar };
        }

        [[nodiscard]] inline AxisSymTensor operator/(CfemReal scalar) const noexcept
        {
            if (std::abs(scalar) < ZERO_TOLERANCE){
                CFEM_WARN << "Division by zero in SymTensor3D::operator/";
            }

            const CfemReal inv = 1.0 / scalar;
            return { xx()*inv, yy()*inv, zz()*inv, xy()*inv };
        }

        [[nodiscard]] constexpr AxisSymTensor operator+(const AxisSymTensor& B) const noexcept
        {
            return { xx() + B.xx(), yy() + B.yy(), zz() + B.zz(), xy() + B.xy() };
        }

        [[nodiscard]] constexpr AxisSymTensor operator-(const AxisSymTensor& B) const noexcept
        {
            return { xx() - B.xx(), yy() - B.yy(), zz() - B.zz(), xy() - B.xy() };
        }

        [[nodiscard]] constexpr Vector3D operator*(const Vector3D &v) const noexcept
        {
            return {
                data[idx_xx] * v.x + data[idx_xy] * v.y, // Row x: [xx, xy, xz]
                data[idx_xy] * v.x + data[idx_yy] * v.y, // Row y: [yx, yy, yz]
                data[idx_zz] * v.z                       // Row z: [zx, zy, zz]
            };
        }

        [[nodiscard]] constexpr Tensor3D operator*(const AxisSymTensor& B) const noexcept
        {
            Tensor3D res;
            res.xx() = xx()*B.xx() + xy()*B.xy();
            res.xy() = xx()*B.xy() + xy()*B.yy();
            res.xz() = 0.0;
            
            res.yx() = xy()*B.xx() + yy()*B.xy();
            res.yy() = xy()*B.xy() + yy()*B.yy();
            res.yz() = 0.0;
            
            res.zx() = 0.0;
            res.zy() = 0.0;
            res.zz() = zz()*B.zz();
            return res;
        }
        
        constexpr AxisSymTensor& operator+=(const AxisSymTensor& S) noexcept
        {
            data[idx_xx] += S.data[idx_xx]; data[idx_yy] += S.data[idx_yy];
            data[idx_zz] += S.data[idx_zz]; data[idx_xy] += S.data[idx_xy];
            return *this;
        }

        constexpr AxisSymTensor& operator-=(const AxisSymTensor& S) noexcept
        {
            data[idx_xx] -= S.data[idx_xx]; data[idx_yy] -= S.data[idx_yy];
            data[idx_zz] -= S.data[idx_zz]; data[idx_xy] -= S.data[idx_xy];
            return *this;
        }

        constexpr AxisSymTensor& operator*=(CfemReal scalar) noexcept
        {
            data[idx_xx] *= scalar; data[idx_yy] *= scalar;
            data[idx_zz] *= scalar; data[idx_xy] *= scalar;
            return *this;
        }

        constexpr AxisSymTensor& operator/=(CfemReal scalar) noexcept
        {
            if (std::abs(scalar) < ZERO_TOLERANCE){
                CFEM_WARN << "Division by zero in SymTensor3D::operator/=";
            }

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

        [[nodiscard]] constexpr AxisSymTensor transpose() const noexcept
        {
            return *this; // Symétrique !
        }

        [[nodiscard]] constexpr AxisSymTensor deviatoric() const noexcept
        {
            const CfemReal mean_tr = trace() / 3.0;
            return { data[idx_xx] - mean_tr, 
                     data[idx_yy] - mean_tr, 
                     data[idx_zz] - mean_tr, 
                     data[idx_xy] };
        }

        [[nodiscard]] constexpr AxisSymTensor inverse() const noexcept
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
        [[nodiscard]] constexpr CfemReal doubleDot(const AxisSymTensor& S) const noexcept
        {
            return data[idx_xx]*S.data[idx_xx] + 
                   data[idx_yy]*S.data[idx_yy] + 
                   data[idx_zz]*S.data[idx_zz] + 
                   2.0 * data[idx_xy]*S.data[idx_xy];
        }

        // Produit doublement contracté (S : T) 
        [[nodiscard]] constexpr CfemReal doubleDot(const Tensor3D& T) const noexcept
        {
            return data[idx_xx]*T.xx() + 
                   data[idx_yy]*T.yy() + 
                   data[idx_zz]*T.zz() + 
                   data[idx_xy]*(T.xy() + T.yx());
        }

        // Produit tensoriel externe S (x) S2 -> Matrice de Voigt 4x4 (souvent K_tangent local)
        [[nodiscard]] constexpr FixedMatrix<CfemReal, 4, 4> outer(const AxisSymTensor& other) const noexcept
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
            return {{{0.0,   xy(), xz()}, 
                     {-xy(), 0.0,  yz()},
                     {-xz(), -yz(), 0.0}}};
        }
        

        inline operator Tensor3D() const noexcept{
            return Tensor3D{ 0.0,   xy(), xz(), 
                            -xy(),  0.0,  yz(),
                            -xz(), -yz(), 0.0 };
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

        [[nodiscard]] inline SkewSymTensor3D operator/(CfemReal scalar) const noexcept
        {
            if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
            {
                CFEM_ERROR << "Division by zero in SkewSymTensor3D::operator/";
            }
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

        [[nodiscard]] constexpr Tensor3D operator*(const SkewSymTensor3D &B) const noexcept
        {
            Tensor3D res;
            const CfemReal Axy = xy(), Axz = xz(), Ayz = yz();
            const CfemReal Bxy = B.xy(), Bxz = B.xz(), Byz = B.yz();

            res.xx() = -Axy * Bxy - Axz * Bxz;
            res.xy() = -Axz * Byz;
            res.xz() =  Axy * Byz;

            res.yx() = -Ayz * Bxz;
            res.yy() = -Axy * Bxy - Ayz * Byz;
            res.yz() = -Axy * Bxz;

            res.zx() =  Ayz * Bxy;
            res.zy() = -Axz * Bxy;
            res.zz() = -Axz * Bxz - Ayz * Byz;

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

        inline SkewSymTensor3D &operator/=(CfemReal scalar) noexcept
        {
            if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
            {
                CFEM_ERROR << "Division by zero in SkewSymTensor3D::operator/=";
            }

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
            const CfemReal m00 = M.xx(), m01 = M.xy(), m02 = M.xz();
            const CfemReal m10 = M.yx(), m11 = M.yy(), m12 = M.yz();
            const CfemReal m20 = M.zx(), m21 = M.zy(), m22 = M.zz();

            // Row 0 of W: [0, xy, xz]
            T.xx() = data[idx_xy] * m10 + data[idx_xz] * m20;
            T.xy() = data[idx_xy] * m11 + data[idx_xz] * m21;
            T.xz() = data[idx_xy] * m12 + data[idx_xz] * m22;

            // Row 1 of W: [-xy, 0, yz]
            T.yx() = -data[idx_xy] * m00 + data[idx_yz] * m20;
            T.yy() = -data[idx_xy] * m01 + data[idx_yz] * m21;
            T.yz() = -data[idx_xy] * m02 + data[idx_yz] * m22;

            // Row 2 of W: [-xz, -yz, 0]
            T.zx() = -data[idx_xz] * m00 - data[idx_yz] * m10;
            T.zy() = -data[idx_xz] * m01 - data[idx_yz] * m11;
            T.zz() = -data[idx_xz] * m02 - data[idx_yz] * m12;

            // 2. res = M^T * T (Only extract the 3 skew components)
            res.data[idx_yz] = m01 * T.xz() + m11 * T.yz() + m21 * T.zz(); // yz
            res.data[idx_xz] = m00 * T.xz() + m10 * T.yz() + m20 * T.zz(); // xz
            res.data[idx_xy] = m00 * T.xy() + m10 * T.yy() + m20 * T.zy(); // xy

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
            return std::sqrt( xy()*xy() + xz()*xz() + yz()*yz() );
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

    [[nodiscard]] constexpr AxisSymTensor operator*(CfemReal scalar, const AxisSymTensor &S) noexcept
    {
        return S * scalar;
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
        const CfemReal Txx = T.xx(), Txy = T.xy(), Txz = T.xz();
        const CfemReal Tyx = T.yx(), Tyy = T.yy(), Tyz = T.yz();
        const CfemReal Tzx = T.zx(), Tzy = T.zy(), Tzz = T.zz();

        const CfemReal Sxx = S.xx(), Syy = S.yy(), Szz = S.zz(), Sxy = S.xy(), Sxz = S.xz(), Syz = S.yz();

        res.xx() = Txx * Sxx + Txy * Sxy + Txz * Sxz;
        res.xy() = Txx * Sxy + Txy * Syy + Txz * Syz;
        res.xz() = Txx * Sxz + Txy * Syz + Txz * Szz;

        res.yx() = Tyx * Sxx + Tyy * Sxy + Tyz * Sxz;
        res.yy() = Tyx * Sxy + Tyy * Syy + Tyz * Syz;
        res.yz() = Tyx * Sxz + Tyy * Syz + Tyz * Szz;

        res.zx() = Tzx * Sxx + Tzy * Sxy + Tzz * Sxz;
        res.zy() = Tzx * Sxy + Tzy * Syy + Tzz * Syz;
        res.zz() = Tzx * Sxz + Tzy * Syz + Tzz * Szz;

        return res;
    }

    /** Tensor3D x AxisSymTensor */
    [[nodiscard]] constexpr Tensor3D operator*(const Tensor3D &T, const AxisSymTensor &S) noexcept
    {
        Tensor3D res;
        const CfemReal Txx = T.xx(), Txy = T.xy(), Txz = T.xz();
        const CfemReal Tyx = T.yx(), Tyy = T.yy(), Tyz = T.yz();
        const CfemReal Tzx = T.zx(), Tzy = T.zy(), Tzz = T.zz();

        const CfemReal Sxx = S.xx(), Syy = S.yy(), Szz = S.zz(), Sxy = S.xy();

        res.xx() = Txx * Sxx + Txy * Sxy;
        res.xy() = Txx * Sxy + Txy * Syy;
        res.xz() = Txz * Szz;
        
        res.yx() = Tyx * Sxx + Tyy * Sxy;
        res.yy() = Tyx * Sxy + Tyy * Syy;
        res.yz() = Tyz * Szz;
        
        res.zx() = Tzx * Sxx + Tzy * Sxy;
        res.zy() = Tzx * Sxy + Tzy * Syy;
        res.zz() = Tzz * Szz;
        return res;
    }

    /** Tensor3D x SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const Tensor3D &T, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res;
        const CfemReal Txx = T.xx(), Txy = T.xy(), Txz = T.xz();
        const CfemReal Tyx = T.yx(), Tyy = T.yy(), Tyz = T.yz();
        const CfemReal Tzx = T.zx(), Tzy = T.zy(), Tzz = T.zz();

        const CfemReal Wxy = W.xy(), Wxz = W.xz(), Wyz = W.yz();

        res.xx() = -Txy * Wxy - Txz * Wxz;
        res.yx() = -Tyy * Wxy - Tyz * Wxz;
        res.zx() = -Tzy * Wxy - Tzz * Wxz;

        res.xy() = Txx * Wxy - Txz * Wyz;
        res.yy() = Tyx * Wxy - Tyz * Wyz;
        res.zy() = Tzx * Wxy - Tzz * Wyz;

        res.xz() = Txx * Wxz + Txy * Wyz;
        res.yz() = Tyx * Wxz + Tyy * Wyz;
        res.zz() = Tzx * Wxz + Tzy * Wyz;

        return res;
    }

    /** AxisSymTensor x Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const AxisSymTensor &S, const Tensor3D &T) noexcept
    {
        Tensor3D res;
        const CfemReal Txx = T.xx(), Txy = T.xy(), Txz = T.xz();
        const CfemReal Tyx = T.yx(), Tyy = T.yy(), Tyz = T.yz();
        const CfemReal Tzx = T.zx(), Tzy = T.zy(), Tzz = T.zz();

        const CfemReal Sxx = S.xx(), Syy = S.yy(), Szz = S.zz(), Sxy = S.xy();

        res.xx() = Sxx * Txx + Sxy * Tyx;
        res.xy() = Sxx * Txy + Sxy * Tyy;
        res.xz() = Sxx * Txz + Sxy * Tyz;
        
        res.yx() = Sxy * Txx + Syy * Tyx;
        res.yy() = Sxy * Txy + Syy * Tyy;
        res.yz() = Sxy * Txz + Syy * Tyz;
        
        res.zx() = Szz * Tzx;
        res.zy() = Szz * Tzy;
        res.zz() = Szz * Tzz;
        return res;
    }

    /** AxisSymTensor x SymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const AxisSymTensor& A, const SymTensor3D& S) noexcept
    {
        Tensor3D R;

        R.xx() = A.xx()*S.xx() + A.xy()*S.xy();
        R.xy() = A.xx()*S.xy() + A.xy()*S.yy();
        R.xz() = A.xx()*S.xz() + A.xy()*S.yz();

        R.yx() = A.xy()*S.xx() + A.yy()*S.xy();
        R.yy() = A.xy()*S.xy() + A.yy()*S.yy();
        R.yz() = A.xy()*S.xz() + A.yy()*S.yz();

        R.zx() = A.zz()*S.xz();
        R.zy() = A.zz()*S.yz();
        R.zz() = A.zz()*S.zz();

        return R;
    }

    /** AxisSymTensor x SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const AxisSymTensor &S, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res;
        const CfemReal Wxy = W.xy(), Wxz = W.xz(), Wyz = W.yz();
        const CfemReal Sxx = S.xx(), Syy = S.yy(), Szz = S.zz(), Sxy = S.xy();

        res.xx() =-Sxy * Wxy;
        res.xy() = Sxx * Wxy;
        res.xz() = Sxx * Wxz + Sxy * Wyz;

        res.yx() =-Syy * Wxy;
        res.yy() = Sxy * Wxy;
        res.yz() = Sxy * Wxz + Syy * Wyz;

        res.zx() =-Szz * Wxz;
        res.zy() =-Szz * Wyz;
        res.zz() = 0.0;

        return res;
    }

    /** SymTensor3D x Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const SymTensor3D &S, const Tensor3D &T) noexcept
    {
        Tensor3D res;
        const CfemReal Txx = T.xx(), Txy = T.xy(), Txz = T.xz();
        const CfemReal Tyx = T.yx(), Tyy = T.yy(), Tyz = T.yz();
        const CfemReal Tzx = T.zx(), Tzy = T.zy(), Tzz = T.zz();

        const CfemReal Sxx = S.xx(), Syy = S.yy(), Szz = S.zz(), Sxy = S.xy(), Sxz = S.xz(), Syz = S.yz();

        res.xx() = Sxx * Txx + Sxy * Tyx + Sxz * Tzx;
        res.xy() = Sxx * Txy + Sxy * Tyy + Sxz * Tzy;
        res.xz() = Sxx * Txz + Sxy * Tyz + Sxz * Tzz;

        res.yx() = Sxy * Txx + Syy * Tyx + Syz * Tzx;
        res.yy() = Sxy * Txy + Syy * Tyy + Syz * Tzy;
        res.yz() = Sxy * Txz + Syy * Tyz + Syz * Tzz;

        res.zx() = Sxz * Txx + Syz * Tyx + Szz * Tzx;
        res.zy() = Sxz * Txy + Syz * Tyy + Szz * Tzy;
        res.zz() = Sxz * Txz + Syz * Tyz + Szz * Tzz;

        return res;
    }

    /** SymTensor3D x SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const SymTensor3D &S, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res;
        const CfemReal Sxx = S.xx(), Syy = S.yy(), Szz = S.zz(), Sxy = S.xy(), Sxz = S.xz(), Syz = S.yz();
        const CfemReal Wxy = W.xy(), Wxz = W.xz(), Wyz = W.yz();

        res.xx() =-Sxy * Wxy - Sxz * Wxz;
        res.xy() = Sxx * Wxy - Sxz * Wyz;
        res.xz() = Sxx * Wxz + Sxy * Wyz;

        res.yx() =-Syy * Wxy - Syz * Wxz;
        res.yy() = Sxy * Wxy - Syz * Wyz;
        res.yz() = Sxy * Wxz + Syy * Wyz;

        res.zx() =-Syz * Wxy - Szz * Wxz;
        res.zy() = Sxz * Wxy - Szz * Wyz;
        res.zz() = Sxz * Wxz + Syz * Wyz;

        return res;
    }

    /** SymTensor3D x AxisSymTensor */
    [[nodiscard]] constexpr Tensor3D operator*(const SymTensor3D &S, const AxisSymTensor& A) noexcept
    {
        Tensor3D R;
        const CfemReal Sxx = S.xx(), Syy = S.yy(), Szz = S.zz(), Sxy = S.xy(), Sxz = S.xz(), Syz = S.yz();
        const CfemReal Axx = A.xx(), Ayy = A.yy(), Azz = A.zz(), Axy = A.xy();

        R.xx() = Sxx * Axx + Sxy * Axy;
        R.xy() = Sxx * Axy + Sxy * Ayy;
        R.xz() = Sxz * Azz;

        R.yx() = Sxy * Axx + Syy * Axy;
        R.yy() = Sxy * Axy + Syy * Ayy;
        R.yz() = Syz * Azz;

        R.zx() = Sxz * Axx + Syz * Axy;
        R.zy() = Sxz * Axy + Syz * Ayy;
        R.zz() = Szz * Azz;

        return R;
    }

    /** SkewSymTensor3D x Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const SkewSymTensor3D &W, const Tensor3D &T) noexcept
    {
        Tensor3D res;
        const CfemReal Txx = T.xx(), Txy = T.xy(), Txz = T.xz();
        const CfemReal Tyx = T.yx(), Tyy = T.yy(), Tyz = T.yz();
        const CfemReal Tzx = T.zx(), Tzy = T.zy(), Tzz = T.zz();

        const CfemReal Wxy = W.xy(), Wxz = W.xz(), Wyz = W.yz();

        res.xx() = Wxy * Tyx + Wxz * Tzx;
        res.xy() = Wxy * Tyy + Wxz * Tzy;
        res.xz() = Wxy * Tyz + Wxz * Tzz;

        res.yx() = -Wxy * Txx + Wyz * Tzx;
        res.yy() = -Wxy * Txy + Wyz * Tzy;
        res.yz() = -Wxy * Txz + Wyz * Tzz;

        res.zx() = -Wxz * Txx - Wyz * Tyx;
        res.zy() = -Wxz * Txy - Wyz * Tyy;
        res.zz() = -Wxz * Txz - Wyz * Tyz;

        return res;
    }

    /** SkewSymTensor3D x SymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator*(const SkewSymTensor3D &W, const SymTensor3D &S) noexcept
    {
        Tensor3D res;
        const CfemReal Sxx = S.xx(), Syy = S.yy(), Szz = S.zz(), Sxy = S.xy(), Sxz = S.xz(), Syz = S.yz();
        const CfemReal Wxy = W.xy(), Wxz = W.xz(), Wyz = W.yz();

        res.xx() = Wxy * Sxy + Wxz * Sxz;
        res.xy() = Wxy * Syy + Wxz * Syz;
        res.xz() = Wxy * Syz + Wxz * Szz;

        res.yx() = -Wxy * Sxx + Wyz * Sxz;
        res.yy() = -Wxy * Sxy + Wyz * Syz;
        res.yz() = -Wxy * Sxz + Wyz * Szz;

        res.zx() = -Wxz * Sxx - Wyz * Sxy;
        res.zy() = -Wxz * Sxy - Wyz * Syy;
        res.zz() = -Wxz * Sxz - Wyz * Syz;

        return res;
    }

    /** SkewSymTensor3D x AxisSymTensor */
    [[nodiscard]] constexpr Tensor3D operator*(const SkewSymTensor3D &W, const AxisSymTensor &S) noexcept
    {
        Tensor3D res;
        const CfemReal Sxx = S.xx(), Syy = S.yy(), Szz = S.zz(), Sxy = S.xy();
        const CfemReal Wxy = W.xy(), Wxz = W.xz(), Wyz = W.yz();

        res.xx() =  Wxy * Sxy;             res.xy() = Wxy * Syy;             res.xz() = Wxz * Szz;
        res.yx() = -Wxy * Sxx;             res.yy() =-Wxy * Sxy;             res.yz() = Wyz * Szz;
        res.zx() = -Wxz * Sxx - Wyz * Sxy; res.zy() =-Wxz * Sxy - Wyz * Syy; res.zz() = 0.0;

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
        res.xx() += S.xx(); res.xy() += S.xy(); res.xz() += S.xz();
        res.yx() += S.xy(); res.yy() += S.yy(); res.yz() += S.yz();
        res.zx() += S.xz(); res.zy() += S.yz(); res.zz() += S.zz();
        return res;
    }

    /** Tensor3D + AxisSymTensor */
    [[nodiscard]] constexpr Tensor3D operator+(const Tensor3D &T, const AxisSymTensor &A) noexcept
    {
        Tensor3D res = T;
        res.xx() += A.xx(); res.xy() += A.xy();
        res.yx() += A.xy(); res.yy() += A.yy(); 
        res.zz() += A.zz();
        return res;
    }

    /** Tensor3D + SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const Tensor3D &T, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res = T;
        res.xy() += W.xy();  res.xz() += W.xz();
        res.yx() -= W.xy();  res.yz() += W.yz();
        res.zx() -= W.xz();  res.zy() -= W.yz();
        return res;
    }

    /** SymTensor3D + Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const SymTensor3D &S, const Tensor3D &T) noexcept
    {
        return T + S;
    }

    /** SymTensor3D + AxisSymTensor */
    [[nodiscard]] constexpr SymTensor3D operator+(const SymTensor3D &S, const AxisSymTensor &A) noexcept
    {
        SymTensor3D res = S;
        res.xx() += A.xx(); res.yy() += A.yy(); res.zz() += A.zz(); res.xy() += A.xy();
        return res;
    }

    /** SymTensor3D + SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const SymTensor3D &S, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res;
        res.xx() = S.xx();          res.xy() = S.xy() + W.xy(); res.xz() = S.xz() + W.xz();
        res.yx() = S.xy() - W.xy(); res.yy() = S.yy();          res.yz() = S.yz() + W.yz();
        res.zx() = S.xz() - W.xz(); res.zy() = S.yz() - W.yz(); res.zz() = S.zz();
        return res;
    }

    /** AxisSymTensor + Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const AxisSymTensor &A, const Tensor3D &T) noexcept
    {
        return T + A;
    }

    /** AxisSymTensor + SymTensor3D */
    [[nodiscard]] constexpr SymTensor3D operator+(const AxisSymTensor &A, const SymTensor3D &S) noexcept
    {
        return S + A;
    }

    /** AxisSymTensor + SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const AxisSymTensor &A, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res;
        res.xx() = A.xx();          res.xy() = A.xy() + W.xy(); res.xz() =          W.xz();
        res.yx() = A.xy() - W.xy(); res.yy() = A.yy();          res.yz() =          W.yz();
        res.zx() =        - W.xz(); res.zy() =        - W.yz(); res.zz() = A.zz();
        return res;
    }

    /** SkewSymTensor3D + Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const SkewSymTensor3D &S, const Tensor3D &T) noexcept
    {
        return T + S;
    }

    /** SkewSymTensor3D + SymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator+(const SkewSymTensor3D &W, const SymTensor3D &S) noexcept
    {
        return S + W; 
    }

    /** SkewSymTensor3D + AxisSymTensor */
    [[nodiscard]] constexpr Tensor3D operator+(const SkewSymTensor3D &W, const AxisSymTensor &A) noexcept
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
        res.xx() -= S.xx(); res.xy() -= S.xy(); res.xz() -= S.xz();
        res.yx() -= S.xy(); res.yy() -= S.yy(); res.yz() -= S.yz();
        res.zx() -= S.xz(); res.zy() -= S.yz(); res.zz() -= S.zz();
        return res;
    }

    /** Tensor3D - AxisSymTensor */
    [[nodiscard]] constexpr Tensor3D operator-(const Tensor3D &T, const AxisSymTensor &S) noexcept
    {
        Tensor3D res = T;
        res.xx() -= S.xx(); res.xy() -= S.xy();
        res.yx() -= S.xy(); res.yy() -= S.yy();
        res.zz() -= S.zz();
        return res;
    }

    /** Tensor3D - SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const Tensor3D &T, const SkewSymTensor3D &S) noexcept
    {
        Tensor3D res = T;
        res.xy() -= S.xy(); res.xz() -= S.xz();
        res.yx() += S.xy(); res.yz() -= S.yz();
        res.zx() += S.xz(); res.zy() += S.yz();

        return res;
    }

    /** SymTensor3D - Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const SymTensor3D &S, const Tensor3D &T) noexcept
    {
        Tensor3D res;
        res.xx() = S.xx() - T.xx(); res.xy() = S.xy() - T.xy(); res.xz() = S.xz() - T.xz();
        res.yx() = S.xy() - T.yx(); res.yy() = S.yy() - T.yy(); res.yz() = S.yz() - T.yz();
        res.zx() = S.xz() - T.zx(); res.zy() = S.yz() - T.zy(); res.zz() = S.zz() - T.zz();
        return res;
    }

    /** SymTensor3D - SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const SymTensor3D &S, const SkewSymTensor3D &W) noexcept
    {
        Tensor3D res;
        res.xx() = S.xx();          res.xy() = S.xy() - W.xy(); res.xz() = S.xz() - W.xz();
        res.yx() = S.xy() + W.xy(); res.yy() = S.yy();          res.yz() = S.yz() - W.yz();
        res.zx() = S.xz() + W.xz(); res.zy() = S.yz() + W.yz(); res.zz() = S.zz();
        return res;
    }

    /** SymTensor3D - AxisSymTensor */
    [[nodiscard]] constexpr SymTensor3D operator-(const SymTensor3D &S, const AxisSymTensor &A) noexcept
    {
        return { S.xx() - A.xx(), S.yy() - A.yy(), S.zz() - A.zz(), S.yz(), S.xz(), S.xy() - A.xy() };
    }

    /** AxisSymTensor - Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const AxisSymTensor &S, const Tensor3D &T) noexcept
    {
        Tensor3D res;
        res.xx() = S.xx() - T.xx(); res.xy() = S.xy() - T.xy(); res.xz() =        - T.xz();
        res.yx() = S.xy() - T.yx(); res.yy() = S.yy() - T.yy(); res.yz() =        - T.yz();
        res.zx() =        - T.zx(); res.zy() =        - T.zy(); res.zz() = S.zz() - T.zz();
        return res;
    }

    /** AxisSymTensor - SymTensor3D */
    [[nodiscard]] constexpr SymTensor3D operator-(const AxisSymTensor &A, const SymTensor3D &S) noexcept
    {
        return { A.xx() - S.xx(), A.yy() - S.yy(), A.zz() - S.zz(), -S.yz(), -S.xz(), A.xy() - S.xy() };
    }

    /** AxisSymTensor - SkewSymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const AxisSymTensor &A, const SkewSymTensor3D &W) noexcept
    {
        return { A.xx(),          A.xy() - W.xy(),         -W.xz(),
                 A.xy() + W.xy(), A.yy(),                  -W.yz(),
                          W.xz(),          W.yz(),  A.zz()          };
    }

    /** SkewSymTensor3D - Tensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const SkewSymTensor3D &S, const Tensor3D &T) noexcept
    {
        Tensor3D res;
        res.xx() =         - T.xx(); res.xy() = S.xy() - T.xy(); res.xz() = S.xz() - T.xz();
        res.yx() = -S.xy() - T.yx(); res.yy() =        - T.yy(); res.yz() = S.yz() - T.yz();
        res.zx() = -S.xz() - T.zx(); res.zy() =-S.yz() - T.zy(); res.zz() =        - T.zz();
        return res;
    }

    /** SkewSymTensor3D - SymTensor3D */
    [[nodiscard]] constexpr Tensor3D operator-(const SkewSymTensor3D &W, const SymTensor3D &S) noexcept
    {
        Tensor3D res;
        res.xx() = -S.xx();          res.xy() =  W.xy() - S.xy(); res.xz() =  W.xz() - S.xz();
        res.yx() = -W.xy() - S.xy(); res.yy() = -S.yy();          res.yz() =  W.yz() - S.yz();
        res.zx() = -W.xz() - S.xz(); res.zy() = -W.yz() - S.yz(); res.zz() = -S.zz();
        return res;
    }

    /** SkewSymTensor3D - AxisSymTensor */
    [[nodiscard]] constexpr Tensor3D operator-(const SkewSymTensor3D &W, const AxisSymTensor &A) noexcept
    {
        Tensor3D res;
        res.xx() =         - A.xx(); res.xy() =  W.xy() - A.xy(); res.xz() =  W.xz();
        res.yx() = -W.xy() - A.xy(); res.yy() =         - A.yy(); res.yz() =  W.yz();
        res.zx() = -W.xz();          res.zy() = -W.yz();          res.zz() =         - A.zz();
        return res;
    }
    /** @} */

    [[nodiscard]] constexpr CfemReal det(const SymTensor3D &S) noexcept
    {
        return S.xx() * (S.yy() * S.zz() - S.yz() * S.yz())
             - S.xy() * (S.xy() * S.zz() - S.yz() * S.xz())
             + S.xz() * (S.xy() * S.yz() - S.yy() * S.xz());
    }

    [[nodiscard]] constexpr CfemReal det(const SkewSymTensor3D &S) noexcept
    {
        return 0.0;
    }

    // [[nodiscard]] inline SymTensor3D inverse(const SymTensor3D &S, CfemReal determinant) noexcept
    // {
    //     CFEM_ASSERT(std::abs(determinant) < cfem::ZERO_TOLERANCE), "Singular 3x3 matrix/Jacobian detected.";
    //     // if (std::abs(determinant) < cfem::ZERO_TOLERANCE)
    //     //     CFEM_ERROR << "Singular 3x3 matrix/Jacobian detected.";

    //     // En pratique, tu peux ajouter ici un assert(std::abs(det) > 1e-14);
    //     const CfemReal inv_det = 1.0 / determinant;
        
    //     SymTensor3D res;
    //     res.xx() = (S.yy() * S.zz() - S.yz() * S.yz()) * inv_det;
    //     res.yy() = (S.xx() * S.zz() - S.xz() * S.xz()) * inv_det;
    //     res.zz() = (S.xx() * S.yy() - S.xy() * S.xy()) * inv_det;
    //     res.xy() = (S.xz() * S.yz() - S.xy() * S.zz()) * inv_det;
    //     res.xz() = (S.xy() * S.yz() - S.xx() * S.xz()) * inv_det;
    //     res.yz() = (S.xy() * S.xz() - S.yy() * S.xz()) * inv_det;
    //     return res;
    // }

    [[nodiscard]] inline SymTensor3D inverse(const SymTensor3D &S, CfemReal determinant)
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

    /**
     * @brief Tries to inverse a skew-symmetric second order tensor (This operation is impossible).
     * @throw std::runtime_error.
     */
    [[nodiscard]] inline SkewSymTensor3D inverse(const SkewSymTensor3D &S) 
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
    [[nodiscard]] inline FixedVector<CfemReal, 3> eigenvalues(const SymTensor3D& S) noexcept
    {
        const CfemReal I1 = S.xx() + S.yy() + S.zz();
        const CfemReal I2 = S.xx()*S.yy() + S.yy()*S.zz() + S.xx()*S.zz() 
                          - S.xy()*S.xy() - S.yz()*S.yz() - S.xz()*S.xz();
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