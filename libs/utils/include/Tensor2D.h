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
 * @file Tensor2D.h
 * @brief Physical tensor definitions for CFEM.
 */

#pragma once

#include "FixedMatrix.h"
#include "Vector2D.h"

namespace cfem
{

    /**
     * @brief General 2x2 tensor (4 components).
     * Inherits from FixedMatrix to add specific physical indices.
     */
    struct Tensor2D : public FixedMatrix<CfemReal, 2, 2>
    {
        using FixedMatrix<CfemReal, 2, 2>::FixedMatrix;

        constexpr Tensor2D(const FixedMatrix<CfemReal, 2, 2>& base) noexcept 
            : FixedMatrix<CfemReal, 2, 2>(base) {}

        /**
         * @brief Storage indices for symmetric 2x2 tensors.
         */
        static constexpr CfemInt idx_xx = 0;
        static constexpr CfemInt idx_xy = 1;
        static constexpr CfemInt idx_yx = 2;
        static constexpr CfemInt idx_yy = 3;
    };

    /**
     * @brief Identity tensor constant.
     */
    [[nodiscard]] constexpr Tensor2D IdentityTensor2D() noexcept
    {
        Tensor2D I;
        I.setIdentity();
        return I;
    }

    /**
     * @brief Stores a symmetric 2x2 tensor using only 3 values.
     * Storage order (Voigt): [0]=xx, [1]=yy, [2]=xy
     */
    struct SymTensor2D
    {
        CfemReal data[3]{};


        /**
         * @brief Storage indices for symmetric 2x2 tensors.
         * 
         * @{
         */
        static constexpr CfemInt idx_xx = 0;
        static constexpr CfemInt idx_yy = 1;
        static constexpr CfemInt idx_xy = 2;
        /** @} */

    private:
        static constexpr CfemInt map[2][2] = { {0, 2}, {2, 1} };

    public:
        constexpr SymTensor2D() noexcept = default;

        constexpr SymTensor2D(CfemReal xx, CfemReal yy, CfemReal xy) noexcept 
        : data{xx, yy, xy} 
        {}

        /** @brief Access by component index (0 to 2). */
        constexpr CfemReal &operator[](CfemInt i)
        {
            CFEM_ASSERT(i >= 0 && i < 3);
            return data[i];
        }

        [[nodiscard]] constexpr const CfemReal &operator[](CfemInt i) const
        {
            CFEM_ASSERT(i >= 0 && i < 3);
            return data[i];
        }

        /**
         * @brief Map 2x2 indices to the symmetric storage.
         *
         * @param i Row index
         * @param j Column index
         * @return CfemReal
         */
        [[nodiscard]] constexpr CfemReal operator()(CfemInt i, CfemInt j) const
        {
            CFEM_ASSERT(i >= 0 && i < 2 && j >= 0 && j < 2);
            return data[map[i][j]];
        }

        /**
         * @brief Negate operator
         *
         * @return SymTensor2D
         */
        [[nodiscard]] constexpr SymTensor2D operator-() const noexcept
        {
            return {-data[idx_xx], -data[idx_yy], -data[idx_xy]};
        }

        [[nodiscard]] constexpr CfemReal xx() const noexcept { return data[idx_xx]; }
        [[nodiscard]] constexpr CfemReal yy() const noexcept { return data[idx_yy]; }
        [[nodiscard]] constexpr CfemReal xy() const noexcept { return data[idx_xy]; }

        [[nodiscard]] constexpr CfemReal &xx() noexcept { return data[idx_xx]; }
        [[nodiscard]] constexpr CfemReal &yy() noexcept { return data[idx_yy]; }
        [[nodiscard]] constexpr CfemReal &xy() noexcept { return data[idx_xy]; }

        // Basic Math Ops

        [[nodiscard]] static constexpr SymTensor2D Identity() noexcept
        {
            return {1, 1, 0};
        }

        [[nodiscard]] static constexpr SymTensor2D Zero() noexcept
        {
            return {0, 0, 0};
        }

        constexpr void setZero() noexcept
        {
            for (CfemInt i = 0; i < 3; ++i)
                data[i] = 0;
        }

        constexpr void setTo(CfemReal value) noexcept
        {
            for (CfemInt i = 0; i < 3; ++i)
                data[i] = value;
        }

        constexpr void toMatrix(CfemReal M[2][2]) const noexcept
        {
            M[0][0] = data[idx_xx];
            M[0][1] = data[idx_xy];
            M[1][0] = data[idx_xy];
            M[1][1] = data[idx_yy];
        }

        [[nodiscard]] constexpr SymTensor2D operator+(const SymTensor2D &other) const noexcept
        {
            return { data[idx_xx] + other.data[idx_xx], 
                     data[idx_yy] + other.data[idx_yy],
                     data[idx_xy] + other.data[idx_xy] };
        }

        [[nodiscard]] constexpr SymTensor2D operator-(const SymTensor2D &other) const noexcept
        {
            return {data[idx_xx] - other.data[idx_xx],
                    data[idx_yy] - other.data[idx_yy],
                    data[idx_xy] - other.data[idx_xy]};
        }

        [[nodiscard]] constexpr SymTensor2D operator*(CfemReal scalar) const noexcept
        {
            return {data[idx_xx] * scalar, 
                    data[idx_yy] * scalar, 
                    data[idx_xy] * scalar};
        }

        [[nodiscard]] inline SymTensor2D operator/(CfemReal scalar) const noexcept
        {
            if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
            {
                CFEM_ERROR << "Division by zero in SymTensor2D::operator/";
            }
            const CfemReal inv = 1.0 / scalar;
            return {data[idx_xx] * inv,
                    data[idx_yy] * inv,
                    data[idx_xy] * inv };
        }

        [[nodiscard]] constexpr Vector2D operator*(const Vector2D &v) const noexcept
        {
            return {
                data[idx_xx] * v.x + data[idx_xy] * v.y, // Row x: [xx, xy, xz]
                data[idx_xy] * v.x + data[idx_yy] * v.y // Row y: [yx, yy, yz]
            };
        }

        constexpr SymTensor2D &operator+=(const SymTensor2D &other) noexcept
        {
            data[idx_xx] += other.data[idx_xx];
            data[idx_yy] += other.data[idx_yy];
            data[idx_xy] += other.data[idx_xy];
            return *this;
        }

        constexpr SymTensor2D &operator-=(const SymTensor2D &other) noexcept
        {
            data[idx_xx] -= other.data[idx_xx];
            data[idx_yy] -= other.data[idx_yy];
            data[idx_xy] -= other.data[idx_xy];
            return *this;
        }

        constexpr SymTensor2D &operator*=(CfemReal scalar) noexcept
        {
            data[idx_xx] *= scalar;
            data[idx_yy] *= scalar;
            data[idx_xy] *= scalar;
            return *this;
        }

        inline SymTensor2D &operator/=(CfemReal scalar) noexcept
        {
            if (std::abs(scalar) < cfem::ZERO_TOLERANCE)
            {
                CFEM_ERROR << "Division by zero in SymTensor2D::operator/=";
            }

            CfemReal inv = 1.0 / scalar;
            data[idx_xx] *= inv;
            data[idx_yy] *= inv;
            data[idx_xy] *= inv;
            return *this;
        }

        [[nodiscard]] constexpr bool operator==(const SymTensor2D &other) const noexcept
        {
            for (CfemInt i = 0; i < 3; ++i)
                if (data[i] != other[i])
                    return false;
            return true;
        }

        [[nodiscard]] constexpr bool operator!=(const SymTensor2D &other) const noexcept
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
        [[nodiscard]] constexpr SymTensor2D quadform(const Tensor2D &M) const noexcept
        {
            SymTensor2D res;

            // On calcule d'abord le produit intermédiaire T = (*this) * M
            // Rappel stockage : [0]=xx, [1]=yy, [2]=xy
            Tensor2D T;

            // Ligne 0 de (*this) : [xx, xy]
            const CfemReal m00 = M.data[Tensor2D::idx_xx], m01 = M.data[Tensor2D::idx_xy];
            const CfemReal m10 = M.data[Tensor2D::idx_yx], m11 = M.data[Tensor2D::idx_yy];

            T.data[Tensor2D::idx_xx] = data[idx_xx] * m00 + data[idx_xy] * m10;
            T.data[Tensor2D::idx_xy] = data[idx_xx] * m01 + data[idx_xy] * m11;

            // Ligne 1 de (*this) : [xy, yy]
            T.data[Tensor2D::idx_yx] = data[idx_xy] * m00 + data[idx_yy] * m10;
            T.data[Tensor2D::idx_yy] = data[idx_xy] * m01 + data[idx_yy] * m11;

            // Compute res = M^T * T
            res[0] = m00 * T.data[Tensor2D::idx_xx] + m10 * T.data[Tensor2D::idx_yx]; // xx
            res[1] = m01 * T.data[Tensor2D::idx_xy] + m11 * T.data[Tensor2D::idx_yy]; // yy
            res[2] = m00 * T.data[Tensor2D::idx_xy] + m10 * T.data[Tensor2D::idx_yy]; // xy

            return res;
        }

        // Physical Invariants

        /** @brief Trace of the tensor (sum of diagonals). */
        [[nodiscard]] constexpr CfemReal trace() const noexcept
        {
            return data[idx_xx] + data[idx_yy];
        }

        [[nodiscard]] constexpr CfemReal determinant() const noexcept
        {
            return data[idx_xx] * data[idx_yy] - data[idx_xy] * data[idx_xy];
        }

        [[nodiscard]] inline SymTensor2D inverse() const noexcept
        {
            const CfemReal det = determinant();
            if (std::abs(det) < cfem::ZERO_TOLERANCE)
                CFEM_ERROR << "Singular 2x2 matrix/Jacobian detected.";

            const CfemReal inv_det = 1.0 / det;
            
            SymTensor2D res;
            res.data[idx_xx] =  data[idx_yy] * inv_det;
            res.data[idx_yy] =  data[idx_xx] * inv_det;
            res.data[idx_xy] = -data[idx_xy] * inv_det;
            return res;
        }

        [[nodiscard]] constexpr SymTensor2D transpose() const noexcept
        {
            return *this;
        }

        [[nodiscard]] constexpr CfemReal doubleDot(const Tensor2D& T) const noexcept
        {
            // S : T = S_xx*T_xx + S_yy*T_yy + S_xy*(T_xy + T_yx)
            return data[idx_xx] * T.data[Tensor2D::idx_xx] + 
                   data[idx_yy] * T.data[Tensor2D::idx_yy] + 
                   data[idx_xy] *(T.data[Tensor2D::idx_xy] + T.data[Tensor2D::idx_yx]);
        }

        [[nodiscard]] constexpr SymTensor2D deviatoric() const noexcept
        {
            // En 2D (déformations planes ou contraintes planes), 
            // la pression moyenne s'écrit souvent p = (T_xx + T_yy) / 2
            // Attention : selon ton hypothèse mécanique (ex: contraintes planes vs def planes), 
            // le terme 1/2 ou 1/3 peut varier si on considère la composante hors-plan !
            const CfemReal mean = 0.5 * (data[idx_xx] + data[idx_yy]);
            SymTensor2D res = *this;
            res.data[idx_xx] -= mean;
            res.data[idx_yy] -= mean;
            return res;
        }

        // Calcul analytique exact des valeurs propres via le discriminant
        [[nodiscard]] inline FixedVector<CfemReal, 2> eigenvalues() const noexcept
        {
            const CfemReal tr = data[idx_xx] + data[idx_yy];
            const CfemReal gap = data[idx_xx] - data[idx_yy];
            const CfemReal desc = std::sqrt(gap * gap + 4.0 * data[idx_xy] * data[idx_xy]);

            FixedVector<CfemReal, 2> eig;
            eig[0] = 0.5 * (tr + desc);
            eig[1] = 0.5 * (tr - desc);
            return eig;
        }

        
    };

    // Global Ops

    [[nodiscard]] constexpr SymTensor2D operator*(CfemReal scalar, const SymTensor2D &S) noexcept
    {
        return S * scalar;
    }

    [[nodiscard]] constexpr Tensor2D operator+(const Tensor2D &T, const SymTensor2D &S) noexcept
    {
        Tensor2D res = T;
        res.data[Tensor2D::idx_xx] += S.data[SymTensor2D::idx_xx];
        res.data[Tensor2D::idx_yy] += S.data[SymTensor2D::idx_yy];
        res.data[Tensor2D::idx_xy] += S.data[SymTensor2D::idx_xy];
        res.data[Tensor2D::idx_yx] += S.data[SymTensor2D::idx_xy];
        return res;
    }

    [[nodiscard]] constexpr Tensor2D operator+(const SymTensor2D &S, const Tensor2D &T) noexcept
    {
        return T + S;
    }

    [[nodiscard]] constexpr Tensor2D operator-(const Tensor2D &T, const SymTensor2D &S) noexcept
    {
        Tensor2D res = T;
        res.data[Tensor2D::idx_xx] -= S.data[SymTensor2D::idx_xx];
        res.data[Tensor2D::idx_yy] -= S.data[SymTensor2D::idx_yy];
        res.data[Tensor2D::idx_xy] -= S.data[SymTensor2D::idx_xy];
        res.data[Tensor2D::idx_yx] -= S.data[SymTensor2D::idx_xy];
        return res;
    }

    [[nodiscard]] constexpr Tensor2D operator-(const SymTensor2D &S, const Tensor2D &T) noexcept
    {
        Tensor2D res;
        res.data[Tensor2D::idx_xx] = S.data[SymTensor2D::idx_xx] - T.data[Tensor2D::idx_xx];
        res.data[Tensor2D::idx_yy] = S.data[SymTensor2D::idx_yy] - T.data[Tensor2D::idx_yy];
        res.data[Tensor2D::idx_xy] = S.data[SymTensor2D::idx_xy] - T.data[Tensor2D::idx_xy];
        res.data[Tensor2D::idx_yx] = S.data[SymTensor2D::idx_xy] - T.data[Tensor2D::idx_yx];
        return res;
    }

    [[nodiscard]] constexpr Tensor2D operator*(const Tensor2D &T, const SymTensor2D &S) noexcept
    {
        Tensor2D res;
        res.data[Tensor2D::idx_xx] = T.data[Tensor2D::idx_xx]*S.data[SymTensor2D::idx_xx] + T.data[Tensor2D::idx_xy]*S.data[SymTensor2D::idx_xy];
        res.data[Tensor2D::idx_xy] = T.data[Tensor2D::idx_xx]*S.data[SymTensor2D::idx_xy] + T.data[Tensor2D::idx_xy]*S.data[SymTensor2D::idx_yy];
        
        res.data[Tensor2D::idx_yx] = T.data[Tensor2D::idx_yx]*S.data[SymTensor2D::idx_xx] + T.data[Tensor2D::idx_yy]*S.data[SymTensor2D::idx_xy];
        res.data[Tensor2D::idx_yy] = T.data[Tensor2D::idx_yx]*S.data[SymTensor2D::idx_xy] + T.data[Tensor2D::idx_yy]*S.data[SymTensor2D::idx_yy];
        return res;
    }

    [[nodiscard]] constexpr Tensor2D operator*(const SymTensor2D &S, const Tensor2D &T) noexcept
    {
        Tensor2D res;
        res.data[Tensor2D::idx_xx] = S.data[SymTensor2D::idx_xx]*T.data[Tensor2D::idx_xx] + S.data[SymTensor2D::idx_xy]*T.data[Tensor2D::idx_yx];
        res.data[Tensor2D::idx_xy] = S.data[SymTensor2D::idx_xx]*T.data[Tensor2D::idx_xy] + S.data[SymTensor2D::idx_xy]*T.data[Tensor2D::idx_yy];
        
        res.data[Tensor2D::idx_yx] = S.data[SymTensor2D::idx_xy]*T.data[Tensor2D::idx_xx] + S.data[SymTensor2D::idx_yy]*T.data[Tensor2D::idx_yx];
        res.data[Tensor2D::idx_yy] = S.data[SymTensor2D::idx_xy]*T.data[Tensor2D::idx_xy] + S.data[SymTensor2D::idx_yy]*T.data[Tensor2D::idx_yy];
        return res;
    }

    [[nodiscard]] constexpr CfemReal det(const SymTensor2D &S) noexcept
    {
        return S.determinant();
    }

    [[nodiscard]] inline SymTensor2D inverse(const SymTensor2D &S) noexcept
    {
        return S.inverse();
    }

    /**
     * @brief Computes the eigenvalues of a symmetric 2x2 tensor analytically.
     *
     * @param S Symmetric tensor.
     * @return Eigenvalues sorted in descending order.
     *
     * @details
     * Uses the closed-form Cardano method for symmetric tensors.
     */
    [[nodiscard]] inline FixedVector<CfemReal, 2> eigenvalues( const SymTensor2D &S ) noexcept
    {
        const CfemReal tr = S.trace();
        const CfemReal determinant = det(S);
        // Discriminant : (xx-yy)^2 + 4*xy^2
        const CfemReal gap = S.data[SymTensor2D::idx_xx] - S.data[SymTensor2D::idx_yy];
        const CfemReal desc = std::sqrt(gap * gap + 4.0 * S.data[SymTensor2D::idx_xy] * S.data[SymTensor2D::idx_xy]);

        FixedVector<CfemReal, 2> eig;
        eig[0] = 0.5 * (tr + desc);
        eig[1] = 0.5 * (tr - desc);
        return eig;
    }

    inline std::ostream &operator<<(std::ostream &os, const SymTensor2D &T)
    {
        os << "SymTensor2[xx:" << T.data[SymTensor2D::idx_xx] << ", yy:" << T.data[SymTensor2D::idx_yy] << ", xy:" << T.data[SymTensor2D::idx_xy] << "]";
        return os;
    }

} // namespace
