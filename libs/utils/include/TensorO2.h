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
 * @brief General second-order tensor (9 components).
 * Used for non-symmetric quantities like the Deformation Gradient (F)
 * or the Velocity Gradient (L).
 */
using TensorO2 = FixedMatrix<CfemReal, 3, 3>;

/**
 * @brief Identity tensor constant.
 */
inline TensorO2 IdentityTensor() {
    TensorO2 I;
    I.setIdentity();
    return I;
}

/**
 * @brief Stores a symmetric 3x3 tensor using only 6 values.
 * Storage order (Voigt): [0]=xx, [1]=yy, [2]=zz, [3]=yz, [4]=xz, [5]=xy
 */
struct SymTensorO2 {
    CfemReal data[6] = {0};

    SymTensorO2() = default; 

    SymTensorO2(std::initializer_list<CfemReal> list){
        CFEM_ASSERT(list.size() == 6);
        int i = 0;
        for (auto val : list) data[i++] = val;
    }
    
    /** @brief Access by component index (0 to 5). */
    inline CfemReal& operator[](int i) {
        CFEM_ASSERT(i >= 0 && i < 6);
        return data[i];
    }

    inline const CfemReal& operator[](int i) const {
        CFEM_ASSERT(i >= 0 && i < 6);
        return data[i];
    }

    /** @brief Map 3x3 indices to the symmetric storage. */
    inline CfemReal operator()(int i, int j) const {
        CFEM_ASSERT(i >= 0 && i < 3 && j >= 0 && j < 3);
        if (i == j) return data[i]; // xx, yy, zz
        if ((i==1 && j==2) || (i==2 && j==1)) return data[3]; // yz
        if ((i==0 && j==2) || (i==2 && j==0)) return data[4]; // xz
        if ((i==0 && j==1) || (i==1 && j==0)) return data[5]; // xy
        return 0.0;
    }

    SymTensorO2 operator-() const {
        return {-data[0], -data[1], -data[2],
                -data[3], -data[4], -data[5]};
    }

    inline CfemReal xx() const{ return data[0];}
    inline CfemReal yy() const{ return data[1];}
    inline CfemReal zz() const{ return data[2];}
    inline CfemReal yz() const{ return data[3];}
    inline CfemReal xz() const{ return data[4];}
    inline CfemReal xy() const{ return data[5];}

    inline CfemReal& xx() { return data[0];}
    inline CfemReal& yy() { return data[1];}
    inline CfemReal& zz() { return data[2];}
    inline CfemReal& yz() { return data[3];}
    inline CfemReal& xz() { return data[4];}
    inline CfemReal& xy() { return data[5];}

    // --- Basic Math Ops ---

    static SymTensorO2 Identity() {
        return {1,1,1,0,0,0};
    }

    static SymTensorO2 Zero() {
        return {0,0,0,0,0,0};
    }

    void setZero(){ for (CfemInt i=0; i<6; ++i) data[i]=0; }
    void setTo(CfemReal value){ for (CfemInt i=0; i<6; ++i) data[i]=value; }

    inline void toMatrix(CfemReal M[3][3]) const {
        M[0][0] = data[0]; M[0][1] = data[5]; M[0][2] = data[4];
        M[1][0] = data[5]; M[1][1] = data[1]; M[1][2] = data[3];
        M[2][0] = data[4]; M[2][1] = data[3]; M[2][2] = data[2];
    }

    SymTensorO2 operator+(const SymTensorO2& other) const {
        return {data[0]+other[0], data[1]+other[1], data[2]+other[2],
                data[3]+other[3], data[4]+other[4], data[5]+other[5]};
    }

    SymTensorO2 operator-(const SymTensorO2& other) const {
        return {data[0]-other[0], data[1]-other[1], data[2]-other[2],
                data[3]-other[3], data[4]-other[4], data[5]-other[5]};
    }

    SymTensorO2 operator*(CfemReal scalar) const {
        return {data[0]*scalar, data[1]*scalar, data[2]*scalar,
                data[3]*scalar, data[4]*scalar, data[5]*scalar};
    }

    SymTensorO2 operator/(CfemReal scalar) const {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        return {data[0]/scalar, data[1]/scalar, data[2]/scalar,
                data[3]/scalar, data[4]/scalar, data[5]/scalar};
    }

    inline Vector3D operator*(const Vector3D& v) const {
        return {
            data[0]*v.x + data[5]*v.y + data[4]*v.z, // Row x: [xx, xy, xz]
            data[5]*v.x + data[1]*v.y + data[3]*v.z, // Row y: [yx, yy, yz]
            data[4]*v.x + data[3]*v.y + data[2]*v.z  // Row z: [zx, zy, zz]
        };
    }

    SymTensorO2& operator+=(const SymTensorO2& other) {
        data[0] += other[0]; data[1] += other[1]; data[2] += other[2];
        data[3] += other[3]; data[4] += other[4]; data[5] += other[5];
        return *this;
    }

    SymTensorO2& operator-=(const SymTensorO2& other) {
        data[0] -= other[0]; data[1] -= other[1]; data[2] -= other[2];
        data[3] -= other[3]; data[4] -= other[4]; data[5] -= other[5];
        return *this;
    }

    SymTensorO2& operator*=(CfemReal scalar) {
        data[0] *= scalar; data[1] *= scalar; data[2] *= scalar;
        data[3] *= scalar; data[4] *= scalar; data[5] *= scalar;
        return *this;
    }

    SymTensorO2& operator/=(CfemReal scalar) {
        CFEM_ASSERT(std::abs(scalar) > cfem::ZERO_TOLERANCE);
        data[0] /= scalar; data[1] /= scalar; data[2] /= scalar;
        data[3] /= scalar; data[4] /= scalar; data[5] /= scalar;
        return *this;
    }

    bool operator==(const SymTensorO2& other) const {
        for (int i = 0; i < 6; ++i)
            if (data[i] != other[i]) return false;
        return true;
    }

    bool operator!=(const SymTensorO2& other) const {
        return !(*this == other);
    }

    /**
     * @brief Calcule la forme quadratique M^T * (*this) * M.
     * @param M La matrice de transformation (ex: invJacobian).
     * @return Un nouveau SymTensorO2 (la symétrie est préservée).
     * * Utilisé pour transformer les Hessiennes : H_phys = J^-T * R_ref * J^-1
     */
    inline SymTensorO2 quadform(const TensorO2& M) const {
        SymTensorO2 res;
        
        // 1. On calcule d'abord le produit intermédiaire T = (*this) * M
        // Rappel stockage : [0]=xx, [1]=yy, [2]=zz, [3]=yz, [4]=xz, [5]=xy
        TensorO2 T;
        
        // Ligne 0 de (*this) : [xx, xy, xz]
        T(0,0) = data[0]*M(0,0) + data[5]*M(1,0) + data[4]*M(2,0);
        T(0,1) = data[0]*M(0,1) + data[5]*M(1,1) + data[4]*M(2,1);
        T(0,2) = data[0]*M(0,2) + data[5]*M(1,2) + data[4]*M(2,2);
        
        // Ligne 1 de (*this) : [xy, yy, yz]
        T(1,0) = data[5]*M(0,0) + data[1]*M(1,0) + data[3]*M(2,0);
        T(1,1) = data[5]*M(0,1) + data[1]*M(1,1) + data[3]*M(2,1);
        T(1,2) = data[5]*M(0,2) + data[1]*M(1,2) + data[3]*M(2,2);
        
        // Ligne 2 de (*this) : [xz, yz, zz]
        T(2,0) = data[4]*M(0,0) + data[3]*M(1,0) + data[2]*M(2,0);
        T(2,1) = data[4]*M(0,1) + data[3]*M(1,1) + data[2]*M(2,1);
        T(2,2) = data[4]*M(0,2) + data[3]*M(1,2) + data[2]*M(2,2);

        // 2. On calcule res = M^T * T
        // On ne calcule que les 6 composantes symétriques
        res[0] = M(0,0)*T(0,0) + M(1,0)*T(1,0) + M(2,0)*T(2,0); // xx
        res[1] = M(0,1)*T(0,1) + M(1,1)*T(1,1) + M(2,1)*T(2,1); // yy
        res[2] = M(0,2)*T(0,2) + M(1,2)*T(1,2) + M(2,2)*T(2,2); // zz
        
        res[3] = M(0,1)*T(0,2) + M(1,1)*T(1,2) + M(2,1)*T(2,2); // yz
        res[4] = M(0,0)*T(0,2) + M(1,0)*T(1,2) + M(2,0)*T(2,2); // xz
        res[5] = M(0,0)*T(0,1) + M(1,0)*T(1,1) + M(2,0)*T(2,1); // xy
        
        return res;
    }

    // --- Physical Invariants ---

    /** @brief Trace of the tensor (sum of diagonals). */
    inline CfemReal trace() const {
        return data[0] + data[1] + data[2];
    }

    /** @brief Mean (hydrostatic) component: p = tr(sigma)/3 */
    inline CfemReal mean() const {
        return trace() / 3.0;
    }

    /** @brief Returns the deviatoric part: s = sigma - p*I */
    SymTensorO2 deviator() const {
        CfemReal p = mean();
        return {data[0]-p, data[1]-p, data[2]-p, data[3], data[4], data[5]};
    }

    /** @brief Double contraction with another symmetric tensor (A : B). */
    inline CfemReal dotDot(const SymTensorO2& other) const {
        // xx*xx + yy*yy + zz*zz + 2*(yz*yz + xz*xz + xy*xy)
        return data[0]*other[0] + data[1]*other[1] + data[2]*other[2] +
               2.0 * (data[3]*other[3] + data[4]*other[4] + data[5]*other[5]);
    }

    /** @brief Second invariant of deviatoric tensor: J2 = 0.5 * (s : s) */
    inline CfemReal J2() const {
        SymTensorO2 s = deviator();
        return 0.5 * s.dotDot(s);
    }

    /** @brief Computes Von Mises equivalent stress (if this is a stress tensor). */
    CfemReal vonMises() const {
        CfemReal t = trace() / 3.0;
        CfemReal s_xx = data[0] - t;
        CfemReal s_yy = data[1] - t;
        CfemReal s_zz = data[2] - t;
        return std::sqrt(1.5 * (s_xx*s_xx + s_yy*s_yy + s_zz*s_zz + 
                         2.0*(data[3]*data[3] + data[4]*data[4] + data[5]*data[5])));
    }

    /** @brief Norm of the tensor (Frobenius norm) */
    inline CfemReal norm() const {
        return std::sqrt(this->dotDot(*this));
    }


};

// --- Global Ops ---

inline SymTensorO2 operator*(CfemReal scalar, const SymTensorO2& T) {
    return T * scalar;
}

/**
 * @brief Produit TenseurO2 * SymTensorO2.
 * Résultat : Un Tenseur d'ordre 2 général (non symétrique).
 * Algorithme : C_ij = A_ik * S_kj
 */
inline TensorO2 operator*(const TensorO2& A, const SymTensorO2& S) {
    TensorO2 res;

    res(0, 0) = A(0, 0) * S.xx() + A(0, 1) * S.xy() + A(0, 2) * S.xz();
    res(0, 1) = A(0, 0) * S.xy() + A(0, 1) * S.yy() + A(0, 2) * S.yz();
    res(0, 2) = A(0, 0) * S.xz() + A(0, 1) * S.yz() + A(0, 2) * S.zz();

    res(1, 0) = A(1, 0) * S.xx() + A(1, 1) * S.xy() + A(1, 2) * S.xz();
    res(1, 1) = A(1, 0) * S.xy() + A(1, 1) * S.yy() + A(1, 2) * S.yz();
    res(1, 2) = A(1, 0) * S.xz() + A(1, 1) * S.yz() + A(1, 2) * S.zz();

    res(2, 0) = A(2, 0) * S.xx() + A(2, 1) * S.xy() + A(2, 2) * S.xz();
    res(2, 1) = A(2, 0) * S.xy() + A(2, 1) * S.yy() + A(2, 2) * S.yz();
    res(2, 2) = A(2, 0) * S.xz() + A(2, 1) * S.yz() + A(2, 2) * S.zz();

    return res;
}

/**
 * @brief Produit SymTensorO2 * TensorO2.
 * Résultat : Un TensorO2 général.
 * Algorithme : C_ij = S_ik * A_kj
 */
inline TensorO2 operator*(const SymTensorO2& S, const TensorO2& A) {
    TensorO2 res;

    // Ligne 0 de S : [xx, xy, xz]
    res(0, 0) = S.xx() * A(0, 0) + S.xy() * A(1, 0) + S.xz() * A(2, 0);
    res(0, 1) = S.xx() * A(0, 1) + S.xy() * A(1, 1) + S.xz() * A(2, 1);
    res(0, 2) = S.xx() * A(0, 2) + S.xy() * A(1, 2) + S.xz() * A(2, 2);

    // Ligne 1 de S : [yx, yy, yz] -> [xy, yy, yz]
    res(1, 0) = S.xy() * A(0, 0) + S.yy() * A(1, 0) + S.yz() * A(2, 0);
    res(1, 1) = S.xy() * A(0, 1) + S.yy() * A(1, 1) + S.yz() * A(2, 1);
    res(1, 2) = S.xy() * A(0, 2) + S.yy() * A(1, 2) + S.yz() * A(2, 2);

    // Ligne 2 de S : [zx, zy, zz] -> [xz, yz, zz]
    res(2, 0) = S.xz() * A(0, 0) + S.yz() * A(1, 0) + S.zz() * A(2, 0);
    res(2, 1) = S.xz() * A(0, 1) + S.yz() * A(1, 1) + S.zz() * A(2, 1);
    res(2, 2) = S.xz() * A(0, 2) + S.yz() * A(1, 2) + S.zz() * A(2, 2);

    return res;
}

/**
 * @brief Produit SymTensorO2 * SymTensorO2.
 * Résultat : Un TensorO2 général.
 */
inline TensorO2 operator*(const SymTensorO2& A, const SymTensorO2& B) {
    TensorO2 res;

    // Ligne 0 de A : [xx, xy, xz] par colonnes de B
    res(0, 0) = A.xx() * B.xx() + A.xy() * B.xy() + A.xz() * B.xz();
    res(0, 1) = A.xx() * B.xy() + A.xy() * B.yy() + A.xz() * B.yz();
    res(0, 2) = A.xx() * B.xz() + A.xy() * B.yz() + A.xz() * B.zz();

    // Ligne 1 de A : [xy, yy, yz] par colonnes de B
    res(1, 0) = A.xy() * B.xx() + A.yy() * B.xy() + A.yz() * B.xz();
    res(1, 1) = A.xy() * B.xy() + A.yy() * B.yy() + A.yz() * B.yz();
    res(1, 2) = A.xy() * B.xz() + A.yy() * B.yz() + A.yz() * B.zz();

    // Ligne 2 de A : [xz, yz, zz] par colonnes de B
    res(2, 0) = A.xz() * B.xx() + A.yz() * B.xy() + A.zz() * B.xz();
    res(2, 1) = A.xz() * B.xy() + A.yz() * B.yy() + A.zz() * B.yz();
    res(2, 2) = A.xz() * B.xz() + A.yz() * B.yz() + A.zz() * B.zz();

    return res;
}

inline std::ostream& operator<<(std::ostream& os, const SymTensorO2& T) {
    os << "SymTensor2[xx:" << T[0] << ", yy:" << T[1] << ", zz:" << T[2] 
       << ", yz:" << T[3] << ", xz:" << T[4] << ", xy:" << T[5] << "]";
    return os;
}

}
#endif