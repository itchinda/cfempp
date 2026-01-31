/**
 * @file SymTensorO2.h
 * @brief Specialized storage for symmetric 3x3 tensors (6 components).
 */

#ifndef CFEM_SYM_TENSOR2
#define CFEM_SYM_TENSOR2

#include <cmath>
#include <ostream>

#include "cfemutils.h"
#include "Vector3D.h"

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

    // --- Accessors ---
    
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

    // --- Basic Math Ops ---

    static SymTensorO2 Identity() {
        return {1,1,1,0,0,0};
    }

    static SymTensorO2 Zero() {
        return {0,0,0,0,0,0};
    }

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

inline std::ostream& operator<<(std::ostream& os, const SymTensorO2& T) {
    os << "SymTensor2[xx:" << T[0] << ", yy:" << T[1] << ", zz:" << T[2] 
       << ", yz:" << T[3] << ", xz:" << T[4] << ", xy:" << T[5] << "]";
    return os;
}

#endif