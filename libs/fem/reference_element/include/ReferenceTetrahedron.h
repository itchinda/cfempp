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

#ifndef CFEM_REFERENCE_TETRAHEDRON
#define CFEM_REFERENCE_TETRAHEDRON

#include "ReferenceElement.h"
#include <cstring>

namespace cfem::reference_elem::utils{

    inline Vector3D projectFacetRefCoordToTetraRefCoord(CfemInt localFacetId, 
                                                        const Vector3D& xiFacet) 
    {
        const CfemReal u = xiFacet.x;
        const CfemReal v = xiFacet.y;

        CfemReal n0 = 1.0 - u - v;
        switch (localFacetId) {
            case 0: return Vector3D(u, v, 0.0);                // Face {0,1,2} : z=0
            case 1: return Vector3D(u, 0.0, v);                // Face {0,1,3} : y=0
            case 2: return Vector3D(n0, u, v);                 // Face {1,2,3} : x+y+z=1
            case 3: return Vector3D(0.0, n0, v);               // Face {2,0,3} : x=0
            default: CFEM_THROW("Invalid facet ID for Tetra");
        }
        return Vector3D{}; // Make the compiler happy
    }

    inline constexpr std::array<std::array<CfemInt,3>,4> tet4_face_nodes = {{
        {0, 2, 1}, 
        {0, 1, 3}, 
        {1, 2, 3}, 
        {2, 0, 3}
    }};

    inline constexpr std::array<std::array<CfemInt,6>,4> tet10_face_nodes = {{
        {0, 2, 1, 6, 5, 4}, 
        {0, 1, 3, 4, 8, 7}, 
        {1, 2, 3, 5, 9, 8}, 
        {2, 0, 3, 6, 7, 9}
    }};

    inline constexpr std::array<std::array<CfemInt,1>,4> tet1_face_nodes = {{
        {0}, {0}, {0}, {0}
    }};

    inline void computeTetraBarycentric(const Vector3D& xi,
                                        CfemReal& L0,
                                        CfemReal& L1,
                                        CfemReal& L2,
                                        CfemReal& L3) noexcept
    {
        L1 = xi.x;
        L2 = xi.y;
        L3 = xi.z;
        L0 = 1.0 - L1 - L2 - L3;
    }

}

namespace cfem::reference_elem::kernels {

    // Note : L'utilisation de __restrict (ou __restrict__ sous GCC/Clang) 
    // garantit au compilateur qu'il n'y a pas d'aliasing entre les pointeurs.
    
    struct P1TetraKernel {
        [[gnu::always_inline]] static inline void values(CfemReal L0, CfemReal L1, CfemReal L2, CfemReal L3, 
                                  CfemReal* __restrict N) noexcept
        {
            N[0] = L0; N[1] = L1; N[2] = L2; N[3] = L3;
        }

        [[gnu::always_inline]] static inline void grads(CfemReal* __restrict dx, 
                                 CfemReal* __restrict dy, 
                                 CfemReal* __restrict dz) noexcept
        {
            dx[0] = -1.0; dx[1] = 1.0; dx[2] = 0.0; dx[3] = 0.0;
            dy[0] = -1.0; dy[1] = 0.0; dy[2] = 1.0; dy[3] = 0.0;
            dz[0] = -1.0; dz[1] = 0.0; dz[2] = 0.0; dz[3] = 1.0;
        }

        [[gnu::always_inline]] static inline void hess(CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept
        {
            for(CfemInt i = 0; i < 4; ++i) {
                d2x[i] = d2y[i] = d2z[i] = dxy[i] = dxz[i] = dyz[i] = 0.0;
            }
        }
    };

    struct P2TetraKernel {
        [[gnu::always_inline]] static inline void values(CfemReal L0, CfemReal L1, CfemReal L2, CfemReal L3, 
                                  CfemReal* __restrict N) noexcept
        {
            N[0]=L0*(2*L0-1); N[1]=L1*(2*L1-1); N[2]=L2*(2*L2-1); N[3]=L3*(2*L3-1);
            N[4]=4*L0*L1;     N[5]=4*L1*L2;     N[6]=4*L2*L0;
            N[7]=4*L0*L3;     N[8]=4*L1*L3;     N[9]=4*L2*L3;
        }

        [[gnu::always_inline]] static inline void grads(CfemReal L0, CfemReal L1, CfemReal L2, CfemReal L3, 
                                 CfemReal* __restrict dx, 
                                 CfemReal* __restrict dy, 
                                 CfemReal* __restrict dz) noexcept
        {
            dx[0]=1-4*L0; dx[1]=4*L1-1; dx[2]=0;      dx[3]=0;      dx[4]=4*(L0-L1); dx[5]=4*L2; dx[6]=-4*L2;     dx[7]=-4*L3;     dx[8]=4*L3; dx[9]=0;
            dy[0]=1-4*L0; dy[1]=0;      dy[2]=4*L2-1; dy[3]=0;      dy[4]=-4*L1;     dy[5]=4*L1; dy[6]=4*(L0-L2); dy[7]=-4*L3;     dy[8]=0;    dy[9]=4*L3;
            dz[0]=1-4*L0; dz[1]=0;      dz[2]=0;      dz[3]=4*L3-1; dz[4]=-4*L1;     dz[5]=0;    dz[6]=-4*L2;     dz[7]=4*(L0-L3); dz[8]=4*L1; dz[9]=4*L2;
        }

        [[gnu::always_inline]] static inline void hess(CfemReal* __restrict d2x, CfemReal* __restrict d2y, CfemReal* __restrict d2z, 
                                CfemReal* __restrict dxy, CfemReal* __restrict dxz, CfemReal* __restrict dyz) noexcept
        {
            // Lignes = composantes (xx, yy, zz, xy, xz, yz)
            // Colonnes = nœuds (0 à 9)
            static constexpr CfemReal H[6][10] = {
                {4.0, 4.0, 0.0, 0.0,-8.0, 0.0, 0.0, 0.0, 0.0, 0.0}, // d2x2 (CORRIGÉ: Noeuds 1, 6, 7)
                {4.0, 0.0, 4.0, 0.0, 0.0, 0.0,-8.0, 0.0, 0.0, 0.0}, // d2y2
                {4.0, 0.0, 0.0, 4.0, 0.0, 0.0, 0.0,-8.0, 0.0, 0.0}, // d2z2
                {4.0, 0.0, 0.0, 0.0,-4.0, 4.0,-4.0, 0.0, 0.0, 0.0}, // dxy
                {4.0, 0.0, 0.0, 0.0,-4.0, 0.0, 0.0,-4.0, 4.0, 0.0}, // dxz
                {4.0, 0.0, 0.0, 0.0, 0.0, 0.0,-4.0,-4.0, 0.0, 4.0}  // dyz
            };
            
            for(CfemInt i = 0; i < 10; ++i) {
                d2x[i] = H[0][i]; 
                d2y[i] = H[1][i]; 
                d2z[i] = H[2][i]; 
                dxy[i] = H[3][i]; 
                dxz[i] = H[4][i]; 
                dyz[i] = H[5][i];
            }
        }
    };
}

namespace cfem {

    /**
     * @brief Tétraèdre Linéaire (4 nœuds)
     */
    class P1ReferenceTetrahedron : public ReferenceElement {
    public:
        P1ReferenceTetrahedron() {
            m_numNodes = 4;
            m_dimension = 3;
            
            m_nodes = {
                Vector3D{0.0, 0.0, 0.0}, 
                Vector3D{1.0, 0.0, 0.0}, 
                Vector3D{0.0, 1.0, 0.0}, 
                Vector3D{0.0, 0.0, 1.0} 
            };
        }

        void evaluateShapesValues(const Vector3D& xi, CfemReal* values) const override {
            CfemReal L0, L1, L2, L3;
            reference_elem::utils::computeTetraBarycentric(xi, L0, L1, L2, L3);
            reference_elem::kernels::P1TetraKernel::values(L0, L1, L2, L3, values);
        }

        void evaluateShapesDerivatives(const Vector3D&, CfemReal* dN_dxi, CfemReal* dN_deta, CfemReal* dN_dzeta) const override {
            reference_elem::kernels::P1TetraKernel::grads(dN_dxi, dN_deta, dN_dzeta);
        }

        void evaluateShapesSecondDerivatives(const Vector3D&,
                                            CfemReal* d2N_dxi2, 
                                            CfemReal* d2N_deta2, 
                                            CfemReal* d2N_dzeta2,
                                            CfemReal* d2N_dxideta, 
                                            CfemReal* d2N_dxidzeta, 
                                            CfemReal* d2N_detadzeta) const override
        {
            reference_elem::kernels::P1TetraKernel::hess(d2N_dxi2, d2N_deta2, d2N_dzeta2, 
                                            d2N_dxideta, d2N_dxidzeta, d2N_detadzeta);
        }

        // =====================================================================
        // INTERFACES BATCH
        // =====================================================================
        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, DynamicMatrix<CfemReal>& values) const override {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                const Vector3D& xi = xiPoints[q];
                CfemReal L0, L1, L2, L3;
                reference_elem::utils::computeTetraBarycentric(xi, L0, L1, L2, L3);
                reference_elem::kernels::P1TetraKernel::values(L0, L1, L2, L3, values[q].data());
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, 
                                            DynamicMatrix<CfemReal>& dN_dxi, 
                                            DynamicMatrix<CfemReal>& dN_deta, 
                                            DynamicMatrix<CfemReal>& dN_dzeta) const override {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());

            for (CfemInt q = 0; q < nQp; ++q) {
                reference_elem::kernels::P1TetraKernel::grads(dN_dxi[q].data(), dN_deta[q].data(), dN_dzeta[q].data());
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D>& xiPoints,
                                                DynamicMatrix<CfemReal>& d2N_dxi2, 
                                                DynamicMatrix<CfemReal>& d2N_deta2, 
                                                DynamicMatrix<CfemReal>& d2N_dzeta2,
                                                DynamicMatrix<CfemReal>& d2N_dxideta, 
                                                DynamicMatrix<CfemReal>& d2N_dxidzeta, 
                                                DynamicMatrix<CfemReal>& d2N_detadzeta) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            if (nQp == 0) return;

            reference_elem::kernels::P1TetraKernel::hess(
                d2N_dxi2[0].data(), d2N_deta2[0].data(), d2N_dzeta2[0].data(),
                d2N_dxideta[0].data(), d2N_dxidzeta[0].data(), d2N_detadzeta[0].data()
            );

            const size_t bytesToCopy = m_numNodes * sizeof(CfemReal);
            for (CfemInt q = 1; q < nQp; ++q) {
                std::memcpy(d2N_dxi2[q].data(),      d2N_dxi2[0].data(),      bytesToCopy);
                std::memcpy(d2N_deta2[q].data(),     d2N_deta2[0].data(),     bytesToCopy);
                std::memcpy(d2N_dzeta2[q].data(),    d2N_dzeta2[0].data(),    bytesToCopy);
                std::memcpy(d2N_dxideta[q].data(),   d2N_dxideta[0].data(),   bytesToCopy);
                std::memcpy(d2N_dxidzeta[q].data(),  d2N_dxidzeta[0].data(),  bytesToCopy);
                std::memcpy(d2N_detadzeta[q].data(), d2N_detadzeta[0].data(), bytesToCopy);
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToTetraRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tet4_face_nodes[iFace];
        }
    };

    /**
     * @brief Tétraèdre Quadratique (10 nœuds)
     */
    class P2ReferenceTetrahedron : public ReferenceElement {
    public:
        P2ReferenceTetrahedron() {
            m_numNodes = 10;
            m_dimension = 3;
            
            m_nodes = {
                Vector3D{0.0, 0.0 ,0.0},
                Vector3D{1.0, 0.0 ,0.0}, 
                Vector3D{0.0, 1.0 ,0.0}, 
                Vector3D{0.0, 0.0 ,1.0},
                // Milieux d'arêtes (Standard convention: 01, 12, 20, 03, 13, 23)
                Vector3D{0.5, 0.0, 0.0}, 
                Vector3D{0.5, 0.5, 0.0}, 
                Vector3D{0.0, 0.5, 0.0},
                Vector3D{0.0, 0.0, 0.5}, 
                Vector3D{0.5, 0.0, 0.5}, 
                Vector3D{0.0, 0.5, 0.5}
            };
        }

        // =========================================================================
        // FIRST DERIVATIVES
        // =========================================================================
        void evaluateShapesValues(const Vector3D& xi, CfemReal* values) const override
        {
            CfemReal L0, L1, L2, L3;
            reference_elem::utils::computeTetraBarycentric(xi, L0, L1, L2, L3);
            reference_elem::kernels::P2TetraKernel::values(L0, L1, L2, L3, values);
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D>& xiPoints, 
                                                            DynamicMatrix<CfemReal>& values) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                const Vector3D& xi = xiPoints[q];
                CfemReal L0, L1, L2, L3;
                reference_elem::utils::computeTetraBarycentric(xi, L0, L1, L2, L3);
                
                reference_elem::kernels::P2TetraKernel::values(L0, L1, L2, L3, values[q].data());
            }
        }
    
        // =========================================================================
        // FIRST DERIVATIVES
        // =========================================================================
        void evaluateShapesDerivatives(const Vector3D& xi, 
                                    CfemReal* dN_dxi, 
                                    CfemReal* dN_deta, 
                                    CfemReal* dN_dzeta) const override
        {
            CfemReal L0, L1, L2, L3;
            reference_elem::utils::computeTetraBarycentric(xi, L0, L1, L2, L3);
            reference_elem::kernels::P2TetraKernel::grads(L0, L1, L2, L3, dN_dxi, dN_deta, dN_dzeta);
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D>& xiPoints, 
                                            DynamicMatrix<CfemReal>& dN_dxi, 
                                            DynamicMatrix<CfemReal>& dN_deta, 
                                            DynamicMatrix<CfemReal>& dN_dzeta) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                const Vector3D& xi = xiPoints[q];
                CfemReal L0, L1, L2, L3;
                reference_elem::utils::computeTetraBarycentric(xi, L0, L1, L2, L3);
                
                reference_elem::kernels::P2TetraKernel::grads(L0, L1, L2, L3, 
                                                    dN_dxi[q].data(), 
                                                    dN_deta[q].data(), 
                                                    dN_dzeta[q].data());
            }
        }

        // =========================================================================
        // SECOND DERIVATIVES
        // =========================================================================
        void evaluateShapesSecondDerivatives(const Vector3D& /*xi*/,
                                            CfemReal* d2N_dxi2, 
                                            CfemReal* d2N_deta2, 
                                            CfemReal* d2N_dzeta2,
                                            CfemReal* d2N_dxideta, 
                                            CfemReal* d2N_dxidzeta, 
                                            CfemReal* d2N_detadzeta) const override
        {
            reference_elem::kernels::P2TetraKernel::hess(d2N_dxi2, d2N_deta2, d2N_dzeta2, 
                                            d2N_dxideta, d2N_dxidzeta, d2N_detadzeta);
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D>& xiPoints,
                                                DynamicMatrix<CfemReal>& d2N_dxi2, 
                                                DynamicMatrix<CfemReal>& d2N_deta2, 
                                                DynamicMatrix<CfemReal>& d2N_dzeta2,
                                                DynamicMatrix<CfemReal>& d2N_dxideta, 
                                                DynamicMatrix<CfemReal>& d2N_dxidzeta, 
                                                DynamicMatrix<CfemReal>& d2N_detadzeta) const override
        {
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            if (nQp == 0) return;

            reference_elem::kernels::P2TetraKernel::hess(
                d2N_dxi2[0].data(), d2N_deta2[0].data(), d2N_dzeta2[0].data(),
                d2N_dxideta[0].data(), d2N_dxidzeta[0].data(), d2N_detadzeta[0].data()
            );

            const size_t bytesToCopy = m_numNodes * sizeof(CfemReal);
            for (CfemInt q = 1; q < nQp; ++q) {
                std::memcpy(d2N_dxi2[q].data(),      d2N_dxi2[0].data(),      bytesToCopy);
                std::memcpy(d2N_deta2[q].data(),     d2N_deta2[0].data(),     bytesToCopy);
                std::memcpy(d2N_dzeta2[q].data(),    d2N_dzeta2[0].data(),    bytesToCopy);
                std::memcpy(d2N_dxideta[q].data(),   d2N_dxideta[0].data(),   bytesToCopy);
                std::memcpy(d2N_dxidzeta[q].data(),  d2N_dxidzeta[0].data(),  bytesToCopy);
                std::memcpy(d2N_detadzeta[q].data(), d2N_detadzeta[0].data(), bytesToCopy);
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToTetraRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tet10_face_nodes[iFace];
        }
    };

    class MiniReferenceTetrahedron : public ReferenceElement {
    public:
        MiniReferenceTetrahedron() {
            m_numNodes = 5; // 4 sommets + 1 bulle
            m_dimension = 3;
            
            m_nodes = { 
                Vector3D{0.0, 0.0 ,0.0},
                Vector3D{1.0, 0.0 ,0.0},
                Vector3D{0.0, 1.0 ,0.0},
                Vector3D{0.0, 0.0 ,1.0},
                Vector3D{0.25, 0.25, 0.25}
            };
        }

        void evaluateShapesValues(const Vector3D& xi, CfemReal* values) const override 
        {
            CfemReal L0, L1, L2, L3;
            reference_elem::utils::computeTetraBarycentric(xi, L0, L1, L2, L3);

            // Bubble function: 256 * L0 * L1 * L2 * L3
            const CfemReal b = 256.0 * L0 * L1 * L2 * L3;

            values[4] = b;
            values[0] = L0 - b * 0.25;
            values[1] = L1 - b * 0.25;
            values[2] = L2 - b * 0.25;
            values[3] = L3 - b * 0.25;
        }

        void evaluateShapesDerivatives(const Vector3D& xi, 
                                 CfemReal* dN_dxi, 
                                 CfemReal* dN_deta, 
                                 CfemReal* dN_dzeta) const override 
        {
            CfemReal L0, L1, L2, L3;
            reference_elem::utils::computeTetraBarycentric(xi, L0, L1, L2, L3);

            // Bubble gradients
            const CfemReal db_dxi  = 256.0 * L2 * L3 * (L0 - L1);
            const CfemReal db_deta = 256.0 * L1 * L3 * (L0 - L2);
            const CfemReal db_dzeta = 256.0 * L1 * L2 * (L0 - L3);

            // Node 4 (Bubble)
            dN_dxi[4]   = db_dxi;
            dN_deta[4]  = db_deta;
            dN_dzeta[4] = db_dzeta;

            // Node 0
            dN_dxi[0]   = -1.0 - db_dxi * 0.25;
            dN_deta[0]  = -1.0 - db_deta * 0.25;
            dN_dzeta[0] = -1.0 - db_dzeta * 0.25;

            // Node 1
            dN_dxi[1]   =  1.0 - db_dxi * 0.25;
            dN_deta[1]  =  0.0 - db_deta * 0.25;
            dN_dzeta[1] =  0.0 - db_dzeta * 0.25;

            // Node 2
            dN_dxi[2]   =  0.0 - db_dxi * 0.25;
            dN_deta[2]  =  1.0 - db_deta * 0.25;
            dN_dzeta[2] =  0.0 - db_dzeta * 0.25;

            // Node 3
            dN_dxi[3]   =  0.0 - db_dxi * 0.25;
            dN_deta[3]  =  0.0 - db_deta * 0.25;
            dN_dzeta[3] =  1.0 - db_dzeta * 0.25;
        }

        void evaluateShapesSecondDerivatives(const Vector3D& xi,
                                       CfemReal* d2N_dxi2, 
                                       CfemReal* d2N_deta2, 
                                       CfemReal* d2N_dzeta2,
                                       CfemReal* d2N_dxideta, 
                                       CfemReal* d2N_dxidzeta, 
                                       CfemReal* d2N_detadzeta) const override 
        {
            CfemReal L0, L1, L2, L3;
            reference_elem::utils::computeTetraBarycentric(xi, L0, L1, L2, L3);

            // Second derivatives of the bubble function
            const CfemReal d2b_dxi2     = -512.0 * L2 * L3;
            const CfemReal d2b_deta2    = -512.0 * L1 * L3;
            const CfemReal d2b_dzeta2   = -512.0 * L1 * L2;
            const CfemReal d2b_dxideta  =  256.0 * L3 * (L0 - L1 - L2);
            const CfemReal d2b_dxidzeta =  256.0 * L2 * (L0 - L1 - L3);
            const CfemReal d2b_detadzeta=  256.0 * L1 * (L0 - L2 - L3);

            // Node 4 (Bubble)
            d2N_dxi2[4]      = d2b_dxi2; 
            d2N_deta2[4]     = d2b_deta2; 
            d2N_dzeta2[4]    = d2b_dzeta2;
            d2N_dxideta[4]   = d2b_dxideta; 
            d2N_dxidzeta[4]  = d2b_dxidzeta; 
            d2N_detadzeta[4] = d2b_detadzeta;

            // Nodes 0 to 3 inherit the bubble's correction divided by -4.0
            for (CfemInt i = 0; i < 4; ++i) {
                d2N_dxi2[i]      = -d2b_dxi2 * 0.25;
                d2N_deta2[i]     = -d2b_deta2 * 0.25;
                d2N_dzeta2[i]    = -d2b_dzeta2 * 0.25;
                d2N_dxideta[i]   = -d2b_dxideta * 0.25;
                d2N_dxidzeta[i]  = -d2b_dxidzeta * 0.25;
                d2N_detadzeta[i] = -d2b_detadzeta * 0.25;
            }
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, DynamicMatrix<CfemReal>& values) const override
        {
            // Constant piecewise field value
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                CfemReal L0, L1, L2, L3;
                reference_elem::utils::computeTetraBarycentric(xiPoints[q], L0, L1, L2, L3);

                // Bubble function: 256 * L0 * L1 * L2 * L3
                const CfemReal b = 256.0 * L0 * L1 * L2 * L3;

                values[q][4] = b;
                values[q][0] = L0 - b * 0.25;
                values[q][1] = L1 - b * 0.25;
                values[q][2] = L2 - b * 0.25;
                values[q][3] = L3 - b * 0.25;
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> & xiPoints, 
                                      DynamicMatrix<CfemReal> &dN_dxi, 
                                      DynamicMatrix<CfemReal> &dN_deta, 
                                      DynamicMatrix<CfemReal> &dN_dzeta) const override
        {
            // Constant function results in flat zero gradients everywhere
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                CfemReal L0, L1, L2, L3;
                reference_elem::utils::computeTetraBarycentric(xiPoints[q], L0, L1, L2, L3);

                // Bubble gradients
                const CfemReal db_dxi  = 256.0 * L2 * L3 * (L0 - L1);
                const CfemReal db_deta = 256.0 * L1 * L3 * (L0 - L2);
                const CfemReal db_dzeta = 256.0 * L1 * L2 * (L0 - L3);

                // Node 4 (Bubble)
                dN_dxi[q][4]   = db_dxi;
                dN_deta[q][4]  = db_deta;
                dN_dzeta[q][4] = db_dzeta;

                // Node 0
                dN_dxi[q][0]   = -1.0 - db_dxi * 0.25;
                dN_deta[q][0]  = -1.0 - db_deta * 0.25;
                dN_dzeta[q][0] = -1.0 - db_dzeta * 0.25;

                // Node 1
                dN_dxi[q][1]   =  1.0 - db_dxi * 0.25;
                dN_deta[q][1]  =  0.0 - db_deta * 0.25;
                dN_dzeta[q][1] =  0.0 - db_dzeta * 0.25;

                // Node 2
                dN_dxi[q][2]   =  0.0 - db_dxi * 0.25;
                dN_deta[q][2]  =  1.0 - db_deta * 0.25;
                dN_dzeta[q][2] =  0.0 - db_dzeta * 0.25;

                // Node 3
                dN_dxi[q][3]   =  0.0 - db_dxi * 0.25;
                dN_deta[q][3]  =  0.0 - db_deta * 0.25;
                dN_dzeta[q][3] =  1.0 - db_dzeta * 0.25;
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints,
                                            DynamicMatrix<CfemReal> &d2N_dxi2, 
                                            DynamicMatrix<CfemReal> &d2N_deta2, 
                                            DynamicMatrix<CfemReal> &d2N_dzeta2,
                                            DynamicMatrix<CfemReal> &d2N_dxideta, 
                                            DynamicMatrix<CfemReal> &d2N_dxidzeta, 
                                            DynamicMatrix<CfemReal> &d2N_detadzeta) const override
        {
            // Flat null second derivatives
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                d2N_dxi2[q][0]      = 0.0;
                d2N_deta2[q][0]     = 0.0;
                d2N_dzeta2[q][0]    = 0.0;
                d2N_dxideta[q][0]   = 0.0;
                d2N_dxidzeta[q][0]  = 0.0;
                d2N_detadzeta[q][0] = 0.0;
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToTetraRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tet4_face_nodes[iFace];
        }
    };

class P0ReferenceTetrahedron : public ReferenceElement {
    public:
        P0ReferenceTetrahedron() {
            m_numNodes = 1;
            m_dimension = 3;
        
            m_nodes = { Vector3D{0.25, 0.25, 0.25} };
        }

        void evaluateShapesValues(const Vector3D&, CfemReal* values) const override 
        {
            values[0] = 1.0;
        }

        void evaluateShapesDerivatives(const Vector3D&, 
                                 CfemReal* dN_dxi, 
                                 CfemReal* dN_deta, 
                                 CfemReal* dN_dzeta) const override 
        {
            dN_dxi[0]   = 0.0;
            dN_deta[0]  = 0.0;
            dN_dzeta[0] = 0.0;
        }

        void evaluateShapesSecondDerivatives(const Vector3D&,
                                       CfemReal* d2N_dxi2, 
                                       CfemReal* d2N_deta2, 
                                       CfemReal* d2N_dzeta2,
                                       CfemReal* d2N_dxideta, 
                                       CfemReal* d2N_dxidzeta, 
                                       CfemReal* d2N_detadzeta) const override 
        {
            d2N_dxi2[0]      = 0.0;
            d2N_deta2[0]     = 0.0;
            d2N_dzeta2[0]    = 0.0;
            d2N_dxideta[0]   = 0.0;
            d2N_dxidzeta[0]  = 0.0;
            d2N_detadzeta[0] = 0.0;
        }

        void evaluateShapesValuesBatch(const DynamicVector<Vector3D> &xiPoints, DynamicMatrix<CfemReal>& values) const override
        {
            // Constant piecewise field value
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                values[q][0] = 1.0;
            }
        }

        void evaluateShapesDerivativesBatch(const DynamicVector<Vector3D> & xiPoints, 
                                      DynamicMatrix<CfemReal> &dN_dxi, 
                                      DynamicMatrix<CfemReal> &dN_deta, 
                                      DynamicMatrix<CfemReal> &dN_dzeta) const override
        {
            // Constant function results in flat zero gradients everywhere
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                dN_dxi[q][0]   = 0.0;
                dN_deta[q][0]  = 0.0;
                dN_dzeta[q][0] = 0.0;
            }
        }

        void evaluateShapesSecondDerivativesBatch(const DynamicVector<Vector3D> &xiPoints,
                                            DynamicMatrix<CfemReal> &d2N_dxi2, 
                                            DynamicMatrix<CfemReal> &d2N_deta2, 
                                            DynamicMatrix<CfemReal> &d2N_dzeta2,
                                            DynamicMatrix<CfemReal> &d2N_dxideta, 
                                            DynamicMatrix<CfemReal> &d2N_dxidzeta, 
                                            DynamicMatrix<CfemReal> &d2N_detadzeta) const override
        {
            // Flat null second derivatives
            const CfemInt nQp = static_cast<CfemInt>(xiPoints.size());
            for (CfemInt q = 0; q < nQp; ++q) {
                d2N_dxi2[q][0]      = 0.0;
                d2N_deta2[q][0]     = 0.0;
                d2N_dzeta2[q][0]    = 0.0;
                d2N_dxideta[q][0]   = 0.0;
                d2N_dxidzeta[q][0]  = 0.0;
                d2N_detadzeta[q][0] = 0.0;
            }
        }

        Vector3D projectFacetRefCoordToCellRefCoord(CfemInt localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToTetraRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(CfemInt iFace) const override {
            return cfem::reference_elem::utils::tet1_face_nodes[iFace];
        }
    };
}
#endif