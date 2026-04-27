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

#ifndef CFEM_REFERENCE_TRIANGLE
#define CFEM_REFERENCE_TRIANGLE

#include "ReferenceElement.h"
#include "Vector3D.h"

namespace cfem::reference_elem::utils{

    inline Vector3D projectFacetRefCoordToTriRefCoord(CfemInt localFacetId, 
                                                      const Vector3D& xiFacet) 
    {
        const CfemReal u = xiFacet.x; // [-1, 1] pour une ligne
        CfemReal n1 = 0.5 * (1.0 - u);
        CfemReal n2 = 0.5 * (1.0 + u);

            switch (localFacetId) {
                case 0: return Vector3D(n2, 0.0, 0.0); // Face {0,1}
                case 1: return Vector3D(n1, n2, 0.0); // Face {1,2}
                case 2: return Vector3D(0.0, n1, 0.0); // Face {2,0}
                default: CFEM_ERROR << "Invalid facet ID for Triangle";
            }
        return Vector3D{}; // Make the compiler happy
    }
   
    static const std::vector<std::vector<int>> tri1_face_nodes = {
        {0}, {0}, {0}
    };

    static const std::vector<std::vector<int>> tri3_face_nodes = {
        {0, 1}, {1, 2}, {2, 0}
    };
    
    static const std::vector<std::vector<int>> tri6_face_nodes = {
        {0, 1, 3}, {1, 2, 4}, {2, 0, 5}
    };

}

namespace cfem {
    class P1ReferenceTriangle : public ReferenceElement {

    public:
        
        P1ReferenceTriangle() {
            m_numNodes = 3;
            m_dimension = 2;
            m_nodes.resize(m_numNodes);

            m_nodes = { {0.0, 0.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 1.0, 0.0} };
        }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = getNumNodes();
    
            if (out.values.size() != n) {
                out.setup(n); 
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            // N0 = 1 - xi - eta, N1 = xi, N2 = eta
            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = 1.0 - xi.x - xi.y;
                out.values[1] = xi.x;
                out.values[2] = xi.y;
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Gradient) || hasFlag(requested, ShapeEvaluationFlags::Jacobian)  || hasFlag(requested, ShapeEvaluationFlags::Normal)) {
                out.gradients[0] = Vector3D{-1.0, -1.0, 0.0};
                out.gradients[1] = Vector3D{ 1.0,  0.0, 0.0};
                out.gradients[2] = Vector3D{ 0.0,  1.0, 0.0};
            }

            // Hessiens par rapport à (xi, eta)
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                for (auto& H : out.hessians) H.setZero();
            }

            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::tri3_face_nodes[iFace];
        }

    };


    class P2ReferenceTriangle : public ReferenceElement {
    public:
        P2ReferenceTriangle() {
            m_numNodes = 6;
            m_dimension = 2;
            m_nodes.resize(m_numNodes);
            // Main nodes
            m_nodes[0] = Vector3D{0.0, 0.0, 0.0};
            m_nodes[1] = Vector3D{1.0, 0.0, 0.0};
            m_nodes[2] = Vector3D{0.0, 1.0, 0.0};
            // Edges Middles
            m_nodes[3] = Vector3D{0.5, 0.0, 0.0};
            m_nodes[4] = Vector3D{0.5, 0.5, 0.0};
            m_nodes[5] = Vector3D{0.0, 0.5, 0.0};
        }
        
        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 6;
            if (out.values.size() != n) out.setup(n);
            
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            const CfemReal L1 = xi.x;
            const CfemReal L2 = xi.y;
            const CfemReal L0 = 1.0 - L1 - L2;

            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = L0 * (2.0 * L0 - 1.0);
                out.values[1] = L1 * (2.0 * L1 - 1.0);
                out.values[2] = L2 * (2.0 * L2 - 1.0);
                out.values[3] = 4.0 * L0 * L1;
                out.values[4] = 4.0 * L1 * L2;
                out.values[5] = 4.0 * L2 * L0;
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                out.gradients[0] = Vector3D{-3.0 + 4.0*L1 + 4.0*L2, -3.0 + 4.0*L1 + 4.0*L2, 0.0};
                out.gradients[1] = Vector3D{4.0*L1 - 1.0, 0.0, 0.0};
                out.gradients[2] = Vector3D{0.0, 4.0*L2 - 1.0, 0.0};
                out.gradients[3] = Vector3D{4.0 - 8.0*L1 - 4.0*L2, -4.0*L1, 0.0};
                out.gradients[4] = Vector3D{4.0*L2, 4.0*L1, 0.0};
                out.gradients[5] = Vector3D{-4.0*L2, 4.0 - 4.0*L1 - 8.0*L2, 0.0};
            }

            // --- Hessiennes (Dérivées secondes constantes) ---
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                for (int i = 0; i < n; ++i) {
                    out.hessians[i].setZero(); 
                }

                out.hessians[0].xx() = 4.0; out.hessians[0].yy() = 4.0; out.hessians[0].xy()=4.0;
                out.hessians[1].xx() = 4.0;
                out.hessians[2].yy() = 4.0;
                out.hessians[3].xx() =-8.0; out.hessians[3].xy() = -4.0;
                out.hessians[4].xy() = 4.0; out.hessians[4].xy() =  4.0;
                out.hessians[5].yy() =-8.0; out.hessians[5].xy() = -4.0;
            }

            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override{
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet); // Make the compiler happy
        }

        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::tri6_face_nodes[iFace];
        }
    };

    class MiniReferenceTriangle : public ReferenceElement {
    public:
        MiniReferenceTriangle() {
            m_numNodes = 4; // 3 sommets + 1 bulle
            m_dimension = 2;
            m_nodes.resize(m_numNodes);

            m_nodes[0] = {0.0, 0.0, 0.0};
            m_nodes[1] = {1.0, 0.0, 0.0};
            m_nodes[2] = {0.0, 1.0, 0.0};
            m_nodes[3] = {1.0/3.0, 1.0/3.0, 0.0};
        }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 4;
            if (out.values.size() != n) out.setup(n);
            if (hasFlag(requested, ShapeEvaluationFlags::Hessian) && out.hessians.size() != n) {
                out.hessians.assign(n, SymTensorO2());
            }

            const CfemReal x = xi.x;
            const CfemReal y = xi.y;
            const CfemReal L0 = 1.0 - x - y;
            const CfemReal L1 = x;
            const CfemReal L2 = y;

            // Fonction Bulle (Normalisée pour valoir 1.0 au centre)
            const CfemReal b = 27.0 * L0 * L1 * L2;

            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[3] = b;
                out.values[0] = L0 - b / 3.0;
                out.values[1] = L1 - b / 3.0;
                out.values[2] = L2 - b / 3.0;
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                // Gradient de la bulle
                CfemReal db_dx = 27.0 * L2 * (1.0 - 2.0*L1 - L2);
                CfemReal db_dy = 27.0 * L1 * (1.0 - L1 - 2.0*L2);

                out.gradients[3] = Vector3D{db_dx, db_dy, 0.0};
                out.gradients[0] = Vector3D{-1.0 - db_dx / 3.0, -1.0 - db_dy / 3.0, 0.0};
                out.gradients[1] = Vector3D{ 1.0 - db_dx / 3.0,  0.0 - db_dy / 3.0, 0.0};
                out.gradients[2] = Vector3D{ 0.0 - db_dx / 3.0,  1.0 - db_dy / 3.0, 0.0};
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                for (auto& H : out.hessians) H.setZero();

                // Hessienne de la bulle
                CfemReal d2b_dxx = -54.0 * L2;
                CfemReal d2b_dyy = -54.0 * L1;
                CfemReal d2b_dxy =  27.0 * (1.0 - 2.0*L1 - 2.0*L2);

                out.hessians[3].xx() = d2b_dxx;
                out.hessians[3].yy() = d2b_dyy;
                out.hessians[3].xy() = d2b_dxy;

                // Les fonctions P1 n'ont pas de Hessienne propre, on hérite juste de la correction bulle
                for (int i = 0; i < 3; ++i) {
                    out.hessians[i].xx() = -d2b_dxx / 3.0;
                    out.hessians[i].yy() = -d2b_dyy / 3.0;
                    out.hessians[i].xy() = -d2b_dxy / 3.0;
                }
            }
            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet);
        }

        // On réutilise tri3_face_nodes car la bulle n'apparaît pas sur les faces !
        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::tri3_face_nodes[iFace];
        }
    };

    class P0ReferenceTriangle : public ReferenceElement {
    public:
        P0ReferenceTriangle() {
            m_numNodes = 1;
            m_dimension = 2;
            m_nodes.resize(m_numNodes);

            m_nodes = { {1.0/3.0, 1.0/3.0, 0.0} };
        }

        void evaluate(const Vector3D& xi, ShapeEvaluationFlags requested, ShapeInfo& out) const override {
            const int n = 1;
            if (out.values.size() != n) {
                out.setup(n); 
            }

            // Valeur constante : N0 = 1
            if (hasFlag(requested, ShapeEvaluationFlags::Value)) {
                out.values[0] = 1.0;
            }

            // Gradient nul : dN0/dxi = 0, dN0/deta = 0
            if (hasFlag(requested, ShapeEvaluationFlags::Gradient)) {
                out.gradients[0] = Vector3D{0.0, 0.0, 0.0};
            }

            if (hasFlag(requested, ShapeEvaluationFlags::Hessian)) {
                if (out.hessians.size() != n) out.hessians.assign(n, SymTensorO2());
                out.hessians[0].setZero();
            }

            out.validFlags = requested;
        }

        Vector3D projectFacetRefCoordToCellRefCoord(int localFacetId, const Vector3D& xiFacet) const override {
            return cfem::reference_elem::utils::projectFacetRefCoordToTriRefCoord(localFacetId, xiFacet);
        }

        std::span<const CfemInt> getFaceNodes(int iFace) const override {
            return cfem::reference_elem::utils::tri1_face_nodes[iFace];
        }
    };


}  // namespace


#endif