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

#include <cmath>
#include <gtest/gtest.h>
#include <memory>
#include <vector>

#include "Quadrature.h"
#include "ReferenceElement.h"
#include "ReferenceQuad.h"
#include "ReferenceSegment.h"
#include "ReferenceTetrahedron.h"
#include "ReferenceTriangle.h"

using namespace cfem;

struct QuadratureTestParams {
    CellType cellType;
    CfemInt numPoints;
    CfemInt maxDegree;         // Degré max de polynôme à tester
    CfemReal expectedVol;
    std::string label;
};

class QuadratureAccuracyTest : public ::testing::TestWithParam<QuadratureTestParams> {
protected:
    /**
     * @brief Calcule l'intégrale analytique de x^a * y^b * z^c sur l'élément de référence.
     */
    CfemReal getAnalyticIntegral(CellType type, CfemInt a, CfemInt b, CfemInt c) {
        if (type == CellType::LINE2 || type == CellType::QUAD4 || type == CellType::HEX8) {
            // Intégrale sur [-1, 1]^n : produit de (x^(a+1)/(a+1)) entre -1 et 1
            auto integrate1D = [](CfemInt p) {
                return (p % 2 != 0) ? 0.0 : 2.0 / (p + 1.0);
            };
            CfemReal res = integrate1D(a);
            if (type == CellType::QUAD4 || type == CellType::HEX8) res *= integrate1D(b);
            if (type == CellType::HEX8) res *= integrate1D(c);
            return res;
        } 
        else if (type == CellType::TRI3 || type == CellType::TET4) {
            // Formule magique pour les simplexes sur [0, 1] : a!b!c! / (a+b+c+dim)! * dim! * Vol
            // Pour un triangle (dim 2) : a!b! / (a+b+2)!
            auto fact = [](CfemInt n) {
                CfemReal res = 1.0;
                for (CfemInt i = 2; i <= n; ++i) res *= i;
                return res;
            };
            if (type == CellType::TRI3) return fact(a) * fact(b) / fact(a + b + 2);
            if (type == CellType::TET4) return fact(a) * fact(b) * fact(c) / fact(a + b + c + 3);
        }
        else if (type == CellType::PRISM6 || type == CellType::PRISM18) {
            // Le prisme est un produit Triangle(x,y) [0,1] x Segment(z) [-1,1]
            auto fact = [](CfemInt n) {
                CfemReal res = 1.0;
                for (CfemInt i = 2; i <= n; ++i) res *= i;
                return res;
            };
            
            // 1. Intégrale sur le triangle de référence (x,y)
            CfemReal triPart = fact(a) * fact(b) / fact(a + b + 2);
            
            // 2. Intégrale sur le segment de référence (z) entre -1 et 1
            CfemReal segmentPart = (c % 2 != 0) ? 0.0 : 2.0 / (c + 1.0);
            
            return triPart * segmentPart;
        }
        return 0.0;
    }
};

TEST_P(QuadratureAccuracyTest, PolynomialIntegration) {
    const auto& params = GetParam();
    IntegrationScheme scheme;
    const QuadratureRule& q = scheme.getRule(params.cellType, params.numPoints);

    // 1. Test du Volume (Polynôme degré 0 : f=1)
    CfemReal computedVol = 0.0;
    for (const auto& pt : q.getPoints()) computedVol += pt.weight;
    EXPECT_NEAR(computedVol, params.expectedVol, 1e-13) << "Volume incorrect pour " << params.label;

    // 2. Test des polynômes jusqu'au degré max d'exactitude
    for (CfemInt d = 1; d <= params.maxDegree; ++d) {
        for (CfemInt a = 0; a <= d; ++a) {
            for (CfemInt b = 0; b <= (d - a); ++b) {
                CfemInt c = d - a - b;
                if (params.cellType == CellType::LINE2 && (b > 0 || c > 0)) continue;
                if ((params.cellType == CellType::TRI3 || params.cellType == CellType::QUAD4) && c > 0) continue;

                CfemReal analytic = getAnalyticIntegral(params.cellType, a, b, c);
                CfemReal numeric = 0.0;

                for (const auto& pt : q.getPoints()) {
                    numeric += std::pow(pt.xi.x, a) * std::pow(pt.xi.y, b) * std::pow(pt.xi.z, c) * pt.weight;
                }

                EXPECT_NEAR(numeric, analytic, 1e-12) 
                    << "Erreur intégration monôme x^" << a << " y^" << b << " z^" << c 
                    << " pour " << params.label << " (Degré " << d << ")";
            }
        }
    }
}


class MappingAccuracyTest : public ::testing::Test {
protected:
    // Helper pour intégrer la fonction constante 1.0 sur un élément physique
    CfemReal computePhysicalVolume(const ReferenceElement& ref, 
                                   const QuadratureRule& quad, 
                                   const DynamicVector<Vector3D>& nodes) {
        ShapeInfo shape;
        MappingInfo map;
        shape.setup(ref.getNumNodes());
        map.setup(ref.getNumNodes());

        CfemReal volume = 0.0;
        for (const auto& pt : quad.getPoints()) {
            // L'appel à computeMapping déclenche toute ta chaîne : 
            // evaluate() -> computeJacobian() -> detJ

            ref.computeMapping(
                pt.xi,
                nodes,
                ShapeEvaluationFlags::FirstDerivatives,
                MappingEvaluationFlags::Determinant,
                shape, 
                map
            );
            
            // Intégrale = Somme( |detJ| * w )
            volume += std::abs(map.detJ) * pt.weight;
        }
        return volume;
    }
};

// --- TEST 1 : QUADRILATÈRE DÉFORMÉ ---
TEST_F(MappingAccuracyTest, DeformedQuad4) {
    Q1ReferenceQuad ref;
    IntegrationScheme scheme;
    const QuadratureRule& q = scheme.getRule(CellType::QUAD4, 4); // Gauss 2x2

    // On crée un trapèze en 2D (Z=0)
    // Base large = 2.0, Base courte = 1.0, Hauteur = 1.0
    // Aire analytique = (B+b)*h / 2 = (2+1)*1 / 2 = 1.5
    DynamicVector<Vector3D> nodes(4);
    nodes[0] = {0.0, 0.0, 0.0};
    nodes[1] = {2.0, 0.0, 0.0};
    nodes[2] = {1.5, 1.0, 0.0};
    nodes[3] = {0.5, 1.0, 0.0};

    CfemReal vol = computePhysicalVolume(ref, q, nodes);
    EXPECT_NEAR(vol, 1.5, 1e-13);
}

// --- TEST 2 : TRIANGLE AVEC ROTATION ET TRANSLATION ---
TEST_F(MappingAccuracyTest, TransformedTriangle3) {
    P1ReferenceTriangle ref;
    IntegrationScheme scheme;
    const QuadratureRule& q = scheme.getRule(CellType::TRI3, 3);

    // Triangle rectangle d'aire 0.5, déplacé et tourné
    // Coordonnées : (10,10), (11,10), (10,11) -> Aire = 0.5
    DynamicVector<Vector3D> nodes(3);
    nodes[0] = {10.0, 10.0, 0.0};
    nodes[1] = {11.0, 10.0, 0.0};
    nodes[2] = {10.0, 11.0, 0.0};

    CfemReal vol = computePhysicalVolume(ref, q, nodes);
    EXPECT_NEAR(vol, 0.5, 1e-13);
}

// --- TEST 3 : CONSISTANCE DE L'INTERPOLATION ---
TEST_F(MappingAccuracyTest, InterpolationConsistency) {
    Q1ReferenceQuad ref;
    ShapeInfo shape;
    shape.setup(ref.getNumNodes());
    //shape.setup(ref.getNumNodes());
    
    
    // On définit un carré de côté 2, centré sur (5,5)
    DynamicVector<Vector3D> nodes(4);
    nodes[0] = {4.0, 4.0, 0.0};
    nodes[1] = {6.0, 4.0, 0.0};
    nodes[2] = {6.0, 6.0, 0.0};
    nodes[3] = {4.0, 6.0, 0.0};

    // On teste au centre de l'élément de référence (0,0)
    Vector3D xi_center{0.0, 0.0, 0.0};
    ref.evaluateShapes(xi_center, ShapeEvaluationFlags::Value, shape);

    // L'interpolation des coordonnées des noeuds doit donner le centre physique (5,5)
    Vector3D x_phys = ref.interpolate(shape, nodes);
    
    EXPECT_NEAR(x_phys.x, 5.0, 1e-14);
    EXPECT_NEAR(x_phys.y, 5.0, 1e-14);
}

// --- TEST 4 : GRADIENT PHYSIQUE (Test de la dérivée d'une fonction linéaire) ---
TEST_F(MappingAccuracyTest, PhysicalGradientLinearField) {
    P1ReferenceTriangle ref;
    IntegrationScheme scheme;
    const QuadratureRule& q = scheme.getRule(CellType::TRI3, 1);
    ShapeInfo shape;
    MappingInfo map;
    shape.setup(ref.getNumNodes());
    map.setup(ref.getNumNodes());

    // Triangle : (0,0), (1,0), (0,1)
    DynamicVector<Vector3D> nodes = {{0,0,0}, {1,0,0}, {0,1,0}};
    
    // Champ scalaire linéaire : u(x,y) = 2x + 3y
    // Gradient physique attendu : [2, 3, 0]
    DynamicVector<CfemReal> nodalValues = {0.0, 2.0, 3.0};

    ref.computeMapping(
        q.getPoints()[0].xi,
        nodes,
        ShapeEvaluationFlags::FirstDerivatives,
        MappingEvaluationFlags::PhysicalFirstDerivatives,
        shape,
        map
    );
    Vector3D gradU = ref.interpolateGradient(map, nodalValues);

    EXPECT_NEAR(gradU.x, 2.0, 1e-13);
    EXPECT_NEAR(gradU.y, 3.0, 1e-13);
}

// --- TEST : TRIANGLE INCLINÉ DANS L'ESPACE 3D (Pseudo-Jacobien) ---
TEST_F(MappingAccuracyTest, SlantedTriangle3D) {
    P1ReferenceTriangle ref;
    IntegrationScheme scheme;
    const QuadratureRule& q = scheme.getRule(CellType::TRI3, 3);
    
    // Triangle rectangle incliné à 45° sur l'axe XZ
    // Noeuds : (0,0,0), (cos(45), 0, sin(45)), (0, 1, 0)
    CfemReal s = std::sqrt(2.0) / 2.0;
    DynamicVector<Vector3D> nodes(3);
    nodes[0] = {0.0, 0.0, 0.0};
    nodes[1] = {s,   0.0, s};   // Longueur de l'arête = 1.0
    nodes[2] = {0.0, 1.0, 0.0}; // Longueur de l'arête = 1.0
    
    // L'aire doit rester 0.5 peu importe l'orientation 3D
    CfemReal vol = computePhysicalVolume(ref, q, nodes);
    EXPECT_NEAR(vol, 0.5, 1e-13);
}

// --- TEST : GRADIENT SUR ÉLÉMENT 3D PLANAIRE ---
TEST_F(MappingAccuracyTest, SlantedGradient3D) {
    P1ReferenceTriangle ref;
    // IntegrationScheme scheme;
    // const QuadratureRule& q = scheme.getRule(CellType::TRI3, 1);
    ShapeInfo shape;
    MappingInfo map;
    shape.setup(3);
    map.setup(3);

    // Même triangle incliné
    CfemReal s = std::sqrt(2.0) / 2.0;
    DynamicVector<Vector3D> nodes = {{0,0,0}, {s,0,s}, {0,1,0}};
    
    // Appliquons un champ u = x + z
    // Sur le noeud 1 (s,0,s), u = s + s = sqrt(2)
    DynamicVector<CfemReal> nodalValues = {0.0, std::sqrt(2.0), 0.0};

    ref.computeMapping(
        ref.getReferenceNodes()[1],
        nodes,
        ShapeEvaluationFlags::FirstDerivatives,
        MappingEvaluationFlags::PhysicalFirstDerivatives,
        shape,
        map
    );
    Vector3D gradU = ref.interpolateGradient(map, nodalValues);

    // Le gradient de x+z est [1, 0, 1]
    // Mais attention : l'interpolation sur l'élément ne peut capter 
    // que les variations dans son plan.
    EXPECT_NEAR(gradU.x, 1.0, 1e-13);
    EXPECT_NEAR(gradU.y, 0.0, 1e-13);
    EXPECT_NEAR(gradU.z, 1.0, 1e-13);
}

TEST_F(MappingAccuracyTest, PolynomialIntegrationOnRealElement) {
    P1ReferenceTriangle ref;
    IntegrationScheme scheme;
    // Une règle à 3 points est largement suffisante pour du degré 1 + detJ constant
    const QuadratureRule& q = scheme.getRule(CellType::TRI3, 3);
    
    ShapeInfo shape;
    MappingInfo map;
    shape.setup(ref.getNumNodes());
    map.setup(ref.getNumNodes());

    // 1. Définir un triangle réel arbitraire (Aire = 2.0)
    // Noeuds : (1,1), (3,1), (1,3)
    DynamicVector<Vector3D> nodes = {
        {1.0, 1.0, 0.0}, 
        {3.0, 1.0, 0.0}, 
        {1.0, 3.0, 0.0}
    };
    CfemReal area = 2.0;

    // 2. Fonction à intégrer : f(x, y) = x + y
    auto f = [](const Vector3D& p) { return p.x + p.y; };
    
    // Valeur analytique : Aire * (f(P0) + f(P1) + f(P2)) / 3
    // 2.0 * ((1+1) + (3+1) + (1+3)) / 3 = 2.0 * (2 + 4 + 4) / 3 = 20 / 3 = 6.666666666666667
    CfemReal analytic = 20.0 / 3.0;

    // 3. Intégration numérique
    CfemReal numeric = 0.0;
    for (const auto& pt : q.getPoints()) {
        // Mapping complet au point de quadrature
        ref.computeMapping(
            pt.xi,
            nodes,
            ShapeEvaluationFlags::Value | ShapeEvaluationFlags::FirstDerivatives,
            MappingEvaluationFlags::PhysicalPosition | MappingEvaluationFlags::Determinant, 
            shape,
            map
        );
        
        // On récupère la position physique x(xi) et le poids detJ * w
        numeric += f(map.physicalPosition) * std::abs(map.detJ) * pt.weight;
    }

    EXPECT_NEAR(numeric, analytic, 1e-13);
}


INSTANTIATE_TEST_SUITE_P(
    StandardRules,
    QuadratureAccuracyTest,
    ::testing::Values(
        // Line
        QuadratureTestParams{CellType::LINE2, 2, 3, 2.0, "Line2_2pts"},
        QuadratureTestParams{CellType::LINE2, 3, 5, 2.0, "Line2_3pts"},
        QuadratureTestParams{CellType::LINE2, 4, 7, 2.0, "Line2_4pts"},
        QuadratureTestParams{CellType::LINE2, 5, 9, 2.0, "Line2_5pts"},
        // Triangle
        QuadratureTestParams{CellType::TRI3, 3, 2, 0.5, "Triangle3_3pts"},
        QuadratureTestParams{CellType::TRI3, 4, 3, 0.5, "Triangle3_4pts"},
        QuadratureTestParams{CellType::TRI3, 6, 4, 0.5, "Triangle3_6pts"},
        QuadratureTestParams{CellType::TRI3, 7, 5, 0.5, "Triangle3_7pts"},
        QuadratureTestParams{CellType::TRI3, 12, 6, 0.5, "Triangle3_12pts"},
        QuadratureTestParams{CellType::TRI3, 16, 8, 0.5, "Triangle3_16pts"},
        QuadratureTestParams{CellType::TRI3, 19, 9, 0.5, "Triangle3_19pts"},
        QuadratureTestParams{CellType::TRI3, 25, 10, 0.5, "Triangle3_25pts"},
        // Quad 
        QuadratureTestParams{CellType::QUAD4, 4, 3, 4.0, "Quad4_4pts"},
        QuadratureTestParams{CellType::QUAD4, 9, 5, 4.0, "Quad4_9pts"},
        QuadratureTestParams{CellType::QUAD4, 16, 7, 4.0, "Quad4_16pts"},
        QuadratureTestParams{CellType::QUAD4, 25, 9, 4.0, "Quad4_25pts"},
        // Tetra
        QuadratureTestParams{CellType::TET4, 4, 2, 1.0/6.0, "Tetra4_4pts"},
        QuadratureTestParams{CellType::TET4, 5, 3, 1.0/6.0, "Tetra4_5pts"},
        QuadratureTestParams{CellType::TET4, 11, 4, 1.0/6.0, "Tetra4_11pts"},
        QuadratureTestParams{CellType::TET4, 15, 5, 1.0/6.0, "Tetra4_15pts"},
        QuadratureTestParams{CellType::TET4, 29, 6, 1.0/6.0, "Tetra4_29pts"},
        QuadratureTestParams{CellType::TET4, 33, 7, 1.0/6.0, "Tetra4_33pts"},
        // Prism
        QuadratureTestParams{CellType::PRISM6, 6, 2, 1.0, "Prism6_6pts"},
        QuadratureTestParams{CellType::PRISM6, 8, 3, 1.0, "Prism6_8pts"},
        QuadratureTestParams{CellType::PRISM6, 18, 4, 1.0, "Prism6_18pts"},
        QuadratureTestParams{CellType::PRISM6, 21, 5, 1.0, "Prism6_21pts"},
        QuadratureTestParams{CellType::PRISM6, 48, 6, 1.0, "Prism6_48pts"},
        QuadratureTestParams{CellType::PRISM6, 80, 8, 1.0, "Prism6_80pts"},
        QuadratureTestParams{CellType::PRISM6, 95, 9, 1.0, "Prism6_95pts"},
        // Hex 
        QuadratureTestParams{CellType::HEX8, 8, 3, 8.0, "Hex8_8pts"},
        QuadratureTestParams{CellType::HEX8, 27, 5, 8.0, "Hex8_27pts"},
        QuadratureTestParams{CellType::HEX8, 64, 7, 8.0, "Hex8_64pts"},
        QuadratureTestParams{CellType::HEX8, 125, 9, 8.0, "Hex8_125pts"}
    ),
    [](const ::testing::TestParamInfo<QuadratureTestParams>& info) { return info.param.label; }
);