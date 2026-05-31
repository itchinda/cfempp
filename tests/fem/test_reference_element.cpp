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

#include "ReferenceElements.h"

using namespace cfem;

/**
 * @brief Fixture pour tester les éléments de référence.
 * On utilise un template ou une interface pour tester tous les éléments.
 */
class ReferenceElementTest : public ::testing::TestWithParam<ReferenceElement*> {
protected:
    ShapeInfo info;
    const CfemReal tol = cfem::ZERO_TOLERANCE;

    DynamicVector<Vector3D> getTestPointsForDim(CfemInt dim){
        DynamicVector<Vector3D> test_points;
        if (dim == 1) {
            test_points.push_back(Vector3D{0.5, 0.0, 0.0}); // Milieu
            test_points.push_back(Vector3D{0.25, 0.0, 0.0});
        } else if (dim == 2) {
            test_points.push_back(Vector3D{1.0/3.0, 1.0/3.0, 0.0}); // Barycentre Triangle
            test_points.push_back(Vector3D{0.2, 0.2, 0.0});
            test_points.push_back(Vector3D{0.1, 0.8, 0.0}); // Près d'un sommet
        } else if (dim == 3) {
            test_points.push_back(Vector3D{0.25, 0.25, 0.25}); // Barycentre Tétraèdre
        } else{
            CFEM_ERROR << "Dim must be within 1 and 3 inclusively";
        }
        return test_points;
    }

    // Petite fonction helper pour accéder proprement aux composantes de Vector3D
    CfemReal getComp(const Vector3D& v, CfemInt dim) {
        if (dim == 0) return v.x;
        if (dim == 1) return v.y;
        return v.z;
    }

    void SetUp() override {
        // Optionnel : configuration avant chaque test
    }



};

/**
 * @brief Fixture étendue pour le Mapping.
 */
class MappingTest : public ::testing::Test {
protected:
    const CfemReal tol = 1e-11;

    // Helper pour un Quad4 parfait dans le plan Z=0
    DynamicVector<Vector3D> createUnitQuad() {
        DynamicVector<Vector3D> coords(4);
        coords[0] = {-1, -1, 0}; coords[1] = { 1, -1, 0};
        coords[2] = { 1,  1, 0}; coords[3] = {-1,  1, 0};
        return coords;
    }

    // Helper pour un Hex8 parfait
    DynamicVector<Vector3D> createUnitHex() {
        DynamicVector<Vector3D> coords(8);
        coords[0] = {-1, -1, -1}; coords[1] = { 1, -1, -1};
        coords[2] = { 1,  1, -1}; coords[3] = {-1,  1, -1};
        coords[4] = {-1, -1,  1}; coords[5] = { 1, -1,  1};
        coords[6] = { 1,  1,  1}; coords[7] = {-1,  1,  1};
        return coords;
    }
};


// TEST DE KRONECKER : N_i(node_j) = delta_ij ---
TEST_P(ReferenceElementTest, KroneckerDeltaProperty) {
    ReferenceElement* el = GetParam();
    const auto& nodes = el->getReferenceNodes();
    CfemInt n = el->getNumNodes();

    for (CfemInt j = 0; j < n; ++j) {
        el->evaluateShapes(nodes[j], ShapeEvaluationFlags::Value, info);
        const auto& values = info.values();
        for (CfemInt i = 0; i < n; ++i) {
            CfemReal expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(values[i], expected, tol) 
                << "Element: " << typeid(*el).name() << " - Node: " << j << " - ShapeFn: " << i;
        }
    }
}

TEST_P(ReferenceElementTest, PartitionOfUnity) {
    ReferenceElement* el = GetParam();
    CfemInt dim = el->getDimension();
    
    DynamicVector<Vector3D> test_points = getTestPointsForDim(dim);

    for (const auto& xi : test_points) {
        el->evaluateShapes(xi, ShapeEvaluationFlags::Value, info);
        const auto& values = info.values();

        CfemReal sum = 0.0;
        for (CfemInt i = 0; i < el->getNumNodes(); ++i) {
            sum += values[i];
        }
        
        EXPECT_NEAR(sum, 1.0, tol) 
            << "Failed Partition of Unity for " << typeid(*el).name() 
            << " at xi(" << xi.x << ", " << xi.y << ", " << xi.z << ")";
    }
}

// TEST DE LA SOMME DES GRADIENTS : Sum(Grad N_i) = 0 ---
TEST_P(ReferenceElementTest, GradientSumIsZero) {
    ReferenceElement* el = GetParam();
    CfemInt dim {el->getDimension()};

    DynamicVector<Vector3D> test_points = getTestPointsForDim(dim);

    for (const auto& xi : test_points){

        el->evaluateShapes(xi, ShapeEvaluationFlags::FirstDerivatives, info);
        
        Vector3D gradSum{0.0, 0.0, 0.0};
        for (CfemInt i = 0; i < el->getNumNodes(); ++i) {
            gradSum += info.gradient(i);
        }

        EXPECT_NEAR(gradSum.x, 0.0, tol);
        EXPECT_NEAR(gradSum.y, 0.0, tol);
        EXPECT_NEAR(gradSum.z, 0.0, tol);
    }
}

TEST_P(ReferenceElementTest, HessianFullConsistency) {
    const CfemReal eps = 1e-6;
    const CfemReal h_tol = 1e-7;

    ReferenceElement* el = GetParam();
    const CfemInt n = el->getNumNodes();
    const CfemInt dim = el->getDimension();
    
    // Réutilisation de tes test_points robustes
    DynamicVector<Vector3D> test_points = getTestPointsForDim(dim); 

    ShapeInfo info, infoPlus, infoMinus;

    for (const auto& xi : test_points) {
        el->evaluateShapes(xi, ShapeEvaluationFlags::FirstDerivatives | ShapeEvaluationFlags::SecondDerivatives, info);

        for (CfemInt i = 0; i < n; ++i) {
            for (CfemInt r = 0; r < dim; ++r) {
                for (CfemInt c = 0; c < dim; ++c) {
                    
                    Vector3D xiPlus = xi;
                    Vector3D xiMinus = xi;
                    
                    if (c == 0) { xiPlus.x += eps; xiMinus.x -= eps; }
                    else if (c == 1) { xiPlus.y += eps; xiMinus.y -= eps; }
                    else if (c == 2) { xiPlus.z += eps; xiMinus.z -= eps; }

                    el->evaluateShapes(xiPlus, ShapeEvaluationFlags::FirstDerivatives, infoPlus);
                    el->evaluateShapes(xiMinus, ShapeEvaluationFlags::FirstDerivatives, infoMinus);

                    // Valeur du gradient au nœud i pour la composante r
                    CfemReal grad_r_plus = getComp(infoPlus.gradient(i), r);
                    CfemReal grad_r_minus = getComp(infoMinus.gradient(i), r);
                    
                    CfemReal d2N_drc_fd = (grad_r_plus - grad_r_minus) / (2.0 * eps);

                    EXPECT_NEAR(info.hessian(i)(r, c), d2N_drc_fd, h_tol)
                        << "Hessian mismatch at node " << i 
                        << " component (" << r << "," << c << ")"
                        << " for element " << typeid(*el).name();
                }
            }
        }
    }
}


// TEST : ÉLÉMENT 2D DANS L'ESPACE 3D ---
TEST_F(MappingTest, QuadPlanarMapping) {
    Q1ReferenceQuad quad;
    auto coords = createUnitQuad();
    MappingInfo map;
    ShapeInfo shape;
    Vector3D xi{0.5, 0.5, 0}; // Un point dans l'élément

    // On vérifie que le mapping ne crash pas malgré le déterminant 3D nul (Z=0)
    
    // Si on simule une épaisseur ou un élément 3D "aplati" :
    for(CfemInt i=0; i<4; ++i) coords[i].z = 0.0; 
    
    // Ici on teste surtout la cohérence des gradients dN/dx et dN/dy
    
    ShapeEvaluationFlags shapeFlags = ShapeEvaluationFlags::FirstDerivatives;
    MappingEvaluationFlags mappingFlags = MappingEvaluationFlags::PhysicalFirstDerivatives;
    quad.computeMapping(
        xi,
        coords,
        shapeFlags,
        mappingFlags,
        shape,
        map
    );
    
    // dN/dz doit être nul car l'élément est dans le plan XY
    for(CfemInt i=0; i<4; ++i) {
        EXPECT_NEAR(map.dN_dz()[i], 0.0, tol);
    }
}

// TEST : DISTORSION AFFINE (Cisaillement) ---
// On teste si J^-T transforme correctement les angles
TEST_F(MappingTest, Hex8ShearMapping) {
    Q1ReferenceHexahedron hex;
    auto coords = createUnitHex();
    
    // Appliquer un cisaillement : x' = x + 0.5*y
    for(auto& p : coords) p.x += 0.5 * p.y;

    MappingInfo map;
    ShapeInfo shape;
    Vector3D xi{0, 0, 0};
    ShapeEvaluationFlags shapeFlags = ShapeEvaluationFlags::FirstDerivatives;
    MappingEvaluationFlags mappingFlags = MappingEvaluationFlags::Jacobian | MappingEvaluationFlags::Determinant;

    hex.computeMapping(
        xi,
        coords,
        shapeFlags,
        mappingFlags,
        shape,
        map
    );

    // Le déterminant doit rester 1.0 (le cisaillement préserve le volume)
    EXPECT_NEAR(map.detJ, 1.0, tol);

    // Vérifier que le Jacobien reflète le cisaillement
    EXPECT_NEAR(map.jacobian(0,0), 1.0, tol);
    EXPECT_NEAR(map.jacobian(0,1), 0.5, tol); // shear: dx/deta = 0.5 (car x dépend de y, et y dépend de eta)
    EXPECT_NEAR(map.jacobian(0,2), 0.0, tol);

    EXPECT_NEAR(map.jacobian(1,0), 0.0, tol);
    EXPECT_NEAR(map.jacobian(1,1), 1.0, tol);
    EXPECT_NEAR(map.jacobian(1,2), 0.0, tol);

    EXPECT_NEAR(map.jacobian(2,0), 0.0, tol);
    EXPECT_NEAR(map.jacobian(2,1), 0.0, tol);
    EXPECT_NEAR(map.jacobian(2,2), 1.0, tol);
}

TEST_F(MappingTest, P2TriangleCurvedHessians) {
    P2ReferenceTriangle tri;
    DynamicVector<Vector3D> coords(6);
    coords[0]={0,0,0}; coords[1]={1,0,0}; coords[2]={0,1,0};
    coords[3]={0.5, 0, 0}; coords[4]={0.5, 0.5, 0}; coords[5]={0, 0.5, 0};
    
    // Courber l'élément : déplacer le noeud milieu de l'hypoténuse
    coords[4].x += 0.2; 
    coords[4].y += 0.2;

    MappingInfo map;
    ShapeInfo shape;
    Vector3D xi{1.0/3.0, 1.0/3.0, 0};
    
    // On demande explicitement les Hessiennes et gradients
    ShapeEvaluationFlags shapeFlags = ShapeEvaluationFlags::FirstDerivatives | ShapeEvaluationFlags::SecondDerivatives;
    MappingEvaluationFlags mappingFlags = MappingEvaluationFlags::PhysicalFirstDerivatives | MappingEvaluationFlags::PhysicalSecondDerivatives;

    tri.computeMapping(
        xi,
        coords, 
        shapeFlags,
        mappingFlags,
        shape,
        map
    );

    EXPECT_GT(map.detJ, 0.0) << "Jacobian determinant should be positive";

    bool hasNonZeroHessian = false;
    for (CfemInt i=0; i<6; ++i) {
        const auto& H = map.hessian(i);

        if (std::abs(H.xx()) > 1e-5 || std::abs(H.yy()) > 1e-5 || std::abs(H.zz()) > 1e-5 ||
            std::abs(H.xy()) > 1e-5 || std::abs(H.xz()) > 1e-5 || std::abs(H.yz()) > 1e-5) {
            hasNonZeroHessian = true;
            break;
        }
    }
    EXPECT_TRUE(hasNonZeroHessian) << "Curvature should generate non-zero physical Hessians";

    // Check that the gradient components are reasonable (not NaN)
    for (CfemInt i=0; i<6; ++i) {
        EXPECT_FALSE(std::isnan(map.dN_dx()[i]));
        EXPECT_FALSE(std::isnan(map.dN_dy()[i]));
        EXPECT_FALSE(std::isnan(map.dN_dz()[i]));
    }
}

// TEST : PHYSICAL GRADIENT SUM = 0 (consistency check in physical space) ---
TEST_F(MappingTest, PhysicalGradientSumIsZero) {
    Q1ReferenceQuad quad;
    auto coords = createUnitQuad();
    MappingInfo map;
    ShapeInfo shape;
    Vector3D xi{0.25, 0.25, 0.0};

    ShapeEvaluationFlags shapeFlags = ShapeEvaluationFlags::FirstDerivatives;
    MappingEvaluationFlags mappingFlags = MappingEvaluationFlags::PhysicalFirstDerivatives;

    quad.computeMapping(
        xi,
        coords,
        shapeFlags,
        mappingFlags,
        shape,
        map
    );

    Vector3D sumGrad{0.0, 0.0, 0.0};
    for (CfemInt i=0; i<4; ++i) sumGrad += map.gradient(i);

    EXPECT_NEAR(sumGrad.x, 0.0, tol);
    EXPECT_NEAR(sumGrad.y, 0.0, tol);
    EXPECT_NEAR(sumGrad.z, 0.0, tol);
}

// J * J_inv ~= Identity (Jacobian inverse correctness) ---
TEST_F(MappingTest, JacobianInverseCheck) {
    Q1ReferenceHexahedron hex;
    auto coords = createUnitHex();
    MappingInfo map;
    ShapeInfo shape;
    Vector3D xi{0.0, 0.0, 0.0};

    ShapeEvaluationFlags shapeFlags = ShapeEvaluationFlags::FirstDerivatives;
    MappingEvaluationFlags mappingFlags = MappingEvaluationFlags::InverseJacobian;

    hex.computeMapping(
        xi,
        coords,
        shapeFlags,
        mappingFlags,
        shape,
        map
    );

    Tensor3D identity;
    identity.setZero();
    for (CfemInt i=0; i<3; ++i) identity(i,i) = 1.0;

    Tensor3D check = map.jacobian * map.invJacobian;

    for (CfemInt i=0; i<3; ++i)
        for (CfemInt j=0; j<3; ++j)
            EXPECT_NEAR(check(i,j), identity(i,j), tol);
}

/**
 * @test Verify that physical Hessians for a curved P2 Tetrahedron are correct.
 * For a P2 element, if we map the element to a specific quadratic transformation,
 * the physical Hessians should remain consistent with the chain rule.
 */
TEST_F(MappingTest, P2TetraCurvedHessians) {
    P2ReferenceTetrahedron tet;
    const CfemInt n = tet.getNumNodes();
    DynamicVector<Vector3D> coords(n);

    // Define a curved physical geometry
    const auto& refNodes = tet.getReferenceNodes();
    for(CfemInt i = 0; i < n; ++i) {
        coords[i] = refNodes[i];
        coords[i].x += 0.2 * refNodes[i].x * refNodes[i].x;
        coords[i].y += 0.1 * refNodes[i].x * refNodes[i].y;
    }

    MappingInfo map;
    ShapeInfo shape;
    Vector3D xi{0.25, 0.25, 0.25}; 
    
    ShapeEvaluationFlags shapeFlags = ShapeEvaluationFlags::FirstDerivatives | ShapeEvaluationFlags::SecondDerivatives;
    MappingEvaluationFlags mappingFlags = MappingEvaluationFlags::PhysicalFirstDerivatives | MappingEvaluationFlags::PhysicalSecondDerivatives;

    ASSERT_NO_THROW(tet.computeMapping(xi, coords, shapeFlags, mappingFlags, shape, map));
    
    // SANITY CHECK
    EXPECT_GT(map.detJ, 0.0) << "Curved element became inverted!";
    
    bool hasNonZeroHessian = false;
    for (CfemInt i = 0; i < n; ++i) {
        const auto& H = map.hessian(i);
        if (std::abs(H.xx())>1e-5 || std::abs(H.yy())>1e-5 || std::abs(H.zz())>1e-5 ||
            std::abs(H.xy())>1e-5 || std::abs(H.xz())>1e-5 || std::abs(H.yz())>1e-5) {
            hasNonZeroHessian = true;
            break;
        }
    }
    ASSERT_TRUE(hasNonZeroHessian) << "All Hessians are zero! Mathematical checks below will falsely pass.";

    // MATHEMATICAL CORRECTNESS
    // The Partition of Unity: Sum of Hessians = 0
    SymTensor3D hessianSum;
    hessianSum.setZero();
    for (CfemInt i = 0; i < n; ++i) {
        hessianSum += map.hessian(i);
    }
    for (CfemInt i = 0; i < 6; ++i) {
        EXPECT_NEAR(hessianSum[i], 0.0, 1e-12)
            << "Physical Hessian sum (Partition of Unity) failed at component " << i;
    }

    // Linear Field Reproduction: u = x + y + z => H(u) = 0
    SymTensor3D fieldHessian;
    fieldHessian.setZero();
    for (CfemInt i = 0; i < n; ++i) {
        fieldHessian += map.hessian(i) * (coords[i].x + coords[i].y + coords[i].z);
    }
    for (CfemInt i = 0; i < 6; ++i) {
        EXPECT_NEAR(fieldHessian[i], 0.0, 1e-11) 
            << "Curvature cancellation failed! Component " << i 
            << " did not account for geometric Hessian correctly.";
    }
}

TEST_F(MappingTest, P2TetraQuadraticFieldReproduction) {
    P2ReferenceTetrahedron tet;
    const CfemInt n = tet.getNumNodes();
    DynamicVector<Vector3D> coords(n);

    const auto& refNodes = tet.getReferenceNodes();

    // Create a distorted element using a linear (affine) transformation.
    // This ensures a constant Jacobian, allowing the P2 element to 
    // perfectly reproduce quadratic fields without mapping-induced rational errors.
    for (CfemInt i = 0; i < n; ++i) {
        coords[i].x = 1.1 * refNodes[i].x + 0.2 * refNodes[i].y;
        coords[i].y = 0.9 * refNodes[i].y - 0.1 * refNodes[i].z;
        coords[i].z = 1.05 * refNodes[i].z + 0.05 * refNodes[i].x;
    }

    MappingInfo map;
    ShapeInfo shape;
    Vector3D xi{0.2, 0.3, 0.1}; 

    
    ShapeEvaluationFlags shapeFlags = ShapeEvaluationFlags::FirstDerivatives | ShapeEvaluationFlags::SecondDerivatives;
    MappingEvaluationFlags mappingFlags = MappingEvaluationFlags::PhysicalFirstDerivatives | MappingEvaluationFlags::PhysicalSecondDerivatives;

    tet.computeMapping(xi, coords, shapeFlags, mappingFlags, shape, map);

    // Quadratic Field: u = x^2 + y^2 + z^2 + xy + xz + yz
    auto fieldFunc = [](const Vector3D& p) {
        return (p.x*p.x + p.y*p.y + p.z*p.z + p.x*p.y + p.x*p.z + p.y*p.z);
    };

    // Accumulate the physical Hessian of the interpolated field
    SymTensor3D interpInHessian;
    interpInHessian.setZero();
    for (CfemInt i = 0; i < n; ++i) {
        interpInHessian += map.hessian(i) * fieldFunc(coords[i]);
    }

    // 4. VALIDATION: Compare against the analytical Hessian
    // Components: xx=2, yy=2, zz=2, xy=1, xz=1, yz=1
    EXPECT_NEAR(interpInHessian.xx(), 2.0, 1e-10);
    EXPECT_NEAR(interpInHessian.yy(), 2.0, 1e-10);
    EXPECT_NEAR(interpInHessian.zz(), 2.0, 1e-10);
    EXPECT_NEAR(interpInHessian.xy(), 1.0, 1e-10);
    EXPECT_NEAR(interpInHessian.xz(), 1.0, 1e-10);
    EXPECT_NEAR(interpInHessian.yz(), 1.0, 1e-10);
}
// --- 8. TEST : DEGENERATE / NEAR-DEGENERATE TETRAHEDRON ---
TEST_F(MappingTest, DegenerateTetrahedron) {
#ifndef NDEBUG
    P1ReferenceTetrahedron tet;
    DynamicVector<Vector3D> coords(4);

    // Standard tetrahedron
    coords[0]={0,0,0}; coords[1]={1,0,0}; coords[2]={0,1,0}; coords[3]={0,0,1};

    MappingInfo map;
    ShapeInfo shape;
    Vector3D xi{0.25,0.25,0.25};
    
    ShapeEvaluationFlags shapeFlags = ShapeEvaluationFlags::FirstDerivatives;
    MappingEvaluationFlags mappingFlags = MappingEvaluationFlags::InverseJacobian;

    // Collapse two nodes -> degenerate
    coords[3] = coords[2];
    EXPECT_THROW({
        tet.computeMapping(xi, coords, shapeFlags, mappingFlags, shape, map);
    }, AssertionException);
#endif
}

TEST_F(MappingTest, P2TetraInterpolationSuite) {
    P2ReferenceTetrahedron tet;
    const CfemInt n = tet.getNumNodes();
    
    // 1. Setup Physical Coordinates (Distorted but Affine for exactness)
    DynamicVector<Vector3D> coords(n);
    const auto& refNodes = tet.getReferenceNodes();
    for(CfemInt i = 0; i < n; ++i) {
        coords[i].x = 1.2 * refNodes[i].x + 0.1 * refNodes[i].y;
        coords[i].y = 0.9 * refNodes[i].y;
        coords[i].z = 1.1 * refNodes[i].z + 0.05 * refNodes[i].x;
    }

    // 2. Define Nodal Values for Scalar Field (u = x + y) and Vector Field (v = [x, y, z])
    DynamicVector<CfemReal> scalarNodalValues(n);
    DynamicVector<Vector3D> vectorNodalValues(n);
    for(CfemInt i = 0; i < n; ++i) {
        scalarNodalValues[i] = coords[i].x + coords[i].y;
        vectorNodalValues[i] = coords[i]; // Identity field: v(x) = x
    }

    // 3. Evaluate Shape Functions at an arbitrary point
    Vector3D xi{0.25, 0.25, 0.25}; // Center-ish of the tet
    ShapeInfo shape;
    MappingInfo map;
    // We need Gradients for the gradient tests later
    
    ShapeEvaluationFlags shapeFlags = ShapeEvaluationFlags::Value | ShapeEvaluationFlags::FirstDerivatives;
    MappingEvaluationFlags mappingFlags = MappingEvaluationFlags::PhysicalPosition | MappingEvaluationFlags::PhysicalFirstDerivatives;

    tet.computeMapping(xi, coords, shapeFlags, mappingFlags, shape, map);

    // --- TEST A: Scalar Value Interpolation ---
    CfemReal expectedScalar = map.physicalPosition.x + map.physicalPosition.y;
    CfemReal interpolatedScalar = tet.interpolate(shape, scalarNodalValues);
    EXPECT_NEAR(interpolatedScalar, expectedScalar, 1e-12);

    // --- TEST B: Vector Value Interpolation ---
    // For the identity field v(x)=x, the interpolated value should be the physical position
    Vector3D expectedVector = map.physicalPosition;
    Vector3D interpolatedVector = tet.interpolate(shape, vectorNodalValues);
    EXPECT_NEAR(interpolatedVector.x, expectedVector.x, 1e-12);
    EXPECT_NEAR(interpolatedVector.y, expectedVector.y, 1e-12);
    EXPECT_NEAR(interpolatedVector.z, expectedVector.z, 1e-12);

    // --- TEST C: Scalar Gradient Interpolation ---
    // For u = x + y, grad(u) = [1, 1, 0]
    Vector3D interpolatedGrad = tet.interpolateGradient(map, scalarNodalValues);
    EXPECT_NEAR(interpolatedGrad.x, 1.0, 1e-12);
    EXPECT_NEAR(interpolatedGrad.y, 1.0, 1e-12);
    EXPECT_NEAR(interpolatedGrad.z, 0.0, 1e-12);
}



// --- INSTANCIATION DES TESTS ---
//1D elements
static P0ReferenceSegment  p0seg;
static P1ReferenceSegment  p1seg;
static P2ReferenceSegment  p2seg;

static P0ReferenceTriangle   p0tri;
static P1ReferenceTriangle   p1tri;
static P2ReferenceTriangle   p2tri;
static MiniReferenceTriangle minitri;

static Q0ReferenceQuad       q0quad;
static Q1ReferenceQuad       q1quad;
static Q2ReferenceQuad       q2quad;
static MiniReferenceQuad     miniquad;

// 3D elements
static P0ReferenceTetrahedron   p0tetra;
static P1ReferenceTetrahedron   p1tetra;
static P2ReferenceTetrahedron   p2tetra;
static MiniReferenceTetrahedron minitetra;

static P0ReferencePrism       p0prism;
static P1ReferencePrism       p1prism;
static P2ReferencePrism       p2prism;
static MiniReferencePrism     miniprism;

static Q0ReferenceHexahedron   q0hex;
static Q1ReferenceHexahedron   q1hex;
static Q2ReferenceHexahedron   q2hex;
static MiniReferenceHexahedron minihex;

INSTANTIATE_TEST_SUITE_P(
    AllElements,
    ReferenceElementTest,
    ::testing::Values(
                    &p0seg,   &p1seg,   &p2seg,  
                    &p0tri,   &p1tri,   &p2tri,   &minitri,
                    &q0quad,  &q1quad,  &q2quad,  &miniquad,
                    &p0tetra, &p1tetra, &p2tetra, &minitetra,
                    &p0prism, &p1prism, &p2prism, &miniprism,
                    &q0hex,   &q1hex,   &q2hex,   &minihex
                )
    // ::testing::Values(
    //                 &p0seg,   &p1seg,   &p2seg,
    //                 &p0tri,   &p1tri,  &minitri
    //             )
);