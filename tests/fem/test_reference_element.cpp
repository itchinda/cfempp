#include <gtest/gtest.h>
#include <cmath>

#include "ReferenceTriangle.h"
#include "ReferenceSegment.h"
#include "ReferenceQuad.h"

#include "ReferenceTetrahedron.h"
#include "ReferenceHexahedron.h"
#include "ReferencePrism.h"

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
            CFEM_ERROR("Dim must be within 1 and 3 inclusively");
        }
        return test_points;
    }

    // Petite fonction helper pour accéder proprement aux composantes de Vector3D
    CfemReal getComp(const Vector3D& v, int dim) {
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


// --- 1. TEST DE KRONECKER : N_i(node_j) = delta_ij ---
TEST_P(ReferenceElementTest, KroneckerDeltaProperty) {
    ReferenceElement* el = GetParam();
    const auto& nodes = el->referenceNodes();
    int n = el->numNodes();

    for (CfemInt j = 0; j < n; ++j) {
        el->evaluate(nodes[j], ShapeActiveField::Value, info);
        for (CfemInt i = 0; i < n; ++i) {
            CfemReal expected = (i == j) ? 1.0 : 0.0;
            EXPECT_NEAR(info.values[i], expected, tol) 
                << "Element: " << typeid(*el).name() << " - Node: " << j << " - ShapeFn: " << i;
        }
    }
}

TEST_P(ReferenceElementTest, PartitionOfUnity) {
    ReferenceElement* el = GetParam();
    int dim = el->dimension();
    
    DynamicVector<Vector3D> test_points = getTestPointsForDim(dim);

    for (const auto& xi : test_points) {
        el->evaluate(xi, ShapeActiveField::Value, info);
        
        CfemReal sum = 0.0;
        for (CfemInt i = 0; i < el->numNodes(); ++i) {
            sum += info.values[i];
        }
        
        EXPECT_NEAR(sum, 1.0, tol) 
            << "Failed Partition of Unity for " << typeid(*el).name() 
            << " at xi(" << xi.x << ", " << xi.y << ", " << xi.z << ")";
    }
}

// --- 3. TEST DE LA SOMME DES GRADIENTS : Sum(Grad N_i) = 0 ---
TEST_P(ReferenceElementTest, GradientSumIsZero) {
    ReferenceElement* el = GetParam();
    CfemInt dim {el->dimension()};

    DynamicVector<Vector3D> test_points = getTestPointsForDim(dim);

    for (const auto& xi : test_points){

        el->evaluate(xi, ShapeActiveField::Gradient, info);
        
        Vector3D gradSum{0.0, 0.0, 0.0};
        for (CfemInt i = 0; i < el->numNodes(); ++i) {
            gradSum.x += info.gradients[i].x;
            gradSum.y += info.gradients[i].y;
            gradSum.z += info.gradients[i].z;
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
    const int n = el->numNodes();
    const int dim = el->dimension();
    
    // Réutilisation de tes test_points robustes
    DynamicVector<Vector3D> test_points = getTestPointsForDim(dim); 

    ShapeInfo info, infoPlus, infoMinus;

    for (const auto& xi : test_points) {
        el->evaluate(xi, ShapeActiveField::Hessian | ShapeActiveField::Gradient, info);

        for (CfemInt i = 0; i < n; ++i) {
            for (CfemInt r = 0; r < dim; ++r) {
                for (CfemInt c = 0; c < dim; ++c) {
                    
                    Vector3D xiPlus = xi;
                    Vector3D xiMinus = xi;
                    
                    if (c == 0) { xiPlus.x += eps; xiMinus.x -= eps; }
                    else if (c == 1) { xiPlus.y += eps; xiMinus.y -= eps; }
                    else if (c == 2) { xiPlus.z += eps; xiMinus.z -= eps; }

                    el->evaluate(xiPlus, ShapeActiveField::Gradient, infoPlus);
                    el->evaluate(xiMinus, ShapeActiveField::Gradient, infoMinus);

                    // Valeur du gradient au nœud i pour la composante r
                    CfemReal grad_r_plus = getComp(infoPlus.gradients[i], r);
                    CfemReal grad_r_minus = getComp(infoMinus.gradients[i], r);
                    
                    CfemReal d2N_drc_fd = (grad_r_plus - grad_r_minus) / (2.0 * eps);

                    EXPECT_NEAR(info.hessians[i](r, c), d2N_drc_fd, h_tol)
                        << "Hessian mismatch at node " << i 
                        << " component (" << r << "," << c << ")"
                        << " for element " << typeid(*el).name();
                }
            }
        }
    }
}


// --- 1. TEST : ÉLÉMENT 2D DANS L'ESPACE 3D ---
TEST_F(MappingTest, QuadPlanarMapping) {
    Q1ReferenceQuad quad;
    auto coords = createUnitQuad();
    MappingInfo map;
    Vector3D xi{0.5, 0.5, 0}; // Un point dans l'élément

    // On vérifie que le mapping ne crash pas malgré le déterminant 3D nul (Z=0)
    // NOTE: Si ton inverse() est purement 3x3, il faudra s'assurer que J(2,2) = 1 
    // ou traiter le cas 2D séparément.
    
    // Si on simule une épaisseur ou un élément 3D "aplati" :
    for(int i=0; i<4; ++i) coords[i].z = 0.0; 
    
    // Ici on teste surtout la cohérence des gradients dN/dx et dN/dy
    
    CfemIntFlag lRequestedFields = ShapeActiveField::Gradient;
    quad.computeMapping(xi, coords, lRequestedFields, map);
    
    // dN/dz doit être nul car l'élément est dans le plan XY
    for(int i=0; i<4; ++i) {
        EXPECT_NEAR(map.dNdz(i), 0.0, tol);
    }
}

// --- 2. TEST : DISTORSION AFFINE (Cisaillement) ---
// On teste si J^-T transforme correctement les angles
TEST_F(MappingTest, Hex8ShearMapping) {
    Q1ReferenceHexahedron hex;
    auto coords = createUnitHex();
    
    // Appliquer un cisaillement : x' = x + 0.5*y
    for(auto& p : coords) p.x += 0.5 * p.y;

    MappingInfo map;
    Vector3D xi{0, 0, 0};
    CfemIntFlag lRequestedFields = ShapeActiveField::Gradient;

    hex.computeMapping(xi, coords, lRequestedFields, map);

    // Le déterminant doit rester 1.0 (le cisaillement préserve le volume)
    EXPECT_NEAR(map.detJ, 1.0, tol);

    // Vérifier que le Jacobien reflète le cisaillement
    // dx/deta = 0.5 (car x dépend de y, et y dépend de eta)
    EXPECT_NEAR(map.jacobian(0,1), 0.5, tol);
}

// --- 3. TEST : COURBURE ET HESSIENNES (Le plus complet) ---
// On utilise un élément quadratique et on déplace un noeud pour créer une courbure
TEST_F(MappingTest, P2TriangleCurvedHessians) {
    P2ReferenceTriangle tri; // Élément à 6 noeuds
    DynamicVector<Vector3D> coords(6);

    // Triangle rectangle isocèle
    coords[0]={0,0,0}; coords[1]={1,0,0}; coords[2]={0,1,0};
    // Noeuds milieux droits
    coords[3]={0.5, 0, 0}; coords[4]={0.5, 0.5, 0}; coords[5]={0, 0.5, 0};
    
    // ON COURBE L'ÉLÉMENT : on déplace le noeud milieu de l'hypoténuse
    coords[4].x += 0.1; 
    coords[4].y += 0.1;

    MappingInfo map;
    Vector3D xi{1.0/3.0, 1.0/3.0, 0}; // Barycentre
    
    // On demande explicitement les Hessiennes
    
    CfemIntFlag lRequestedFields = ShapeActiveField::Gradient;

    lRequestedFields |= ShapeActiveField::Hessian;

    tri.computeMapping(xi, coords, lRequestedFields, map);

    // 1. Le déterminant ne doit pas être nul ou délirant
    EXPECT_GT(map.detJ, 0.0);

    // 2. Les Hessiennes physiques ne doivent PAS être nulles 
    // à cause de la courbure géométrique induite par le noeud déplacé.
    bool hasNonZeroHessian = false;
    for(int i=0; i<6; ++i) {
        if(std::abs(map.d2Ndxdx(i)) > 1e-5) hasNonZeroHessian = true;
    }
    EXPECT_TRUE(hasNonZeroHessian) << "La courbure géométrique devrait générer des Hessiennes physiques non nulles.";
}

// --- INSTANCIATION DES TESTS ---
// On crée les instances des éléments à tester
//1D elements
static P1ReferenceSegment  p1seg;
static P2ReferenceSegment  p2seg;

static P1ReferenceTriangle p1tri;
static P2ReferenceTriangle p2tri;

static Q1ReferenceQuad     q1quad;
static Q2ReferenceQuad     q2quad;

// 3D elements
static P1ReferenceTetrahedron p1tetra;
static P2ReferenceTetrahedron p2tetra;

static P1ReferencePrism       p1prism;
static P2ReferencePrism       p2prism;

static Q1ReferenceHexahedron  q1hex;
static Q2ReferenceHexahedron  q2hex;




INSTANTIATE_TEST_SUITE_P(
    AllElements,
    ReferenceElementTest,
    ::testing::Values(&p1seg, &p2seg, &p1tri, &p2tri, &q1quad, &q2quad, &p1tetra, &p2tetra, &p1prism, &p2prism, &q1hex, &q2hex)
);