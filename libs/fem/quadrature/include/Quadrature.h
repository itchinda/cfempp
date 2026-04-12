#ifndef CFEM_QUADRATURE_H
#define CFEM_QUADRATURE_H

#include <map>
#include <utility>
#include <string>
#include "cfemutils.h"
#include "Vector3D.h"
#include "DynamicVector.h"
#include "cfem_topology.h"

namespace cfem {

/**
 * @brief Represents a single integration point in the reference domain.
 */
struct QuadraturePoint {
    Vector3D xi;    ///< Reference coordinates (\f$ \xi, \eta, \zeta \f$)
    CfemReal weight;  ///< Integration weight associated with this point
};

/**
 * @brief A collection of quadrature points and weights forming a numerical integration rule.
 */
class QuadratureRule
{
public:

    const DynamicVector<QuadraturePoint>& getPoints() const { 
        CFEM_ASSERT(m_points.size()>0, "le vecteur des points de quadrature est vide.");
        return m_points;
     }


    CfemInt size() const
    {
        return static_cast<CfemInt>(m_points.size());
    }

private:
    DynamicVector<QuadraturePoint> m_points;

    friend class QuadratureFactory;
};

/**
 * @brief Static factory for generating standard Gaussian quadrature rules.
 * Rules are defined on reference elements (e.g., unit triangle with vertices 
 * (0,0), (1,0), (0,1)).
 */
class QuadratureFactory {
public:

    enum class GaussLegendre1DRule{
        GAUSS_1Pt = 1,
        GAUSS_2Pt,
        GAUSS_3Pt,
        GAUSS_4Pt,
        GAUSS_5Pt
    };

    /** @brief Supported integration rules for triangular elements. */
    enum class TriangleRule { 
        GAUSS_1Pt = 1, ///< 1-point rule (Centroid), exact for linear polynomials
        GAUSS_3Pt = 3, ///< 3-point rule, exact for quadratic polynomials
        GAUSS_4Pt = 4, ///< 4-point rule, exact for 3rd order polynomials
        GAUSS_6Pt = 6, ///< 6-point rule, exact for 4th order polynomials
        GAUSS_7Pt = 7, ///< 7-point rule, exact for 5th order polynomials
        GAUSS_12Pt = 12, ///< 12-point rule, exact for 6th order polynomials
        GAUSS_16Pt = 16, ///< 16-point rule, exact for 8th order polynomials
        GAUSS_19Pt = 19, ///< 19-point rule, exact for 9th order polynomials
        GAUSS_25Pt = 25 ///< 25-point rule, exact for 10th order polynomials
        //TODO: ADD MORE
    };

    enum class TetRule { 
        GAUSS_1Pt  = 1, ///< 1-point rule (Centroid), exact for linear polynomials
        GAUSS_4Pt  = 4, ///< 4-point rule, exact for second order polynomials
        GAUSS_5Pt  = 5,  ///< 5-point rule, exact for third order polynomials
        GAUSS_11Pt = 11,  ///< 11-point rule, exact for fourth order polynomials
        GAUSS_15Pt = 15,  ///< 15-point rule, exact for fith order polynomials
        GAUSS_29Pt = 29,  ///< 29-point rule, exact for sixth order polynomials
        GAUSS_33Pt = 33  ///< 33-point rule, exact for seventh order polynomials
        // TODO: ADD MORE
    };

    static QuadratureRule getGaussLengendre1DRule(GaussLegendre1DRule rule);

    /**
     * @brief Retrieves a specific triangle quadrature rule.
     * @param rule The requested TriangleRule enum.
     * @return A QuadratureRule containing the reference points and weights.
     */
    static QuadratureRule getTriangleRule(TriangleRule rule);

    /**
     * @brief Retrieves a specific triangle quadrature rule.
     * @param rule The requested TriangleRule enum.
     * @return A QuadratureRule containing the reference points and weights.
     */
    static QuadratureRule getTetRule(TetRule rule);

    

    static QuadratureRule getQuadRule(int nPoints1D);

    static QuadratureRule getHexRule(int nPoints1D);

    static QuadratureRule getPrismRule(TriangleRule triRule, GaussLegendre1DRule element1DRule);



private:

    /** @brief Midpoint rule: 1 point at (1/3, 1/3), weight = 0.5 (Area of ref triangle). */
    static QuadratureRule create1DRuleGaussLegendre1();

    /** @brief Hammer's rule: 3 points at edge midpoints, weights = 1/6. */
    static QuadratureRule create1DRuleGaussLegendre2();

    /** @brief 4-point rule: includes a central point with a negative weight for higher accuracy. */
    static QuadratureRule create1DRuleGaussLegendre3();

    static QuadratureRule create1DRuleGaussLegendre4();
    static QuadratureRule create1DRuleGaussLegendre5();

    /** @brief Midpoint rule: 1 point at (1/3, 1/3), weight = 0.5 (Area of ref triangle). */
    static QuadratureRule createTriangleRuleGauss1();

    /** @brief Hammer's rule: 3 points at edge midpoints, weights = 1/6. */
    static QuadratureRule createTriangleRuleGauss3();

    /** @brief 4-point rule: includes a central point with a negative weight for higher accuracy. */
    static QuadratureRule createTriangleRuleGauss4();
    static QuadratureRule createTriangleRuleGauss6();
    static QuadratureRule createTriangleRuleGauss7();
    static QuadratureRule createTriangleRuleGauss12();
    static QuadratureRule createTriangleRuleGauss16();
    static QuadratureRule createTriangleRuleGauss19();
    static QuadratureRule createTriangleRuleGauss25();

    /**
     * @brief 1-point rule (Centroid). 
     * **Precision**: Degree 1 (Linear).
     * **Coordinates**: Barycenter of the tetrahedron.
     * **Weight**: Reference volume \f$ V_{ref} = 1/6 \f$.
     */
    static QuadratureRule createTetRuleGauss1();

    /**
     * @brief 4-point Keast rule.
     * **Precision**: Degree 2 (Quadratic).
     * **Formula**: \f$ a = \frac{5 + 3\sqrt{5}}{20}, b = \frac{5 - \sqrt{5}}{20} \f$.
     * **Usage**: Ideal for standard assembly of linear elements.
     */
    static QuadratureRule createTetRuleGauss4_deg2();


    /**
     * @brief 5-point Hammer-Stroud rule.
     * **Precision**: Degree 3 (Cubic).
     * **Features**: Central point with negative weight to achieve higher order.
     * **Usage**: Recommended for P2 (quadratic) elements or complex material models.
     */
    static QuadratureRule createTetRuleGauss5_deg3();
    static QuadratureRule createTetRuleGauss11_deg4();
    static QuadratureRule createTetRuleGauss15_deg5();
    static QuadratureRule createTetRuleGauss29_deg6();
    static QuadratureRule createTetRuleGauss33_deg7();
};

/**
 * @brief Manages and caches quadrature rules to prevent redundant allocations during assembly.
 *        This class acts as the central registry for numerical integration across the mesh.
 */
class IntegrationScheme {
public:
    /**
     * @brief Retrieves a cached quadrature rule for a given cell type and point count.
     * @param type The geometric type of the element (e.g., TRI3, TET4).
     * @param pointCount The number of points requested (must match a factory enum).
     * @return const reference to the cached QuadratureRule.
     */
    const QuadratureRule& getRule(CellType type, int pointCount) {
        auto key = std::make_pair(type, pointCount);
        
        auto it = m_cache.find(key);
        if (it != m_cache.end()) {
            return it->second;
        }

        m_cache[key] = fetchFromFactory(type, pointCount);
        return m_cache[key];
    }

private:
    using Key = std::pair<CellType,int>;
    std::map<Key, QuadratureRule> m_cache;

    QuadratureRule fetchFromFactory(CellType type, int count);

};

/**
 * @brief Configuration for quadrature orders per element type.
 */
struct IntegrationSchemeParams {
    int NumIntgPt1DElem    = 3;   /// 5th order accuracy
    int NumIntgPtTriElem   = 6;   /// 4th order accuracy
    int NumIntgPtQuadElem  = 9;   /// 5th order accuracy
    int NumIntgPtTetElem   = 11;  /// 4th order accuracy
    int NumIntgPtPrismElem = 18;  /// 4th order accuracy
    int NumIntgPtHexElem   = 27;  /// 5th order accuracy

    /**
     * @brief Helper to get the order for a specific element type.
     */
    int getNumIntgPt(CellType type) const {
        switch (type) {
            case CellType::LINE2:  
            case CellType::LINE3: return NumIntgPt1DElem;
            case CellType::TRI3:  
            case CellType::TRI6: return NumIntgPtTriElem;
            case CellType::QUAD4:  
            case CellType::QUAD9: return NumIntgPtQuadElem;
            case CellType::TET4:  
            case CellType::TET10: return NumIntgPtTetElem;
            case CellType::PRISM6:  
            case CellType::PRISM18: return NumIntgPtPrismElem;
            case CellType::HEX8:  
            case CellType::HEX27: return NumIntgPtHexElem;
            default: return 1;
        }
    }

    bool operator==(const IntegrationSchemeParams& other ) const {
        return NumIntgPt1DElem    == other.NumIntgPt1DElem    &&
               NumIntgPtTriElem   == other.NumIntgPtTriElem   &&
               NumIntgPtQuadElem  == other.NumIntgPtQuadElem  &&
               NumIntgPtTetElem   == other.NumIntgPtTetElem   &&
               NumIntgPtPrismElem == other.NumIntgPtPrismElem &&
               NumIntgPtHexElem   == other.NumIntgPtHexElem;

    }

    bool operator!=(const IntegrationSchemeParams& other ) const {
        return !(*this == other);
    }

};

}

#include "Quadrature.tpp"

#endif