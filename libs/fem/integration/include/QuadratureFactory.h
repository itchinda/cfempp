#ifndef CFEM_QUADRATURE_RULE_
#define CFEM_QUADRATURE_RULE_

#include <vector>
#include "cfemutils.h"
#include "Vector3D.h"

struct QuadraturePoint {
    Vector3D xi;
    double weight;
};

class QuadratureRule {
public:
    std::vector<QuadraturePoint> points;
    
    // Utilitaire pour récupérer uniquement les poids ou points si besoin
    CfemInt size() const { return static_cast<CfemInt>(points.size()); }
};

#endif