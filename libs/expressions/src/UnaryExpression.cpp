#include "MathOperators.h"
#include "ZeroExpression.h"
#include "UnaryExpression.h"

namespace cfem::sym
{

    std::shared_ptr<const CfemExpression> UnaryExpression::grad() const {
        auto u = m_expr;
        auto du = m_expr->grad(); // Appel récursif (Chain Rule)

        // Optimisation cruciale : si la dérivée interne est nulle, tout est nul.
        // Cela empêche l'arbre de grandir inutilement (ex: grad(sin(5)) -> 0).
        if (du->isZero()) {
            return du; 
        }

        // --- Constantes Utiles ---
        auto c_one  = std::make_shared<ConstantExpression>(1.0);
        auto c_two  = std::make_shared<ConstantExpression>(2.0);
        auto c_half = std::make_shared<ConstantExpression>(0.5);

        switch (m_op) {
            // ==========================================
            // Algèbre de Base & Matrices
            // ==========================================
            case MathOp::Negate:    
                return -du;
            case MathOp::Transpose: 
                return transpose(du);
            case MathOp::Abs:       
                return sign(u) * du;
            case MathOp::Sqrt:      
                // d/dx sqrt(u) = 1 / (2 * sqrt(u)) * du
                return du / (c_two * sqrt(u));
            case MathOp::Cbrt:      
                // d/dx cbrt(u) = 1/3 * u^(-2/3) * du
                return (std::make_shared<ConstantExpression>(1.0/3.0) * du) / (cbrt(u) * cbrt(u));

            // ==========================================
            // Fonctions Logiques / Discontinues
            // ==========================================
            case MathOp::Sign:
            case MathOp::Heaviside:
                // La dérivée (Dirac) est nulle presque partout. 
                // En FEM, on renvoie 0 (sauf si tu utilises une version régularisée).
                return zeros::Zero3D; 

            // ==========================================
            // Trigonométrie Classique
            // ==========================================
            case MathOp::Sin: return cos(u) * du;
            case MathOp::Cos: return -sin(u) * du;
            case MathOp::Tan: 
                // d/dx tan(u) = sec^2(u) * du = (1 + tan^2(u)) * du
                return (c_one + tan(u) * tan(u)) * du;
            case MathOp::Cot: 
                // d/dx cot(u) = -csc^2(u) * du = -(1 + cot^2(u)) * du
                return -(c_one + cot(u) * cot(u)) * du;
            case MathOp::Sec: 
                return (sec(u) * tan(u)) * du;
            case MathOp::Csc: 
                return -(csc(u) * cot(u)) * du;

            // ==========================================
            // Trigonométrie Inverse
            // ==========================================
            case MathOp::Asin: 
                // d/dx asin(u) = 1 / sqrt(1 - u^2) * du
                return du / sqrt(c_one - u * u);
            case MathOp::Acos: 
                // d/dx acos(u) = -1 / sqrt(1 - u^2) * du
                return -du / sqrt(c_one - u * u);
            case MathOp::Atan: 
                // d/dx atan(u) = 1 / (1 + u^2) * du
                return du / (c_one + u * u);
            case MathOp::Acot: 
                // d/dx acot(u) = -1 / (1 + u^2) * du
                return -du / (c_one + u * u);
            case MathOp::Asec: 
                // d/dx asec(u) = 1 / (|u| * sqrt(u^2 - 1)) * du
                return du / (abs(u) * sqrt(u * u - c_one));
            case MathOp::Acsc: 
                // d/dx acsc(u) = -1 / (|u| * sqrt(u^2 - 1)) * du
                return -du / (abs(u) * sqrt(u * u - c_one));

            // ==========================================
            // Fonctions Hyperboliques
            // ==========================================
            case MathOp::Sinh: return cosh(u) * du;
            case MathOp::Cosh: return sinh(u) * du;
            case MathOp::Tanh: 
                // d/dx tanh(u) = sech^2(u) * du = (1 - tanh^2(u)) * du
                return (c_one - tanh(u) * tanh(u)) * du;

            // ==========================================
            // Fonctions Hyperboliques Inverses
            // ==========================================
            case MathOp::Asinh: 
                // d/dx asinh(u) = 1 / sqrt(u^2 + 1) * du
                return du / sqrt(u * u + c_one);
            case MathOp::Acosh: 
                // d/dx acosh(u) = 1 / sqrt(u^2 - 1) * du
                return du / sqrt(u * u - c_one);
            case MathOp::Atanh: 
                // d/dx atanh(u) = 1 / (1 - u^2) * du
                return du / (c_one - u * u);

            // ==========================================
            // Exponentielles et Logarithmes
            // ==========================================
            case MathOp::Exp: return exp(u) * du;
            case MathOp::Log: return du / u;

            // ==========================================
            // Fonctions Spéciales
            // ==========================================
            case MathOp::Erf: {
                // d/dx erf(u) = 2/sqrt(pi) * exp(-u^2) * du
                CfemReal factor = 2.0 / std::sqrt(M_PI);
                auto c_erf = std::make_shared<ConstantExpression>(factor);
                return c_erf * exp(-(u * u)) * du;
            }
            case MathOp::Erfc: {
                // d/dx erfc(u) = -2/sqrt(pi) * exp(-u^2) * du
                CfemReal factor = -2.0 / std::sqrt(M_PI);
                auto c_erfc = std::make_shared<ConstantExpression>(factor);
                return c_erfc * exp(-(u * u)) * du;
            }
            case MathOp::Gamma:
                // Mathématiquement: d/dx Gamma(u) = Gamma(u) * Digamma(u) * du
                // Le C++ standard ne fournit pas std::digamma. 
                CFEM_ERROR("Gradient de Gamma non supporté (nécessite l'implémentation de la fonction Digamma).");
                return nullptr;
            case MathOp::Zeta:
                // La dérivée de la fonction Zeta de Riemann (Zeta') n'a pas de forme analytique simple
                // exprimable avec les autres fonctions de base.
                CFEM_ERROR("Gradient de Zeta non supporté analytiquement.");
                return nullptr;

            default:
                CFEM_ERROR("Opérateur unitaire non reconnu lors du calcul du gradient.");
                return nullptr;
        }
    }
}