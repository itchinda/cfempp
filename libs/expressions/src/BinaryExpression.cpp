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

#include "BinaryExpression.h"
#include "ConstantExpression.h"
#include "UnaryExpression.h"
#include "VectorExpression.h"
#include "MathOperators.h"

namespace cfem::sym
{    
    BinaryExpression::BinaryExpression(std::shared_ptr<const CfemExpression> lhs,
                                   std::shared_ptr<const CfemExpression> rhs,
                                   MathOp op)
    : m_lhs(std::move(lhs)), m_rhs(std::move(rhs)), m_op(op)
    {
        CFEM_ASSERT(m_lhs && m_rhs, "Operands cannot be null.");
        m_dimL = m_lhs->getDimension();
        m_dimR = m_rhs->getDimension();

        // Dimension Logic & Validation
        switch (m_op)
        {
        case MathOp::Add:
        case MathOp::Sub:
            // Rigid Physics: Only allow same-dimension operations
            if (m_dimL != m_dimR)
            {
                CFEM_ERROR << "Addition/Subtraction requires identical dimensions (L:" << std::to_string(m_dimL) << " R:" << std::to_string(m_dimR) << ").";
            }
            m_dimOut = m_dimL;
            break;

        case MathOp::Mul:
            if (m_dimL == 1 || m_dimR == 1)
            {
                m_dimOut = std::max(m_dimL, m_dimR);
            }
            else if (m_dimL == 9 && m_dimR == 9)
            {
                m_dimOut = 9;
            }
            else if (m_dimL == 9 && m_dimR == 3)
            {
                m_dimOut = 3;
            }
            // Cas Matching (Hadamard / Element-wise) : 3 * 3 -> 3
            else if (m_dimL == m_dimR)
            {
                m_dimOut = m_dimL;
            }
            else
            {
                CFEM_ERROR << "Multiplication: dimensions incompatibles (L:" << m_dimL << " R:" << m_dimR << ")";
            }
            break;
        case MathOp::Div:
        case MathOp::Pow:
        case MathOp::Max:
        case MathOp::Min:
        case MathOp::Atan2:
            // Broadcasting: Scalar with Vector/Tensor, or element-wise matching
            if (m_dimL != m_dimR && m_dimL != 1 && m_dimR != 1)
            {
                CFEM_ERROR << "Operation requires matching dimensions or a scalar operand.";
            }
            m_dimOut = std::max(m_dimL, m_dimR);
            break;
            
        case MathOp::Dot:
            if (m_dimL == 3 && m_dimR == 3) {      m_dimOut = 1;  } 
            else if (m_dimL == 9 && m_dimR == 3) { m_dimOut = 3; }
            else if (m_dimL == 3 && m_dimR == 9) { m_dimOut = 3; }
            else if (m_dimL == 9 && m_dimR == 9) { m_dimOut = 9; }
            else if (m_dimL == 1 || m_dimR == 1) { m_dimOut = std::max(m_dimL, m_dimR); }
            else {
                CFEM_ERROR << "Dot product: dimensions incompatibles (L:" << m_dimL << " R:" << m_dimR << ")";
            }
            break;

        case MathOp::Cross:
            if (m_dimL != 3 || m_dimR != 3)
            {
                CFEM_ERROR << "Cross product is only defined for 3D vectors.";
            }
            m_dimOut = 3;
            break;

        case MathOp::Outer:
            CFEM_ASSERT(m_dimL == 3 && m_dimR == 3, "Outer product requires two vectors of dimension 3.");
            m_dimOut = 9;
            break;

        case MathOp::Inner:
            // Double contraction (Frobenius inner product). 
            // Dimensions must match to perform component-wise sum.
            if (m_dimL != m_dimR)
            {
                CFEM_ERROR << "Inner product requires matching dimensions.";
            }
            m_dimOut = 1; 
            break;

        default:
            CFEM_ERROR << "Unsupported operation for BinaryExpression.";
        }
    }

    std::shared_ptr<const CfemExpression> BinaryExpression::diff(CfemInt var_idx) const {
        auto dL = m_lhs->diff(var_idx);
        auto dR = m_rhs->diff(var_idx);

        if (!dL || !dR) {
            CFEM_ERROR << "One of the children couldn't be differentiated";
            return nullptr;
        }

        // Optimization: If both derivatives are zero, the whole derivative is zero.
        // Assuming you have an isZero() virtual method on CfemExpression.
        if (dL->isZero() && dR->isZero()) {
            return std::make_shared<ZeroExpression>(m_dimOut);
        }

        switch (m_op) {
            case MathOp::Add:
                return dL + dR;

            case MathOp::Sub:
                return std::make_shared<const BinaryExpression>(dL, dR, MathOp::Sub);

            case MathOp::Mul:
                // Product Rule: (u*v)' = u'*v + u*v'
                if (dL->isZero()) return m_lhs * dR;
                if (dR->isZero()) return dL * m_rhs;
                return (dL * m_rhs) + (m_lhs * dR);

            case MathOp::Div:
                // Quotient Rule: (u/v)' = (u'*v - u*v') / v^2
                if (dR->isZero()) return dL / m_rhs; // Simplified if denominator is constant
                return ((dL * m_rhs) - (m_lhs * dR)) / (m_rhs * m_rhs);

            case MathOp::Dot:
                // Dot Product Rule: (u . v)' = u' . v + u . v'
                if (dL->isZero()) return dot(m_lhs, dR);
                if (dR->isZero()) return dot(dL, m_rhs);
                return dot(dL, m_rhs) + dot(m_lhs, dR);

            case MathOp::Inner:
                // Double contraction (A : B)' = A' : B + A : B'
                if (dL->isZero()) return inner(m_lhs, dR);
                if (dR->isZero()) return inner(dL, m_rhs);
                return inner(dL, m_rhs) + inner(m_lhs, dR);

            case MathOp::Cross:
                // Cross Product Rule: (u x v)' = u' x v + u x v'
                if (dL->isZero()) return cross(m_lhs, dR);
                if (dR->isZero()) return cross(dL, m_rhs);
                return cross(dL, m_rhs) + cross(m_lhs, dR);

            case MathOp::Pow: {
                // Power Rule: d(u^v) = u^v * (v' * ln(u) + v * u' / u)
                if (dR->isZero()) {
                    // If the exponent is a constant (e.g., x^2), optimize to: v * u^(v-1) * u'
                    auto v_minus_1 = m_rhs - std::make_shared<ConstantExpression>(1.0);
                    return m_rhs * pow(m_lhs, v_minus_1) * dL;
                }
                auto term1 = dR * log(m_lhs);
                auto term2 = m_rhs * dL / m_lhs;
                return pow(m_lhs, m_rhs) * (term1 + term2);
            }

            case MathOp::Max: {
                // max(u, v)' = u' * H(u - v) + v' * H(v - u)
                // Assuming you have a helper for Heaviside, e.g., heaviside(expr)
                auto H_u_gt_v = heaviside(m_lhs - m_rhs);
                auto H_v_gt_u = heaviside(m_rhs - m_lhs);
                return (dL * H_u_gt_v) + (dR * H_v_gt_u);
            }

            case MathOp::Min: {
                // min(u, v)' = u' * H(v - u) + v' * H(u - v)
                auto H_v_gt_u = heaviside(m_rhs - m_lhs);
                auto H_u_gt_v = heaviside(m_lhs - m_rhs);
                return (dL * H_v_gt_u) + (dR * H_u_gt_v);
            }

            case MathOp::Atan2: {
                // atan2(y, x)' = (x * y' - y * x') / (x^2 + y^2)
                // Note: lhs = y, rhs = x
                auto num = (m_rhs * dL) - (m_lhs * dR);
                auto den = (m_lhs * m_lhs) + (m_rhs * m_rhs);
                return num / den;
            }

            default:
                return std::make_shared<ZeroExpression>(m_dimOut);
        }
    }

    std::shared_ptr<const CfemExpression> BinaryExpression::createGradient() const {
        auto L = m_lhs;
        auto R = m_rhs;
        auto gL = L->grad();
        auto gR = R->grad();

        // Protection : si une branche n'est pas dérivable, l'opération globale ne l'est pas.
        if (!gL || !gR) {
            CFEM_ERROR << "One of the operands in BinaryExpression is not differentiable.";
            return nullptr;
        }

        switch (m_op) {
            // Linéarité : (f ± g)' = f' ± g'
            case MathOp::Add:
            case MathOp::Sub:
                return std::make_shared<BinaryExpression>(gL, gR, m_op);

            // Règle du Produit (Leibniz) : (f*g)' = f'*g + f*g'
            // S'applique à Mul, Dot (produit scalaire), Inner (contraction), Cross (produit vectoriel)
            case MathOp::Mul:
            case MathOp::Dot:
            case MathOp::Inner:
            case MathOp::Cross: {
                auto term1 = std::make_shared<BinaryExpression>(gL, R, m_op);
                auto term2 = std::make_shared<BinaryExpression>(L, gR, m_op);
                return std::make_shared<BinaryExpression>(term1, term2, MathOp::Add);
            }

            // Règle du Quotient : (f/g)' = (f'*g - f*g') / g²
            case MathOp::Div: {
                auto t1 = std::make_shared<BinaryExpression>(gL, R, MathOp::Mul);
                auto t2 = std::make_shared<BinaryExpression>(L, gR, MathOp::Mul);
                auto num = std::make_shared<BinaryExpression>(t1, t2, MathOp::Sub);
                auto den = std::make_shared<BinaryExpression>(R, R, MathOp::Mul);
                return std::make_shared<BinaryExpression>(num, den, MathOp::Div);
            }

            // Puissance Générale : (f^g)' = f^g * (g' * ln(f) + g * f' / f)
            case MathOp::Pow: {
                // Cas n * f^(n-1) * f' (Terme de base)
                auto un = std::make_shared<ConstantExpression>(1.0);
                auto gm1 = std::make_shared<BinaryExpression>(R, un, MathOp::Sub);
                auto f_pow_gm1 = std::make_shared<BinaryExpression>(L, gm1, MathOp::Pow);
                auto coeff = std::make_shared<BinaryExpression>(R, f_pow_gm1, MathOp::Mul);
                auto term1 = std::make_shared<BinaryExpression>(coeff, gL, MathOp::Mul);

                // Cas ln(f) * f^g * g' (Si l'exposant R est aussi une variable)
                auto ln_f = std::make_shared<UnaryExpression>(L, MathOp::Log);
                auto f_pow_g = std::make_shared<BinaryExpression>(L, R, MathOp::Pow);
                auto ln_f_f_g = std::make_shared<BinaryExpression>(ln_f, f_pow_g, MathOp::Mul);
                auto term2 = std::make_shared<BinaryExpression>(ln_f_f_g, gR, MathOp::Mul);

                return std::make_shared<BinaryExpression>(term1, term2, MathOp::Add);
            }

            // Atan2 : grad(atan2(y, x)) = (x*gy - y*gx) / (x^2 + y^2)
            case MathOp::Atan2: {
                auto x_gy = std::make_shared<BinaryExpression>(R, gL, MathOp::Mul);
                auto y_gx = std::make_shared<BinaryExpression>(L, gR, MathOp::Mul);
                auto num = std::make_shared<BinaryExpression>(x_gy, y_gx, MathOp::Sub);
                
                auto x2 = std::make_shared<BinaryExpression>(R, R, MathOp::Mul);
                auto y2 = std::make_shared<BinaryExpression>(L, L, MathOp::Mul);
                auto den = std::make_shared<BinaryExpression>(x2, y2, MathOp::Add);
                
                return std::make_shared<BinaryExpression>(num, den, MathOp::Div);
            }

            // Opérateurs Logiques / Comparaison
            // La dérivée est nulle presque partout (on considère ces fonctions comme des constantes par morceaux)

            // Max / Min : Non-lisse
            case MathOp::Max:
            case MathOp::Min:
                // Mathématiquement : (f > g) ? f' : g'
                // Nécessite un nœud de type "ConditionalExpression" que nous n'avons pas encore détaillé.
                CFEM_ERROR << "Gradient of Max/Min is non-smooth and requires conditional node support.";
                return nullptr;

            default:
                CFEM_ERROR << "Binary operation not supported for differentiation.";
                return nullptr;
        }
    }

    std::shared_ptr<const CfemExpression> BinaryExpression::div() const {
        using namespace cfem::sym::utils;

        if (m_dimOut == 1) {
            CFEM_ERROR << "Impossible de calculer la divergence d'une expression binaire scalaire : " 
                    << opToSym(m_op);
            throw std::logic_error("La divergence d'un scalaire n'est pas définie.");
        }

        switch (m_op) {
            // Distributivité Linéaire
            case MathOp::Add: 
                return m_lhs->div() + m_rhs->div();
                
            case MathOp::Sub: 
                return m_lhs->div() - m_rhs->div();

            // Règle du produit pour la multiplication
            // div(a * V) = a * div(V) + grad(a) . V
            case MathOp::Mul: {
                if (m_dimL == 1 && m_dimR == 3) {
                    auto term1 = m_lhs * m_rhs->div();
                    auto term2 = std::make_shared<const BinaryExpression>(m_lhs->grad(), m_rhs, MathOp::Dot);
                    return term1 + term2;
                }
                if (m_dimR == 1 && m_dimL == 3) {
                    auto term1 = m_rhs * m_lhs->div();
                    auto term2 = std::make_shared<const BinaryExpression>(m_rhs->grad(), m_lhs, MathOp::Dot);
                    return term1 + term2;
                }
                // Si c'est une multiplication composante par composante (ex: vecteur * vecteur),
                // on laisse le fallback universel ci-dessous s'en occuper !
                break; 
            }

            // 3. Cas impossibles pour la divergence (toujours scalaires)
            case MathOp::Dot:
            case MathOp::Inner:
                // Ces opérations retournent toujours une dimension de 1.
                // On ne devrait jamais arriver ici grâce au `if (m_dimOut == 1)` au tout début.
                break;

            // 4. Fallback Universel via diff()
            // Pour : Div, Cross, Pow, Max, Min, Atan2
            // Ton moteur symbolique sait déjà dériver ces opérations via diff().
            // On s'appuie dessus pour éviter de dupliquer de grosses formules mathématiques !
            default:
                break;
        }

        // LE MOTEUR DE DIVERGENCE GÉNÉRIQUE
        // On demande les dérivées spatiales exactes de l'expression binaire complète.
        auto dx = this->diff(0); 
        auto dy = this->diff(1); 
        auto dz = this->diff(2); 

        if (m_dimOut == 3) {
            // Le résultat est un Vecteur. div(V) = dVx/dx + dVy/dy + dVz/dz
            return extractComponent(dx, 0) + 
                extractComponent(dy, 1) + 
                extractComponent(dz, 2);
        }

        if (m_dimOut == 9) {
            // Le résultat est un Tenseur. La divergence donne un vecteur (somme par ligne).
            auto row1 = extractComponent(dx, 0) + extractComponent(dy, 1) + extractComponent(dz, 2);
            auto row2 = extractComponent(dx, 3) + extractComponent(dy, 4) + extractComponent(dz, 5);
            auto row3 = extractComponent(dx, 6) + extractComponent(dy, 7) + extractComponent(dz, 8);

            return std::make_shared<const VectorExpression>(row1, row2, row3);
        }

        CFEM_ERROR << "Erreur de dimension lors du calcul de la divergence.";
        return nullptr;
    }

    std::shared_ptr<const CfemExpression> BinaryExpression::createTimeDerivative() const {
        auto L = m_lhs;
        auto R = m_rhs;
        auto dL = L->timeDerivative();
        auto dR = R->timeDerivative();

        // Protection: if a branch isn't differentiable over time, the global op isn't either.
        if (!dL || !dR) {
            CFEM_ERROR << "One of the operands in BinaryExpression is not differentiable wrt time.";
            return nullptr;
        }

        switch (m_op) {
            // Linearity: d(f ± g)/dt = df/dt ± dg/dt
            case MathOp::Add:
            case MathOp::Sub:
                return std::make_shared<BinaryExpression>(dL, dR, m_op);

            // Product Rule: d(f*g)/dt = (df/dt)*g + f*(dg/dt)
            // Applies to all bilinear multiplications
            case MathOp::Mul:
            case MathOp::Dot:
            case MathOp::Inner:
            case MathOp::Cross: {
                auto term1 = std::make_shared<BinaryExpression>(dL, R, m_op);
                auto term2 = std::make_shared<BinaryExpression>(L, dR, m_op);
                return std::make_shared<BinaryExpression>(term1, term2, MathOp::Add);
            }

            // Quotient Rule: d(f/g)/dt = ((df/dt)*g - f*(dg/dt)) / g^2
            case MathOp::Div: {
                auto t1 = std::make_shared<BinaryExpression>(dL, R, MathOp::Mul);
                auto t2 = std::make_shared<BinaryExpression>(L, dR, MathOp::Mul);
                auto num = std::make_shared<BinaryExpression>(t1, t2, MathOp::Sub);
                auto den = std::make_shared<BinaryExpression>(R, R, MathOp::Mul);
                return std::make_shared<BinaryExpression>(num, den, MathOp::Div);
            }

            // General Power Rule: d(f^g)/dt = f^g * ( (dg/dt)*ln(f) + g*(df/dt)/f )
            case MathOp::Pow: {
                auto f_pow_g = std::make_shared<BinaryExpression>(L, R, MathOp::Pow);
                
                // Term A: g * f^(g-1) * (df/dt)
                auto un = std::make_shared<ConstantExpression>(1.0, 1);
                auto gm1 = std::make_shared<BinaryExpression>(R, un, MathOp::Sub);
                auto f_pow_gm1 = std::make_shared<BinaryExpression>(L, gm1, MathOp::Pow);
                auto term1 = std::make_shared<BinaryExpression>(
                    std::make_shared<BinaryExpression>(R, f_pow_gm1, MathOp::Mul), 
                    dL, MathOp::Mul);

                // Term B: ln(f) * f^g * (dg/dt)
                auto ln_f = std::make_shared<UnaryExpression>(L, MathOp::Log);
                auto term2 = std::make_shared<BinaryExpression>(
                    std::make_shared<BinaryExpression>(ln_f, f_pow_g, MathOp::Mul), 
                    dR, MathOp::Mul);

                return std::make_shared<BinaryExpression>(term1, term2, MathOp::Add);
            }

            // Atan2: d(atan2(y, x))/dt = (x*(dy/dt) - y*(dx/dt)) / (x^2 + y^2)
            case MathOp::Atan2: {
                auto x_dy = std::make_shared<BinaryExpression>(R, dL, MathOp::Mul);
                auto y_dx = std::make_shared<BinaryExpression>(L, dR, MathOp::Mul);
                auto num = std::make_shared<BinaryExpression>(x_dy, y_dx, MathOp::Sub);
                
                auto x2 = std::make_shared<BinaryExpression>(R, R, MathOp::Mul);
                auto y2 = std::make_shared<BinaryExpression>(L, L, MathOp::Mul);
                auto den = std::make_shared<BinaryExpression>(x2, y2, MathOp::Add);
                
                return std::make_shared<BinaryExpression>(num, den, MathOp::Div);
            }

            // Non-Smooth Operations
            case MathOp::Max:
            case MathOp::Min:
                CFEM_ERROR << "Time derivative of Max/Min requires conditional node support.";
                return nullptr;

            default:
                CFEM_ERROR << "Binary operation not supported for time differentiation.";
                return nullptr;
        }
    }

}