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

/**
 * @file Elasticity.h
 * @brief Symbolic implementation of Elasticity Models such as the isotropic Hooke's Law.
 */

#pragma once

#include "CfemExpression.h"
#include "MathOperators.h"

namespace cfem::sym::models {

    /**
     * @class LinearElasticityModel
     * @brief Implements linear isotropic elasticity.
     * @f[ \mathbf{\sigma} = 2\mu \mathbf{\varepsilon} + \lambda \mathrm{tr}(\mathbf{\varepsilon})\mathbf{I} @f]
     */
    class LinearElasticityModel {
        CfemReal m_E;  ///< Young's Modulus
        CfemReal m_nu; ///< Poisson's ratio
        
        // Lame parameters computed at construction
        CfemReal m_mu; 
        CfemReal m_lambda;

    public:
        LinearElasticityModel(CfemReal E, CfemReal nu) : m_E(E), m_nu(nu) {
            m_mu = E / (2.0 * (1.0 + nu));
            m_lambda = (E * nu) / ((1.0 + nu) * (1.0 - 2.0 * nu));
        }

        /**
         * @brief Returns the Stress Tensor expression given a Strain Tensor eps.
         * @param eps The symmetric strain tensor @f$ \mathbf{\varepsilon} @f$.
         */
        std::shared_ptr<const CfemExpression> stress(std::shared_ptr<const CfemExpression> eps) const {
            CFEM_ASSERT(eps && eps->getDimension() == 9, "`eps` must be a second order tensor." );
            auto I = identityO2TensorExpr(); // Assuming 3D
            return 2.0 * m_mu * eps + m_lambda * trace(eps) * I;
        }

        /**
         * @brief The Gateaux derivative of the stress.
         * Since Hooke's law is linear, the derivative in direction de is simply the model itself.
         * @f[ D\mathbf{\sigma}(\mathbf{\varepsilon})[\delta\mathbf{\varepsilon}] = 2\mu \delta\mathbf{\varepsilon} + \lambda \mathrm{tr}(\delta\mathbf{\varepsilon})\mathbf{I} @f]
         */
        std::shared_ptr<const CfemExpression> gateauDerivative( std::shared_ptr<const CfemExpression> de) const 
        {
            CFEM_ASSERT(de && de->getDimension() == 9, "`de` must be a second order tensor." );
            auto I = identityO2TensorExpr();
            return 2.0 * m_mu * de + m_lambda * trace(de) * I;
        }
    };


    /**
     * @class NeoHookeanModel
     * @brief Represents a Neo-Hookean hyperelastic material for large strains.
     * * The First Piola-Kirchhoff stress @f$ \mathbf{P} @f$ is derived from the deformation gradient @f$ \mathbf{F} @f$:
     * @f[ \mathbf{P} = \mu (\mathbf{F} - \mathbf{F}^{-T}) + \lambda \ln(J) \mathbf{F}^{-T} @f]
     */
    class NeoHookeanModel {
        CfemReal m_mu;
        CfemReal m_lambda;

        /**
         * @brief Computes the First Piola-Kirchhoff stress tensor @f$ \mathbf{P} @f$.
         * * The formula for a compressible Neo-Hookean material is:
         * @f[
         * \mathbf{P} = \mu (\mathbf{F} - \mathbf{F}^{-T}) + \lambda \ln(J) \mathbf{F}^{-T}
         * @f]
         * where @f$ J = \det(\mathbf{F}) @f$.
         * * @param F The deformation gradient @f$ \mathbf{F} = \mathbf{I} + \nabla \mathbf{u} @f$.
         * @return std::shared_ptr<const CfemExpression> PK1 stress tensor expression.
         */
        std::shared_ptr<const CfemExpression> firstPiolaStress(std::shared_ptr<const CfemExpression> F) const {
            CFEM_ASSERT(F && F->getDimension() == 9, "Deformation gradient F must be valid.");

            auto J = det(F);
            auto FinvT = transpose(inv(F));
            
            // PK1 = mu*(F - F^-T) + lambda*ln(J)*F^-T
            return m_mu * (F - FinvT) + m_lambda * log(J) * FinvT;
        }

    public:
        NeoHookeanModel(CfemReal mu, CfemReal lambda) : m_mu(mu), m_lambda(lambda) {}

        /**
         * @brief Computes the First Piola-Kirchhoff stress tensor.
         * @param F The deformation gradient @f$ \mathbf{F} = \mathbf{I} + \nabla \mathbf{u} @f$.
         */
        std::shared_ptr<const CfemExpression> firstPiolaStress(std::shared_ptr<const CfemExpression> F) const {
            auto J = det(F);
            auto FinvT = transpose(inv(F));
            return m_mu * (F - FinvT) + m_lambda * log(J) * FinvT;
        }

        /**
         * @brief Gâteaux derivative of the PK1 stress (The Tangent Stiffness).
         * This is significantly more complex as it involves the derivative of the inverse transpose.
         */
        std::shared_ptr<const CfemExpression> gateauDerivative(
            std::shared_ptr<const CfemExpression> F,
            std::shared_ptr<const CfemExpression> dF) const 
        {
            // Implementation requires symbolic Fourth-order tensors or flattened indices.
            // For now, we rely on the symbolic chain rule of F.
            auto PK1 = firstPiolaStress(F);
            return PK1->diff(F, dF); 
        }
    };
}