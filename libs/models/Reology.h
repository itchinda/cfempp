/**
 * @file CarreauYasudaModel.h
 * @brief High-performance implementation of the Carreau-Yasuda constitutive law.
 */

#pragma once
#include "CfemExpression.h"
#include "cfemutils.h"
#include "MathOperators.h"

namespace cfem::sym::models {

    /**
     * @class CarreauYasudaModel
     * @brief Represents the Carreau-Yasuda rheological model for non-Newtonian fluids.
     * * This class provides symbolic expressions for the dynamic viscosity @f$ \eta @f$ 
     * and its linearization (Gâteaux derivative). It is optimized for Finite Element assembly 
     * by decoupling the constitutive math from the kinematic strain rate calculations.
     * * @section Model Definition
     * The viscosity is defined as:
     * @f[
     * \eta(\dot{\gamma}) = \eta_\infty + (\eta_0 - \eta_\infty) \left[ 1 + (\lambda \dot{\gamma})^a \right]^{\frac{n-1}{a}}
     * @f]
     * where @f$ \dot{\gamma} @f$ is the equivalent shear rate.
     */
    class CarreauYasudaModel {
        CfemReal m_eta0;    ///< Zero-shear-rate viscosity @f$ \eta_0 @f$ [Pa·s]
        CfemReal m_etaInf;  ///< Infinite-shear-rate viscosity @f$ \eta_\infty @f$ [Pa·s]
        CfemReal m_lambda;  ///< Relaxation time constant @f$ \lambda @f$ [s]
        CfemReal m_n;       ///< Power-law index @f$ n @f$
        CfemReal m_a;       ///< Yasuda exponent @f$ a @f$

    public:
        /**
         * @brief Constructor for the Carreau-Yasuda parameters.
         * @param eta0     Zero-shear-rate viscosity.
         * @param etaInf   Infinite-shear-rate viscosity.
         * @param lambda   Time constant.
         * @param n        Flow behavior index (n < 1 for shear-thinning).
         * @param a        Yasuda exponent (controls transition width).
         */
        CarreauYasudaModel(CfemReal eta0, CfemReal etaInf, CfemReal lambda, CfemReal n, CfemReal a = 2.0)
            : m_eta0(eta0), m_etaInf(etaInf), m_lambda(lambda), m_n(n), m_a(a) {}

        /**
         * @brief Computes the symbolic viscosity @f$ \eta(\dot{\gamma}) @f$ from a pre-computed shear rate.
         * * @param g Symbolic expression of the shear rate @f$ \dot{\gamma} = \sqrt{2 \mathbf{D}:\mathbf{D}} @f$.
         * @return std::shared_ptr<const CfemExpression> Viscosity expression.
         */
        std::shared_ptr<const CfemExpression> modelFromShear(std::shared_ptr<const CfemExpression> g) const {
            auto term = pow(1.0 + pow(m_lambda * g, m_a), (m_n - 1.0) / m_a);
            return m_etaInf + (m_eta0 - m_etaInf) * term;
        }

        /**
         * @brief Computes the Gâteaux derivative @f$ D\eta(\mathbf{u})[\mathbf{v}] @f$ using pre-computed kinematics.
         * * This method implements the linearized viscosity for the Newton-Raphson tangent matrix.
         * The derivation follows the chain rule:
         * @f[
         * D\eta(\mathbf{u})[\mathbf{v}] = \frac{\partial \eta}{\partial \dot{\gamma}} \cdot D\dot{\gamma}(\mathbf{u})[\mathbf{v}]
         * @f]
         * * The viscosity sensitivity @f$ \frac{\partial \eta}{\partial \dot{\gamma}} @f$ is:
         * @f[
         * \frac{\partial \eta}{\partial \dot{\gamma}} = (\eta_0 - \eta_\infty) (n-1) \lambda^a \dot{\gamma}^{a-1} \left[ 1 + (\lambda \dot{\gamma})^a \right]^{\frac{n-1-a}{a}}
         * @f]
         * * The directional derivative of the shear rate @f$ D\dot{\gamma}(\mathbf{u})[\mathbf{v}] @f$ is:
         * @f[
         * D\dot{\gamma}(\mathbf{u})[\mathbf{v}] = \frac{2 \mathbf{D}(\mathbf{u}) : \mathbf{D}(\mathbf{v})}{\dot{\gamma}(\mathbf{u})}
         * @f]
         * @param g  Pre-computed equivalent shear rate @f$ \dot{\gamma}(\mathbf{u}) @f$.
         * @param Du Pre-computed solution strain rate @f$ \mathbf{D}(\mathbf{u}) @f$.
         * @param Dv Pre-computed test function strain rate @f$ \mathbf{D}(\mathbf{v}) @f$.
         * @return std::shared_ptr<const CfemExpression> The linearized viscosity expression.
         */
        std::shared_ptr<const CfemExpression> gateauDerivative(
            std::shared_ptr<const CfemExpression> g,
            std::shared_ptr<const CfemExpression> Du,
            std::shared_ptr<const CfemExpression> Dv) const 
        {
            CFEM_ASSERT(g && Du && Dv, "Inputs g, Du, and Dv must not be null.");

            // Constants for the power term: (n-1)/a - 1 = (n-1-a)/a
            auto power_term = (m_n - 1.0) / m_a - 1.0;
            
            // d_eta / d_gamma_dot
            auto d_eta_d_g = (m_eta0 - m_etaInf) * (m_n - 1.0) * std::pow(m_lambda, m_a) 
                           * pow(g, m_a - 1.0) 
                           * pow(1.0 + pow(m_lambda * g, m_a), power_term);

            // Gateaux derivative: (d_eta/dg) * (2 * Du : Dv / g)
            return d_eta_d_g * (2.0 * inner(Du, Dv) / g);
        }
    };

    /**
     * @class PowerLawModel
     * @brief Represents the Power Law rheological model for non-Newtonian fluids.
     * * The model defines viscosity @f$ \eta @f$ as a power function of the shear rate:
     * @f[
     * \eta(\dot{\gamma}) = K \dot{\gamma}^{n-1}
     * @f]
     * where @f$ K @f$ is the consistency index and @f$ n @f$ is the flow behavior index.
     */
    class PowerLawModel {
        CfemReal m_K; ///< Consistency index @f$ K @f$ [Pa·s^n]
        CfemReal m_n; ///< Power-law index @f$ n @f$

    public:
        /**
         * @brief Constructor for the Power Law parameters.
         * @param K Consistency index.
         * @param n Flow behavior index (n < 1: shear-thinning, n > 1: shear-thickening).
         */
        PowerLawModel(CfemReal K, CfemReal n) : m_K(K), m_n(n) {}

        /**
         * @brief Computes symbolic viscosity @f$ \eta(\dot{\gamma}) @f$.
         * @param g Symbolic expression of the shear rate @f$ \dot{\gamma} @f$.
         */
        std::shared_ptr<const CfemExpression> modelFromShear(std::shared_ptr<const CfemExpression> g) const {
            return m_K * pow(g, m_n - 1.0);
        }

        /**
         * @brief Computes the Gâteaux derivative @f$ D\eta(\mathbf{u})[\mathbf{v}] @f$.
         * * The viscosity sensitivity @f$ \frac{\partial \eta}{\partial \dot{\gamma}} @f$ is:
         * @f[
         * \frac{\partial \eta}{\partial \dot{\gamma}} = K (n-1) \dot{\gamma}^{n-2}
         * @f]
         * @param g  Pre-computed equivalent shear rate @f$ \dot{\gamma}(\mathbf{u}) @f$.
         * @param Du Pre-computed solution strain rate @f$ \mathbf{D}(\mathbf{u}) @f$.
         * @param Dv Pre-computed test function strain rate @f$ \mathbf{D}(\mathbf{v}) @f$.
         * @return std::shared_ptr<const CfemExpression> The linearized viscosity expression.
         */
        std::shared_ptr<const CfemExpression> gateauDerivative(
            std::shared_ptr<const CfemExpression> g,
            std::shared_ptr<const CfemExpression> Du,
            std::shared_ptr<const CfemExpression> Dv) const 
        {
            CFEM_ASSERT(g && Du && Dv, "Inputs g, Du, and Dv must not be null.");

            // d_eta / d_gamma_dot = K * (n-1) * g^(n-2)
            auto d_eta_d_g = m_K * (m_n - 1.0) * pow(g, m_n - 2.0);

            return d_eta_d_g * (2.0 * inner(Du, Dv) / g);
        }
    };


    /**
     * @class BinghamPapanastasiouModel
     * @brief Represents a regularized Bingham plastic model.
     * 
     * To handle the non-differentiability of the ideal Bingham model at @f$ \dot{\gamma} = 0 @f$, 
     * we use the Papanastasiou regularization:
     * @f[
     * \eta(\dot{\gamma}) = \mu_p + \frac{\tau_y}{\dot{\gamma}} \left( 1 - e^{-m\dot{\gamma}} \right)
     * @f]
     * where @f$ \tau_y @f$ is the yield stress, @f$ \mu_p @f$ is the plastic viscosity, 
     * and @f$ m @f$ is the regularization parameter.
     */
    class BinghamPapanastasiouModel {
        CfemReal m_tauY; ///< Yield stress @f$ \tau_y @f$ [Pa]
        CfemReal m_muP;  ///< Plastic viscosity @f$ \mu_p @f$ [Pa·s]
        CfemReal m_m;    ///< Regularization parameter @f$ m @f$ [s]

    public:
        /**
         * @brief Constructor for Bingham-Papanastasiou parameters.
         * @param tauY Yield stress.
         * @param muP  Plastic viscosity (viscosity after yielding).
         * @param m    Regularization parameter (typically 1e4 to 1e6).
         */
        BinghamPapanastasiouModel(CfemReal tauY, CfemReal muP, CfemReal m = 1000.0)
            : m_tauY(tauY), m_muP(muP), m_m(m) {}

        /**
         * @brief Computes symbolic viscosity @f$ \eta(\dot{\gamma}) @f$.
         * @param g Symbolic expression of the shear rate @f$ \dot{\gamma} @f$.
         */
        std::shared_ptr<const CfemExpression> modelFromShear(std::shared_ptr<const CfemExpression> g) const {
            auto exponential_term = 1.0 - exp(-m_m * g);
            return m_muP + (m_tauY / g) * exponential_term;
        }

        /**
         * @brief Computes the Gâteaux derivative @f$ D\eta(\mathbf{u})[\mathbf{v}] @f$.
         * * The viscosity sensitivity @f$ \frac{\partial \eta}{\partial \dot{\gamma}} @f$ is:
         * @f[
         * \frac{\partial \eta}{\partial \dot{\gamma}} = \frac{\tau_y}{\dot{\gamma}^2} \left[ (m\dot{\gamma} + 1) e^{-m\dot{\gamma}} - 1 \right]
         * @f]
         * * @param g  Pre-computed equivalent shear rate @f$ \dot{\gamma}(\mathbf{u}) @f$.
         * @param Du Pre-computed solution strain rate @f$ \mathbf{D}(\mathbf{u}) @f$.
         * @param Dv Pre-computed test function strain rate @f$ \mathbf{D}(\mathbf{v}) @f$.
         * @return std::shared_ptr<const CfemExpression> The linearized viscosity expression.
         */
        std::shared_ptr<const CfemExpression> gateauDerivative(
            std::shared_ptr<const CfemExpression> g,
            std::shared_ptr<const CfemExpression> Du,
            std::shared_ptr<const CfemExpression> Dv) const 
        {
            CFEM_ASSERT(g && Du && Dv, "Inputs g, Du, and Dv must not be null.");

            // d_eta / d_gamma_dot derivation:
            // d/dg [ (tauY/g) * (1 - exp(-m*g)) ]
            // = tauY * [ (m * exp(-m*g))/g - (1 - exp(-m*g))/g^2 ]
            // = (tauY / g^2) * [ (m*g + 1) * exp(-m*g) - 1 ]
            
            auto d_eta_d_g = (m_tauY / (g * g)) * ((m_m * g + 1.0) * exp(-m_m * g) - 1.0);

            return d_eta_d_g * (2.0 * inner(Du, Dv) / g);
        }
    };

}