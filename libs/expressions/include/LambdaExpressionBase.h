#pragma once

#include "CfemExpression.h"
#include "ZeroExpression.h"

namespace cfem::sym
{
    /**
     * @brief Classe de base non-template pour toutes les LambdaExpressions.
     * Elle permet d'identifier l'expression comme étant issue d'une fonction 
     * utilisateur personnalisée via dynamic_pointer_cast, sans avoir besoin 
     * de connaître le type exact (opaque) de la lambda fournie par l'utilisateur.
     */
    class LambdaExpressionBase : public CfemExpression
    {
    protected:
        CfemInt m_dimension;

        /**
         * @brief Retourne le gradient analytique fourni par l'utilisateur (ou nullptr).
         */
        std::shared_ptr<const CfemExpression> createGradient() const override {
            if (m_gradient) {
                return m_gradient;
            }
            CFEM_ERROR << "You requested the gradient, but this lambda expression " \
                            "was not provided with its gradient.";

            return std::make_shared<const ZeroExpression>(3*m_dimension);
        }

        /**
         * @brief Retourne la dérivée temporelle analytique fournie par l'utilisateur.
         * Si non fournie, retourne une expression nulle (Zero).
         */
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override {
            if (m_timeDeriv) {
                return m_timeDeriv;
            }
            CFEM_ERROR << "You requested the time derivative, but this lambda expression " \
                            "was not provided with its time derivative.";
            return std::make_shared<const ZeroExpression>(m_dimension);
        }


        /**
         * @brief Constructeur de la base.
         * @param dimension Dimension de l'espace d'arrivée (1 pour scalaire, 3 pour vecteur, etc.)
         * @param gradient Expression analytique du gradient (optionnelle)
         * @param timeDeriv Expression analytique de la dérivée temporelle (optionnelle)
         */
        LambdaExpressionBase(CfemInt dimension,
                             std::shared_ptr<const CfemExpression> gradient = nullptr,
                             std::shared_ptr<const CfemExpression> timeDeriv = nullptr);

    public:

        virtual ~LambdaExpressionBase() override = default;

        CfemInt getDimension() const override { return m_dimension; }
        ExpressionType getType() const override { return ExpressionType::Lambda; }

        bool isPointwise() const override {  return true; }

        bool dependsOnFEField() const override {
            // Uniquement si le gradient ou la dérivée temporelle fournis en dépendent
            bool dep = false;
            if (m_gradient) dep = dep || m_gradient->dependsOnFEField();
            if (m_timeDeriv) dep = dep || m_timeDeriv->dependsOnFEField();
            return dep;
        }

        bool dependsOnShape() const override { return false; }

        ShapeEvaluationFlags requiredShapeFlags() const override {
            return ShapeEvaluationFlags::None;
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            if (m_gradient) m_gradient->collectRequiredSpaces(spaces);
            if (m_timeDeriv) m_timeDeriv->collectRequiredSpaces(spaces);
        }

        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            if (m_gradient) m_gradient->bind(spaceToId);
            if (m_timeDeriv) m_timeDeriv->bind(spaceToId);
        }
        
        std::string toString() const override {
            std::string typeStr = (m_dimension == 1) ? "Scalar" : 
                                  (m_dimension == 3) ? "Vector" : "Tensor";
            return "LambdaExpression<" + typeStr + ">";
        }
    };

} // namespace cfem::sym