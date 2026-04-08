#pragma once

#include "cfemutils.h"
#include "Vector3D.h"
#include "CfemExpression.h"
#include "EvalContext.h"
#include <span>
#include <memory>
#include <sstream>

namespace cfem::sym
{

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class VectorExpression : public CfemExpression
    {
    private:
        std::shared_ptr<const CfemExpression> m_components[3];

    public:
        VectorExpression (std::shared_ptr<const CfemExpression> expr_x,
                          std::shared_ptr<const CfemExpression> expr_y,
                          std::shared_ptr<const CfemExpression> expr_z)
        { 
            m_components[0] = std::move(expr_x);
            m_components[1] = std::move(expr_y);
            m_components[2] = std::move(expr_z);

            for (int i = 0; i < 3; ++i) {
                if (!m_components[i] || m_components[i]->getDimension() != 1) {
                    CFEM_ERROR << "VectorExpression component " << i << " must be a valid scalar expression. Got dim = " << m_components[i]->getDimension();
                }
            }

        }

        CfemInt getDimension() const override { return 3; }

        bool isPointwise() const override{
            return m_components[0]->isPointwise() && m_components[1]->isPointwise() && m_components[2]->isPointwise();
        }

        bool isZero() const override{
            return m_components[0]->isZero() && m_components[1]->isZero() && m_components[2]->isZero();
        }

        bool isConstant() const override{
            return m_components[0]->isConstant() && m_components[1]->isConstant() && m_components[2]->isConstant();
        };

        inline void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override
        {
            CFEM_ASSERT(out.size() >= 3, "Expected a buffer of size >=9.");

            for (CfemInt i = 0; i < 3; ++i){
                m_components[i]->evaluate(cellId, xi, out.subspan(i, 1), ctx);
            }
        }

        inline std::shared_ptr<const CfemExpression> getComponent(int i) const
        {
            CFEM_ASSERT(i < 3 && i >= 0, "Index out of bounds.");
            return m_components[i];
        }

        bool dependsOnFEField() const override { 
            return m_components[0]->dependsOnFEField() || 
                   m_components[1]->dependsOnFEField() ||
                   m_components[2]->dependsOnFEField(); 
        }

        bool dependsOnShape() const override { 
            return m_components[0]->dependsOnShape() || 
                   m_components[1]->dependsOnShape() ||
                   m_components[2]->dependsOnShape(); 
        }

        ShapeEvaluationFlags requiredShapeFlags() const override {
            return m_components[0]->requiredShapeFlags() | 
                   m_components[1]->requiredShapeFlags() | 
                   m_components[2]->requiredShapeFlags();
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            m_components[0]->collectRequiredSpaces(spaces); 
            m_components[1]->collectRequiredSpaces(spaces); 
            m_components[2]->collectRequiredSpaces(spaces);   
        }

        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            m_components[0]->bind(spaceToId);
            m_components[1]->bind(spaceToId);
            m_components[2]->bind(spaceToId);
        }

        ExpressionType getType() const override { return ExpressionType::VectorExpression; }

        std::string toString() const override {
            std::ostringstream os;
            os << "[" << m_components[0]->toString() << ", " << m_components[1]->toString() << ", " << m_components[2]->toString() << "]";
            return os.str();
        }

        std::shared_ptr<const CfemExpression> timeDerivative() const override{
           return std::make_shared<const VectorExpression>( m_components[0]->timeDerivative(), m_components[1]->timeDerivative(), m_components[2]->timeDerivative() );
        }

        std::shared_ptr<const CfemExpression> grad() const override;
        std::shared_ptr<const CfemExpression> div() const override;
        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const {
            return std::make_shared<const VectorExpression>(m_components[0]->diff(var_idx), m_components[1]->diff(var_idx), m_components[2]->diff(var_idx));
        }

         
    };

}