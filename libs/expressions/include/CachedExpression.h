
#pragma once

#include "CfemExpression.h"
#include <memory>
#include <utility>

namespace cfem::sym {

    /**
     * @brief A wrapper expression that memoizes results in EvalContext.
     * Compatible with std::shared_ptr<const CfemExpression>.
     */
    class CachedExpression : public CfemExpression {
    public:
        // Use 'const CfemExpression' to match the library's return types
        explicit CachedExpression(std::shared_ptr<const CfemExpression> expr)
            : m_inner(std::move(expr)) {
            if (!m_inner) {
                throw std::runtime_error("CachedExpression: Inner expression cannot be null.");
            }
        }

        void prepare() const override {
            m_inner->prepare();
        }

        virtual CfemInt getDimension() const override {
            return m_inner->getDimension();
        }

        void evaluate(CfemInt cellId, const Vector3D &xi, 
                              std::span<CfemReal> out, EvalContext &ctx) const override {
            
            // Key is the address of this wrapper instance
            if (ctx.tryGetCache(this, out)) {
                return; 
            }

            m_inner->evaluate(cellId, xi, out, ctx);
            ctx.storeInCache(this, out);
        }

        bool isConstant() const override { return m_inner->isConstant(); }
        
        std::shared_ptr<const CfemExpression> grad() const override {
            return std::make_shared<const CachedExpression>(m_inner->grad());
        }

        std::shared_ptr<const CfemExpression> timeDerivative() const override {
            return std::make_shared<const CachedExpression>(m_inner->timeDerivative());
        }

        std::string toString() const override { return m_inner->toString(); }
        ExpressionType getType() const override { return m_inner->getType(); }

        

    private:
        std::shared_ptr<const CfemExpression> m_inner;
    };

    /**
     * @brief Helper function to wrap any expression shared_ptr in a cache layer.
     */
    inline std::shared_ptr<const CfemExpression> cache(std::shared_ptr<const CfemExpression> e) {
        return std::make_shared<const CachedExpression>(std::move(e));
    }
}