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

        virtual CfemInt getNumComponents() const override {
            return m_inner->getNumComponents();
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
        bool isContinuous() const override {
            return  m_inner->isContinuous();
        }
        bool dependsOnFEField() const override { return m_inner->dependsOnFEField(); }
        bool dependsOnShape() const override { return m_inner->dependsOnShape(); }
        ShapeEvaluationFlags requiredShapeFlags() const override {
            return m_inner->requiredShapeFlags();
        }

        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const override {
            if (m_inner) m_inner->collectRequiredSpaces(spaces);
        }

        void bind(const std::map<const FESpace*, int>& spaceToId) const override {
            if (m_inner) m_inner->bind(spaceToId);
        }

        std::string toString() const override { return m_inner->toString(); }
        ExpressionType getType() const override { return m_inner->getType(); }

        

    private:
        std::shared_ptr<const CfemExpression> m_inner;

        std::shared_ptr<const CfemExpression> createGradient() const override {
            return std::make_shared<const CachedExpression>(m_inner->grad());
        }

        std::shared_ptr<const CfemExpression> createTimeDerivative() const override {
            return std::make_shared<const CachedExpression>(m_inner->timeDerivative());
        }
    };

    /**
     * @brief Helper function to wrap any expression shared_ptr in a cache layer.
     */
    inline std::shared_ptr<const CfemExpression> cache(std::shared_ptr<const CfemExpression> e) {
        return std::make_shared<const CachedExpression>(std::move(e));
    }
}