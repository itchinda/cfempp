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

#include <memory>
#include <regex>
#include <span>
#include <symengine/basic.h>
#include <symengine/expression.h>
#include <symengine/functions.h>
#include <symengine/lambda_double.h>
#include <symengine/symbol.h>
#include <vector>

#include "CfemExpression.h"
#include "cfemutils.h"
#include "EvalContext.h"
#include "SymbolicEvaluator.h"
#include "Vector3D.h"

namespace cfem::sym
{

    inline std::string finalizeString(std::string str) {
        // This regex looks for a decimal point followed by zeros at the end of a number
        // It turns 1.0 -> 1, 4.0 -> 4, and 3.14000 -> 3.14
        str = std::regex_replace(str, std::regex(R"(\.0+(?![0-9]))"), "");
        str = std::regex_replace(str, std::regex(R"((\.[0-9]*[1-9])0+(?![0-9]))"), "$1");
        return str;
    }
    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class ScalarSymbolicExpression : public CfemExpression
    {
    public:

        ScalarSymbolicExpression(const SymEngine::Expression                               &ex, 
                                 const std::vector<SymEngine::RCP<const SymEngine::Symbol>> &vars);
                                 
        CfemInt getNumComponents() const override { return 1; }
        bool isPointwise() const override { return true; }

        const SymEngine::Expression &get_expression() const { return m_expression; }
        const std::vector<SymEngine::RCP<const SymEngine::Symbol>> &get_vars() const { return m_vars; }

        void evaluate(CfemInt, const Vector3D &, std::span<CfemReal> out, EvalContext &ctx) const{
            const CfemReal input[4] = { ctx.singleWS.physPos.x, ctx.singleWS.physPos.y, ctx.singleWS.physPos.z, ctx.currentTime};
            out[0] = static_cast<CfemReal>(m_evaluator->evaluate(input));
        }


        void evaluateBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const override
        {
            const CfemReal t = ctx.currentTime;
            const auto* physPos_ptr = ctx.batchWS.physPos; 
            
            CFEM_ASSERT(physPos_ptr != nullptr && (*physPos_ptr).size() >= nQp, "Batch physical positions not properly sized in EvalContext");

            const auto& physPos = *physPos_ptr;

            for (CfemInt q = 0; q < nQp; ++q) {
                const CfemReal input[4] = { physPos[q].x, physPos[q].y, physPos[q].z, t };
                
                outValues[q] = m_evaluator->evaluate(input);
            }
        }

        std::shared_ptr<const ScalarSymbolicExpression> diff( const SymEngine::RCP<const SymEngine::Symbol> &var) const {
            return std::make_shared<const ScalarSymbolicExpression>(m_expression.diff(var), m_vars);
        }

        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx ) const override;

        std::shared_ptr<const CfemExpression> div() const override;

        std::string toString() const override { return finalizeString( m_expression.get_basic()->__str__() ); } // e.g., "x"
        ExpressionType getType() const override { return ExpressionType::ScalarSymbolic; }


    private:
        SymEngine::Expression                                m_expression;
        std::vector<SymEngine::RCP<const SymEngine::Symbol>> m_vars;
        std::unique_ptr<SymbolicEvaluator>                   m_evaluator;   //! Compiled evaluator owned by this expression instance
        
        std::shared_ptr<const CfemExpression> createGradient() const override;
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override;
    };




}