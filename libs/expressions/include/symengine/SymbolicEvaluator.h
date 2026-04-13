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

#include <symengine/symengine_config.h>
#include <symengine/lambda_double.h>
#include <vector>
#include <memory>

#include "cfemutils.h"

#ifdef HAVE_SYMENGINE_LLVM
#include <symengine/llvm_double.h>
#endif

namespace cfem::sym
{

    /**
     * @ingroup Expressions
     * @class SymbolicEvaluator
     * @brief A lightweight wrapper to evaluate SymEngine symbolic expressions numerically.
     *
     * This class provides a uniform interface to evaluate a SymEngine expression
     * at a given set of numerical values. Depending on the build configuration,
     * it uses either LLVM-based JIT evaluation (fast) or LambdaDouble evaluation (portable).
     *
     * The evaluator supports the following features:
     * - Evaluate a symbolic expression at a given array or vector of variables.
     * - Enforces runtime size checks for variable arrays.
     * - Non-copyable, movable only.
     *
     * @note CfemReal must be `double` when using this class.
     */
    class SymbolicEvaluator
    {
    public:
        /**
         * @brief Construct a symbolic evaluator for a given expression and variable set.
         *
         * @param expr The SymEngine expression to evaluate.
         * @param vars The vector of SymEngine variable symbols corresponding to `expr`.
         *
         * Depending on the build, the evaluator is either JIT-compiled using LLVM or
         * interpreted using LambdaDoubleVisitor.
         */
        SymbolicEvaluator(const SymEngine::RCP<const SymEngine::Basic> &expr,
                          const std::vector<SymEngine::RCP<const SymEngine::Basic>> &vars)
            : m_expr(expr), m_nvars(vars.size())
        {
#ifdef HAVE_SYMENGINE_LLVM
            m_evaluator = std::make_unique<SymEngine::LLVMDoubleVisitor>();
            m_evaluator->init(vars, {expr});
#else
            m_evaluator = std::make_unique<SymEngine::LambdaDoubleVisitor<CfemReal>>();
            m_evaluator->init(vars, *expr);
#endif
        }

        /**
         * @brief Evaluate the expression at a given array of values.
         *
         * @param values Pointer to an array of `CfemReal` with length equal to the number of variables.
         * @return The computed value of the expression.
         */
        CfemReal evaluate(const CfemReal *values) const noexcept
        {
            CfemReal out = 0.0;
            m_evaluator->call(&out, values);
            return out;
        }

        /**
         * @brief Evaluate the expression at a given vector of values.
         *
         * Performs a runtime check that the vector size matches the number of variables.
         *
         * @param values A vector of `CfemReal` with size equal to the number of variables.
         * @return The computed value of the expression.
         *
         * @throws std::runtime_error if the input size does not match.
         */
        CfemReal evaluate(const std::vector<CfemReal> &values) const
        {
            CFEM_ASSERT(values.size() == m_nvars, "Evaluation point dimension mismatch.");
            return evaluate(values.data());
        }

        // Non-copyable
        SymbolicEvaluator(const SymbolicEvaluator &) = delete;
        SymbolicEvaluator &operator=(const SymbolicEvaluator &) = delete;

        // Movable
        SymbolicEvaluator(SymbolicEvaluator &&) noexcept = default;
        SymbolicEvaluator &operator=(SymbolicEvaluator &&) noexcept = default;

    private:
        SymEngine::RCP<const SymEngine::Basic> m_expr; ///< Stored expression (useful for debugging)
        std::size_t m_nvars;                           ///< Number of variables in the expression

#ifdef HAVE_SYMENGINE_LLVM
        std::unique_ptr<SymEngine::LLVMDoubleVisitor> m_evaluator; ///< LLVM-based evaluator
#else
        std::unique_ptr<SymEngine::LambdaDoubleVisitor<CfemReal>> m_evaluator; ///< Lambda-based evaluator
#endif
    };

} // namespace cfem