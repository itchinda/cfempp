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

#include <type_traits>

#include "LambdaExpressionBase.h"

namespace cfem::sym
{
    /**
     * @brief Final template implementation of Lambda expressions.
     * @tparam F The specific type of the lambda/functor (deduced by the compiler).
     */
    template <typename F>
    class LambdaExpression : public LambdaExpressionBase
    {
    private:
        F m_f; // The callable is stored directly (potentially inlined)

    public:
        LambdaExpression(F&& f, 
                         CfemInt dimension,
                         std::shared_ptr<const CfemExpression> gradient = nullptr,
                         std::shared_ptr<const CfemExpression> timeDeriv = nullptr)
            : LambdaExpressionBase(dimension, std::move(gradient), std::move(timeDeriv)),
              m_f(std::forward<F>(f))
        {}

        /**
         * @brief Core evaluation logic. 
         * Uses compile-time dispatch (if constexpr) to handle different lambda signatures.
         */
        void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override
        {
            CFEM_ASSERT(out.size() >= m_dimension);
            
            // 1. Invoke the lambda based on its signature
            auto result = [&]() {
                // Signature: f(Vector3D, CfemReal) -> Spatiotemporal
                if constexpr (std::is_invocable_v<F, const Vector3D&, CfemReal>) {
                    return m_f(ctx.singleWS.physPos, ctx.currentTime);
                } 
                // Signature: f(Vector3D) -> Pure Spatial
                else if constexpr (std::is_invocable_v<F, const Vector3D&>) {
                    return m_f(ctx.singleWS.physPos);
                }
                else {
                    static_assert(sizeof(F) == 0, "Invalid Lambda signature. Use f(X) or f(X, t).");
                }
            }();

            // 2. Map the result to the output span
            using ResultType = decltype(result);
            if constexpr (std::is_arithmetic_v<ResultType>) {
                out[0] = static_cast<CfemReal>(result);
            } 
            else if constexpr (std::is_same_v<ResultType, Vector3D>) {
                out[0] = result.x; out[1] = result.y; out[2] = result.z;
            }
            else if constexpr (std::is_same_v<ResultType, Tensor3D>) {
                const CfemReal* data = result.data;
                for(int i = 0; i < 9; ++i) out[i] = data[i];
            }
        }

        // ====================================================================
        // --- Vectorized Batch (SoA) Evaluation
        // ====================================================================
        void evaluateBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const override
        {
            CFEM_ASSERT(outValues.size() >= static_cast<size_t>(m_dimension * nQp), 
                        "LambdaExpression: Batch buffer too small.");
            
            // Accès direct au tableau des positions physiques
            const auto* physPos_ptr = ctx.batchWS.physPos; 
            CFEM_ASSERT(physPos_ptr != nullptr && (*physPos_ptr).size() >= nQp, "LambdaExpression requires physical positions in batchWS.");

            const auto& physPos = *physPos_ptr;

            
            const CfemReal t = ctx.currentTime;

            // Boucle sur les points de quadrature
            for (CfemInt q = 0; q < nQp; ++q) 
            {
                const Vector3D& pos = physPos[q];

                // 1. Appel du functor à la compilation
                auto result = [&]() {
                    if constexpr (std::is_invocable_v<F, const Vector3D&, CfemReal>) {
                        return m_f(pos, t);
                    } 
                    else if constexpr (std::is_invocable_v<F, const Vector3D&>) {
                        return m_f(pos);
                    }
                    else {
                        static_assert(sizeof(F) == 0, "Invalid Lambda signature.");
                    }
                }();

                // 2. Écriture au format SoA: out[composante * nQp + q]
                using ResultType = decltype(result);
                if constexpr (std::is_arithmetic_v<ResultType>) {
                    outValues[q] = static_cast<CfemReal>(result);
                } 
                else if constexpr (std::is_same_v<ResultType, Vector3D>) {
                    outValues[0 * nQp + q] = result.x;
                    outValues[1 * nQp + q] = result.y;
                    outValues[2 * nQp + q] = result.z;
                }
                else if constexpr (std::is_same_v<ResultType, Tensor3D>) {
                    // Le compilateur va dérouler (unroll) cette petite boucle
                    for(int i = 0; i < 9; ++i) {
                        outValues[i * nQp + q] = result.data[i];
                    }
                }
            }
        }

        // --- Static Factories ---

        static auto fromScalar(F&& f, std::shared_ptr<const CfemExpression> g = nullptr, std::shared_ptr<const CfemExpression> t = nullptr) {
            return std::make_shared<const LambdaExpression<std::decay_t<F>>>(std::forward<F>(f), 1, std::move(g), std::move(t));
        }

        static auto fromVector(F&& f, std::shared_ptr<const CfemExpression> g = nullptr, std::shared_ptr<const CfemExpression> t = nullptr) {
            return std::make_shared<const LambdaExpression<std::decay_t<F>>>(std::forward<F>(f), 3, std::move(g), std::move(t));
        }

        static auto fromTensor(F&& f, std::shared_ptr<const CfemExpression> g = nullptr, std::shared_ptr<const CfemExpression> t = nullptr) {
            return std::make_shared<const LambdaExpression<std::decay_t<F>>>(std::forward<F>(f), 9, std::move(g), std::move(t));
        }
    };
}