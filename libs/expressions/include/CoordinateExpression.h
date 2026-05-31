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
#include "ConstantExpression.h"
#include "ZeroExpression.h"

namespace cfem::sym
{

    class CoordinateExpression : public CfemExpression
    {
        CfemInt m_dimIndex; // 0 for x, 1 for y, 2 for z and 3 for t 
        std::string m_name; // e.g. x, y, z, t

        std::shared_ptr<const CfemExpression> createGradient() const override
        {            
            Vector3D g{0, 0, 0};
            if (m_dimIndex == 3) return std::make_shared<ConstantExpression>(g); // grad(t) = [0,0,0]

            g[m_dimIndex] = 1.0;
            return std::make_shared<ConstantExpression>(g);
        }

        std::shared_ptr<const CfemExpression> createTimeDerivative() const override {
            if (m_name == "t" || m_dimIndex == 3) {
                return std::make_shared<ConstantExpression>(1.0);
            }
            return std::make_shared<ZeroExpression>(1);
        }

    public:
        CoordinateExpression(const std::string& name, CfemInt index) : m_name(name), m_dimIndex(index) {
            CFEM_ASSERT(m_dimIndex >=0 && m_dimIndex<4, "Coordinate index should be within 0, 1, 2, and 4" );
            NamingPolicy::validate(name);
            VariableRegistry::registerVariable(name);
        }

        void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override
        {
            if (m_dimIndex == 3) 
            {
                out[0] = ctx.currentTime; return;
            }

            const Vector3D& globalX = ctx.singleWS.physPos;
            out[0] = globalX[m_dimIndex];
        }

        inline void evaluateBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const override
        {
            CFEM_ASSERT(outValues.size() >= static_cast<size_t>(nQp), 
                        "CoordinateExpression: Batch output buffer too small.");

            if (m_dimIndex == 3) 
            {
                const CfemReal t = ctx.currentTime;
                for (CfemInt q = 0; q < nQp; ++q) {
                    outValues[q] = t;
                }
                return;
            }

            // Récupération du tableau contigu des positions physiques (généré à l'Étape 3.1)
            const auto* physPos_ptr = ctx.batchWS.physPos; 
            CFEM_ASSERT(physPos_ptr != nullptr && (*physPos_ptr).size() >= nQp, "CoordinateExpression requires physical positions in batchWS.");

            const auto& physPos = *physPos_ptr;


            for (CfemInt q = 0; q < nQp; ++q) {
                outValues[q] = physPos[q][m_dimIndex];
            }
        }

        CfemInt getNumComponents() const override { return 1;}
        CfemInt getIndex() const { return m_dimIndex; }
        
        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const override {
            if (var_idx == m_dimIndex) {
                return std::make_shared<ConstantExpression>(1.0);
            }
            return std::make_shared<ZeroExpression>(1);
        }

        std::string toString() const override { return m_name; }
        ExpressionType getType() const override { return ExpressionType::ScalarSymbolic; }
    };

}