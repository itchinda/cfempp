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

#include <functional>
#include <memory>
#include <span>
#include <vector>

#include "CfemExpression.h"
#include "cfemutils.h"
#include "EvalContext.h"
#include "Tensor3D.h"
#include "Vector3D.h"
#include "ZeroExpression.h"


namespace cfem::sym
{

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class ConstantExpression : public CfemExpression
    {
    public:
        ConstantExpression(std::vector<CfemReal> val)
            : m_value(std::move(val)) 
        {}

        ConstantExpression(CfemReal val)
            : m_value({val}) 
        {}
        
        ConstantExpression(const Vector3D& val)
            : m_value({val.x, val.y, val.z}) 
        {}
        
        ConstantExpression(const Tensor3D& val)
            : m_value({val(0,0), val(0,1), val(0,2),
                       val(1,0), val(1,1), val(1,2),
                       val(2,0), val(1,2), val(2,2)}) 
        {}

        ConstantExpression(CfemReal val, CfemInt dim)
            : m_value(dim, val)
        {}

        CfemInt getNumComponents() const override { return m_value.size(); }
        bool isPointwise() const override { return true; }
        bool isConstant() const override {return true; }
        bool isZero() const override {
            for (auto val : m_value) {
                if (! NumericPolicy::isZero(std::abs(val)) ) return false;
            }
            return true;
        }

        const std::vector<CfemReal>& getValues() const { return m_value; }

        void evaluate(CfemInt, const Vector3D &, std::span<CfemReal> out, EvalContext &) const override {
            std::copy(m_value.begin(), m_value.end(), out.begin());
        }

        inline void evaluateBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const override
        {
            CfemInt nComp = m_value.size();
            
            CFEM_ASSERT(outValues.size() >= static_cast<size_t>(nComp * nQp), 
                        "Batch evaluation output buffer too small for ConstantExpression");

            // =========================================================================
            // FORMAT SoA: outValues[composante * nQp + q]
            // Pour une constante, la même valeur est répétée sur tous les points q
            // =========================================================================
            for (CfemInt i = 0; i < nComp; ++i) {
                const CfemReal val = m_value[i];
                
                for (CfemInt q = 0; q < nQp; ++q) {
                    outValues[i * nQp + q] = val;
                }
            }
        }

        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const override {
            return std::make_shared<ZeroExpression>(this->getNumComponents());
        }

        std::shared_ptr<const CfemExpression> div() const override;

        CfemInt precedence() const override{ return 10; }

        ExpressionType getType() const override { return ExpressionType::Constant; }
        std::string toString() const override;
        
    private:
        std::vector<CfemReal> m_value;  //!
        std::shared_ptr<const CfemExpression> createGradient() const override {
            return std::make_shared<ZeroExpression>(getNumComponents() * 3);
        }
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override{
            return std::make_shared<ZeroExpression>(getNumComponents());
        }

    };

}