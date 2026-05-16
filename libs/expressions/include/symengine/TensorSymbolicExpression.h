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

#include "cfemutils.h"
#include "Vector3D.h"
#include "CfemExpression.h"
#include "ScalarSymbolicExpression.h"
#include "EvalContext.h"
#include <span>
#include <memory>

#include <symengine/expression.h>
#include <symengine/lambda_double.h>

#include "SymbolicEvaluator.h"


namespace cfem::sym
{

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class TensorSymbolicExpression : public CfemExpression
    {
    public:
        TensorSymbolicExpression(
            std::shared_ptr<const ScalarSymbolicExpression> dxx, std::shared_ptr<const ScalarSymbolicExpression> dxy, std::shared_ptr<const ScalarSymbolicExpression> dxz,
            std::shared_ptr<const ScalarSymbolicExpression> dyx, std::shared_ptr<const ScalarSymbolicExpression> dyy, std::shared_ptr<const ScalarSymbolicExpression> dyz,
            std::shared_ptr<const ScalarSymbolicExpression> dzx, std::shared_ptr<const ScalarSymbolicExpression> dzy, std::shared_ptr<const ScalarSymbolicExpression> dzz);

        CfemInt getNumComponents() const override { return 9; }
        bool isPointwise() const override { return true; }

        void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const
        {
            CFEM_ASSERT(out.size() >= 9, "TensorSymbolicExpression: 'out' span is too small!");
            for (CfemInt i = 0; i < 9; ++i) m_components[i]->evaluate(cellId, xi, out.subspan(i, 1), ctx);
        }

        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const override;

        // TODO Implement grad()

        std::shared_ptr<const CfemExpression> div() const override;

        std::shared_ptr<const ScalarSymbolicExpression> getComponent(CfemInt i) const
        {
            CFEM_ASSERT(i>=0 && i< 9, "TensorSymbolicExpression: Index out of bounds.");
            return m_components[i];
        }

        ExpressionType getType() const override { return ExpressionType::TensorSymbolic; }

        std::string toString() const override {
            std::ostringstream os;
            os << "[[" << m_components[0]->toString() << ", " << m_components[1]->toString() << ", " << m_components[2]->toString() << "], "
               <<  "[" << m_components[3]->toString() << ", " << m_components[4]->toString() << ", " << m_components[5]->toString() << "], "
               <<  "[" << m_components[6]->toString() << ", " << m_components[7]->toString() << ", " << m_components[8]->toString() << "]]";
            return os.str();
        }
    private:
        std::shared_ptr<const ScalarSymbolicExpression> m_components[9];
        std::shared_ptr<const CfemExpression> createGradient() const override { 
            CFEM_ERROR << "Gradient is not yet supported for tensor expressions. It would required a 3rd order tensor which is not implemented.";
            return nullptr; 
        }
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override;
    };

}
