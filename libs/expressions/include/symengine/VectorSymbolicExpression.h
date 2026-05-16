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


namespace cfem::sym
{

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class VectorSymbolicExpression : public CfemExpression
    {
    public:
        VectorSymbolicExpression(
            std::shared_ptr<const CfemExpression> dx,
            std::shared_ptr<const CfemExpression> dy,
            std::shared_ptr<const CfemExpression> dz
        );

        CfemInt getNumComponents() const override { return 3; }
        bool isPointwise() const override { return true; }
        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const override;
        std::shared_ptr<const CfemExpression> div() const override;

        inline void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const override
        {
            m_components[0]->evaluate(cellId, xi, out.subspan(0, 1), ctx);
            m_components[1]->evaluate(cellId, xi, out.subspan(1, 1), ctx);
            m_components[2]->evaluate(cellId, xi, out.subspan(2, 1), ctx);
        }

        inline std::shared_ptr<const ScalarSymbolicExpression> getComponent(CfemInt i) const
        {
            CFEM_ASSERT(i < 3 && i >= 0, "Index out of bounds.");
            return m_components[i];
        }

        ExpressionType getType() const override { return ExpressionType::VectorSymbolic; }

        std::string toString() const override {
            std::ostringstream os;
            os << "[" << m_components[0]->toString() << ", " << m_components[1]->toString() << ", " << m_components[2]->toString() << "]";
            return os.str();
        }
    private:
        std::shared_ptr<const ScalarSymbolicExpression> m_components[3];
        std::shared_ptr<const CfemExpression> createGradient() const override;
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override;
    };

}