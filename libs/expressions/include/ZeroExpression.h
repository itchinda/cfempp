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
#include "EvalContext.h"
#include <span>
#include <memory>

namespace cfem::sym
{

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class ZeroExpression : public CfemExpression
    {
    private:
        std::shared_ptr<const CfemExpression> createGradient() const override {
            return std::make_shared<ZeroExpression>(m_dim * 3);
        }
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override{
            return std::make_shared<ZeroExpression>(m_dim);
        }

        CfemInt m_dim {0};

    public:
        ZeroExpression (CfemInt dim = 0): m_dim(dim) {
        }

        int getNumComponents() const override { return m_dim; }
        bool isPointwise() const override { return true; }
        bool isZero() const override { return true; };
        bool isConstant() const override {return true; }

        void evaluate(CfemInt, const Vector3D &, std::span<CfemReal> out, EvalContext &) const override{
            std::ranges::fill(out, 0.0);
        }

        std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const override {
            return std::make_shared<ZeroExpression>(this->getNumComponents());
        }

        std::shared_ptr<const CfemExpression> div() const override;

        ExpressionType getType() const override { return ExpressionType::Zero; }

        std::string toString() const override;
         
    };

    namespace zeros{
        inline const auto Zero1D = std::make_shared<ZeroExpression>(1);
        inline const auto Zero2D = std::make_shared<ZeroExpression>(2);
        inline const auto Zero3D = std::make_shared<ZeroExpression>(3);
        inline const auto ZeroTensor = std::make_shared<ZeroExpression>(9);

        inline const std::shared_ptr<const CfemExpression>& Zero = Zero1D;
    }

}