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


#include "BinaryExpression.h"
#include "ConstantExpression.h"
#include "LambdaExpression.h"
#include "UnaryExpression.h"

namespace cfem::sym::utils
{

    // Utilitaire interne pour appliquer une fonction mathématique sur une constante
    template <typename Func>
    inline std::shared_ptr<const CfemExpression> foldUnary(
        std::shared_ptr<const CfemExpression> u, 
        MathOp op, 
        Func func) 
    {
        if (auto cU = std::dynamic_pointer_cast<const ConstantExpression>(u)) {
            const auto& v = cU->getValues();
            std::vector<CfemReal> res(v.size());
            for(size_t i = 0; i < v.size(); ++i) res[i] = func(v[i]);
            return std::make_shared<ConstantExpression>(std::move(res));
        }
        return std::make_shared<UnaryExpression>(u, op);
    }

    // Helper pour plier (fold) deux expressions constantes
    template <typename Func>
    inline std::shared_ptr<const CfemExpression> foldBinary(
        std::shared_ptr<const CfemExpression> lhs, 
        std::shared_ptr<const CfemExpression> rhs, 
        MathOp op, 
        CfemInt outDim,
        Func func) 
    {
        // On vérifie si lhs et rhs sont des constantes
        if (lhs->isConstant() && rhs->isConstant()) {
            auto cLhs = std::static_pointer_cast<const ConstantExpression>(lhs);
            auto cRhs = std::static_pointer_cast<const ConstantExpression>(rhs);
            
            const auto& vL = cLhs->getValues();
            const auto& vR = cRhs->getValues();
            std::vector<CfemReal> res(outDim);
            
            // On exécute la fonction lambda qui contient la logique mathématique
            func(vL, vR, res);
            
            return std::make_shared<ConstantExpression>(std::move(res));
        }
        
        // Si au moins un des deux n'est pas constant, on crée le noeud binaire
        return std::make_shared<BinaryExpression>(lhs, rhs, op);
    }
}

namespace cfem::sym {

    template <typename Func>
    requires std::is_invocable_r_v<CfemReal, Func, const Vector3D&>
    inline std::shared_ptr<const CfemExpression> scalar(
                                                    Func&& f,
                                                    std::shared_ptr<const CfemExpression> gradient)
    {
        return LambdaExpression<Func>::fromScalar(std::forward<Func>(f), std::move(gradient));
    }

    template <typename Func>
    requires std::is_invocable_r_v<CfemReal, Func, const Vector3D&, CfemReal> 
    inline std::shared_ptr<const CfemExpression> scalar(
                                                    Func&& f,
                                                    std::shared_ptr<const CfemExpression> gradient,
                                                    std::shared_ptr<const CfemExpression> timeDeriv)
    {
        return LambdaExpression<Func>::fromScalar(std::forward<Func>(f), std::move(gradient), std::move(timeDeriv));
    }

    template <typename Func>
    requires std::is_invocable_r_v<Vector3D, Func, const Vector3D&>
    inline std::shared_ptr<const CfemExpression> vector(
                                                    Func&& f,
                                                    std::shared_ptr<const CfemExpression> gradient)
    {
        return LambdaExpression<Func>::fromVector(std::forward<Func>(f), std::move(gradient));
    }

    template <typename Func>
    requires std::is_invocable_r_v<Vector3D, Func, const Vector3D&, CfemReal>
    inline std::shared_ptr<const CfemExpression> vector(
                                                    Func&& f,
                                                    std::shared_ptr<const CfemExpression> gradient,
                                                    std::shared_ptr<const CfemExpression> timeDeriv)
    {
        return LambdaExpression<Func>::fromVector(std::forward<Func>(f), std::move(gradient), std::move(timeDeriv));
    }

    template <typename Func>
    requires std::is_invocable_r_v<Tensor3D, Func, const Vector3D&>
    inline std::shared_ptr<const CfemExpression> tensor(
                                                    Func&& f,
                                                    std::shared_ptr<const CfemExpression> gradient)
    {
        return LambdaExpression<Func>::fromTensor(std::forward<Func>(f), std::move(gradient));
    }

    template <typename Func>
    requires std::is_invocable_r_v<Tensor3D, Func, const Vector3D&, CfemReal>
    inline std::shared_ptr<const CfemExpression> tensor(
                                                    Func&& f,
                                                    std::shared_ptr<const CfemExpression> gradient,
                                                    std::shared_ptr<const CfemExpression> timeDeriv)
    {
        return LambdaExpression<Func>::fromTensor(std::forward<Func>(f), std::move(gradient), std::move(timeDeriv));
    }

}