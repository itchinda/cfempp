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

#include "GradExpression.h"
#include "BinaryExpression.h"
#include "BinaryExpression.h"
#include "MagnitudeSqExpression.h"
#include "FEField.h"


namespace cfem::sym
{

    void GradExpression::evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const 
    {
        m_expr->evaluateGradient(cellId, xi, out, ctx);
        // const auto type = m_expr->getType();

        // // // 1. Chemin Analytique (Symbolique)
        // // // On ne le fait que si ce n'est pas un champ/trial/test pour éviter les boucles infinies
        // if (type != ExpressionType::FEField && type != ExpressionType::Trial && type != ExpressionType::Test) 
        // {
        //     if (auto analytical = m_expr->grad()) {
        //         if (analytical->getType() != ExpressionType::Gradient) {
        //             analytical->evaluate(cellId, xi, out, ctx);
        //             return;
        //         }
        //     }
        // }

        // if (type == ExpressionType::FEField) {
        //     static_cast<const FEField *>(m_expr.get())->evaluateGradient(cellId, xi, out, ctx);
        //     return;
        // }

        // 5. Cas Trial / Test (Fonctions de forme)
        // if (type == ExpressionType::Trial || type == ExpressionType::Test) {
        //     const bool isTrial = (type == ExpressionType::Trial);
        //     const auto* shape = isTrial ? ctx.trialShape : ctx.testShape;
        //     const auto* map   = isTrial ? ctx.trialMapping : ctx.testMapping;
        //     const int dof     = isTrial ? ctx.trialDof : ctx.testDof;

        //     const int numNodes = shape->values.size();
        //     const int comp = dof / numNodes;
        //     const int node = dof % numNodes;

        //     const auto& g = map->physicalShapeGradients[node];

        //     std::fill(out.begin(), out.end(), 0.0);
        //     out[comp*3 + 0] = g.x;
        //     out[comp*3 + 1] = g.y;
        //     out[comp*3 + 2] = g.z;
        //     return;
        // }

        // CFEM_ERROR << "GradExpression: Failed to differentiate expression tree.";
    }

}