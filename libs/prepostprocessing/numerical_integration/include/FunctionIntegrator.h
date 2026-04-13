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

#ifndef CFEM_FUNCTIONINTEGRATOR_H
#define CFEM_FUNCTIONINTEGRATOR_H

#include "CfemExpression.h"
#include "Quadrature.h"
#include "cfem_topology.h"
#include "FESpace.h"
#include <memory>
#include <cmath>
#include <stdexcept>
#include <set>

namespace cfem
{
    class Mesh;
    class MeshEntity;
    
    class FunctionIntegrator
    {
    private:
        std::shared_ptr<const Mesh> m_mesh;
        IntegrationSchemeParams m_schemes_params;
        bool m_cacheInitialized = false;

        using EvaluationTarget = std::pair<std::shared_ptr<const sym::CfemExpression>, CfemReal*>;
        
        void integrate(const std::vector<EvaluationTarget> &targets, const MeshEntity &domain);
        
        void integrateBoundary(const std::vector<EvaluationTarget> &targets,
                               const MeshEntity                    &boundary,
                               const MeshEntity                    &volume);

        void dispatchIntegration(const std::vector<EvaluationTarget> &targets,
                                 const MeshEntity &entity,
                                 const MeshEntity &);

        inline bool needsMapping(ShapeEvaluationFlags flags)
        {
            return hasFlag(flags, ShapeEvaluationFlags::Gradient) ||
                   hasFlag(flags, ShapeEvaluationFlags::Hessian);
        }

    public:
        explicit FunctionIntegrator(std::shared_ptr<const Mesh> mesh,
                                    const IntegrationSchemeParams& params = {});
        
        void setIntegrationSchemeParams(const IntegrationSchemeParams& params){
            if (m_schemes_params != params) {
                m_schemes_params = params;
                m_cacheInitialized = false;
            }
        }

        void L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr, CfemReal* l2Out, CfemReal* h1SemiOut);
        void L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr, CfemReal* l2Out, CfemReal* h1SemiOut, const std::string& entity_name);
        void L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr, CfemReal* l2Out, CfemReal* h1SemiOut, CfemInt entity_id);
        void L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr, CfemReal* l2Out, CfemReal* h1SemiOut, const MeshEntity& entity);

        CfemReal computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr);
        CfemReal computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr, const std::string& entity_name);
        CfemReal computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr, CfemInt entity_id);
        CfemReal computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr, const MeshEntity& entity);
    };
}
#endif