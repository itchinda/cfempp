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
    
    /**
     * @class FunctionIntegrator
     * @brief Engine for high-performance numerical integration of FE expressions.
     *
     * This class handles the reduction of symbolic expressions into scalar values (integrals).
     * It supports multi-target integration, allowing the calculation of several values 
     * (like L2 and H1 norms) in a single mesh sweep to maximize cache efficiency.
     *
     * Features:
     * - Hybrid OpenMP parallelism for large-scale meshes.
     * - Multi-level caching of shape functions and geometric mappings.
     * - Integration over volumes (Cells) and boundaries (Facets/Physical Groups).
     */
    class FunctionIntegrator
    {
    private:
        std::shared_ptr<const Mesh> m_mesh;
        IntegrationSchemeParams m_schemes_params;
        bool m_cacheInitialized = false;

        using EvaluationTarget = std::pair<std::shared_ptr<const sym::CfemExpression>, CfemReal*>;
        
        /**
         * @brief Core volume integration loop (OpenMP accelerated).
         */
        void integrate(const std::vector<EvaluationTarget> &targets, const MeshEntity &domain);
        
        /**
         * @brief Core boundary integration loop (OpenMP accelerated).
         * @details Handles facet-to-parent mapping for field evaluation on surfaces.
         */
        void integrateBoundary(const std::vector<EvaluationTarget> &targets,
                               const MeshEntity                    &boundary,
                               const MeshEntity                    &volume);

        /**
         * @brief Internal dispatcher that selects the appropriate integration strategy.
         * @details Decides whether to use volume integration or boundary integration 
         * based on the topological dimension of the provided entity.
         */
        void dispatchIntegration(const std::vector<EvaluationTarget> &targets,
                                 const MeshEntity &entity,
                                 const MeshEntity &);

        inline bool needsMapping(ShapeEvaluationFlags flags)
        {
            return hasFlag(flags, ShapeEvaluationFlags::Gradient) ||
                   hasFlag(flags, ShapeEvaluationFlags::Hessian);
        }

    public:
        /**
         * @brief Constructs a FunctionIntegrator for a specific mesh.
         * @param mesh Shared pointer to the mesh to integrate over.
         * @param params Optional parameters to define quadrature orders per element type.
         */
        explicit FunctionIntegrator(std::shared_ptr<const Mesh> mesh,
                                    const IntegrationSchemeParams& params = {});
        
        /**
         * @brief Updates the integration parameters (quadrature orders).
         * @details Changing params automatically invalidates existing caches to ensure 
         * that subsequent integrations use the correct number of quadrature points.
         */
        void setIntegrationSchemeParams(const IntegrationSchemeParams& params){
            if (m_schemes_params != params) {
                m_schemes_params = params;
                m_cacheInitialized = false;
            }
        }

        /** @name L2 and H1 Norm Calculations */
        ///@{
        /**
         * @brief Computes L2 and H1-semi norms of an expression in one pass.
         * @param expr The symbolic expression to evaluate (e.g., a FEField).
         * @param l2Out Pointer to store the L2 norm result: @f[ \sqrt{\int_{E} |u|^2 d\mathbf{x}} @f].
         * @param h1SemiOut Pointer to store the H1-semi norm result: @f[ \sqrt{\int_{E} |\nabla u|^2 d\mathbf{x}} @f].
         * @param entity (Optional) The specific mesh entity/domain to integrate over.
         * @note If entity is not provided, integration is performed over the whole mesh volume.
         */
        void L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr, CfemReal* l2Out, CfemReal* h1SemiOut);
        void L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr, CfemReal* l2Out, CfemReal* h1SemiOut, const std::string& entity_name);
        void L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr, CfemReal* l2Out, CfemReal* h1SemiOut, CfemInt entity_id);
        void L2Norm(const std::shared_ptr<const sym::CfemExpression> &expr, CfemReal* l2Out, CfemReal* h1SemiOut, const MeshEntity& entity);
        ///@}

        /** @name General Integral Computations */
        ///@{
        /**
         * @brief Computes the integral of a symbolic expression over a domain.
         * @param expr The expression to integrate.
         * @param entity (Optional) The domain (volume or boundary) defined by name, ID, or object.
         * @return The scalar result of the integration.
         */
        CfemReal computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr);
        CfemReal computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr, const std::string& entity_name);
        CfemReal computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr, CfemInt entity_id);
        CfemReal computeIntegral(const std::shared_ptr<const sym::CfemExpression> &expr, const MeshEntity& entity);
        ///@}
    };
}
#endif