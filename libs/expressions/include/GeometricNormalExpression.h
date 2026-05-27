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
#include "EvalContext.h"
#include "Vector3D.h"
#include "EvaluationFlags.h"

namespace cfem::sym
{
    /**
     * @class GeometricNormalExpression
     * @brief Symbolic expression representing the local unit normal vector to a geometric entity.
     * 
     * This class acts as a terminal leaf in the expression tree. It retrieves the physical 
     * normal vector computed by the geometric mapping (MappingInfo) at the current integration point.
     * 
     * Key characteristics:
     * - Output Dimension: 3 (Vector components x, y, z).
     * - Dependency: Requires the 'Normal' flag during shape evaluation.
     * - Applicability: Surface manifolds (2D in 3D) or domain boundaries (d-1).
     */
    class GeometricNormalExpression : public CfemExpression
    {
    public:
        /**
         * @brief Default constructor.
         */
        GeometricNormalExpression() : CfemExpression() {}

        /** 
         * @brief Returns the output dimension of the expression.
         * @return Always 3 for the normal vector components.
         */
        CfemInt getNumComponents() const override { 
            return 3; 
        }

        /** 
         * @brief Returns the expression category.
         * @return ExpressionType::Variable as the normal changes across the geometry.
         */
        ExpressionType getType() const override { 
            return ExpressionType::Variable; 
        }

        virtual bool dependsOnGeometricNormal() const override { return true; }

        /** 
         * @brief Evaluates the normal vector at a specific quadrature point.
         * 
         * Extracts the pre-computed normal vector from the EvalContext's geometric mapping.
         * The mapping must have been initialized with the ShapeEvaluationFlags::Normal flag.
         * 
         * @param cellId Local ID of the element (unused).
         * @param xi Reference coordinates (unused).
         * @param out Output span to be filled with (nx, ny, nz).
         * @param ctx Evaluation context containing geometric mapping data.
         */
        void evaluate([[maybe_unused]] CfemInt cellId, 
                      [[maybe_unused]] const Vector3D &xi, 
                      std::span<CfemReal> out, 
                      EvalContext &ctx) const override 
        {
            CFEM_ASSERT(ctx.geometricMapping, "Geometric mapping must be set in EvalContext to evaluate normal.");
            
            const Vector3D& n = ctx.geometricMapping->normal;
            out[0] = n.x;
            out[1] = n.y;
            out[2] = n.z;
        }

        inline void evaluateBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const override
        {
            CFEM_ASSERT(outValues.size() >= static_cast<size_t>(3 * nQp), 
                        "GeometricNormalExpression: Batch output buffer too small for 3D normal vectors.");

            const auto* geoBatch = ctx.geometricMappingBatch;
            CFEM_ASSERT(geoBatch != nullptr, "Geometric mapping batch must be set in EvalContext to evaluate normal.");

            // =========================================================================
            // FORMAT SoA: [n_x0..n_xn, n_y0..n_yn, n_z0..n_zn]
            // Extraction directe depuis le tableau pré-calculé geoBatch->normals
            // =========================================================================
            for (CfemInt q = 0; q < nQp; ++q) 
            {
                const Vector3D& n = geoBatch->normals[q];
                
                outValues[0 * nQp + q] = n.x;
                outValues[1 * nQp + q] = n.y;
                outValues[2 * nQp + q] = n.z;
            }
        }

    protected:
        /** 
         * @brief Computes the gradient of the normal vector.
         * 
         * The gradient of the normal corresponds to the Shape Operator (curvature tensor).
         * @note Currently not implemented. Returns nullptr.
         */
        std::shared_ptr<const CfemExpression> createGradient() const override {
            return nullptr; 
        }

        /** 
         * @brief Computes the time derivative of the normal vector.
         * @return nullptr as the geometry is assumed stationary.
         */
        std::shared_ptr<const CfemExpression> createTimeDerivative() const override {
            return nullptr; 
        }

    public:
        /** 
         * @brief Defines the geometric data required by this expression.
         * @return ShapeEvaluationFlags::Normal to ensure the mapping computes the normal vector.
         */
        ShapeEvaluationFlags requiredShapeFlags() const override {
            return ShapeEvaluationFlags::FirstDerivatives;
        }

        /** 
         * @brief Indicates if the expression depends on the mesh geometry.
         * @return True, as the normal is derived from the mapping.
         */
        bool dependsOnShape() const override { 
            return true; 
        }

        /** 
         * @brief Returns a string representation for debugging.
         * @return "n_geom"
         */
        std::string toString() const override { 
            return "Mesh::geometric_normal"; 
        }

        /** 
         * @brief Operator precedence for string formatting.
         * @return 10 (High precedence as a terminal variable).
         */
        CfemInt precedence() const override { 
            return 10; 
        }

        /**
         * @brief Does not depend on any Finite Element field. 
         */
        bool dependsOnFEField() const override { return false; }
        
        /** 
         * @brief No specific FE spaces are required for evaluation. 
         */
        void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>&) const override {}
        
        /** 
         * @brief No binding required as it is space-independent. 
         */
        void bind(const std::map<const FESpace*, int>&) const override {}

        /** 
         * @brief Evaluated point-by-point at quadrature points. 
         */
        bool isPointwise() const override { return true; }
        
        /** 
         * @brief Assumed continuous within the element.
         */
        bool isContinuous() const override { return false; }
    };

} // namespace cfem::sym