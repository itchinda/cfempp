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

#include <cmath>
#include <memory>
#include <set>
#include <stdexcept>

#include "cfem_topology.h"
#include "CfemExpression.h"
#include "FESpace.h"
#include "Quadrature.h"

namespace cfem
{
    class Mesh;
    class MeshEntity;
    
    /**
     * @class ScalarFunctionalEvaluator
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
    class ScalarFunctionalEvaluator
    {
    private:
        std::shared_ptr<const Mesh> m_mesh;
        IntegrationSchemeParams m_schemes_params;
        bool m_cacheInitialized = false;

        using EvaluationTarget = std::pair<std::shared_ptr<const sym::CfemExpression>, CfemReal*>;
        
        /**
         * @brief Numerically integrates a set of scalar expressions over a mesh entity.
         *
         * @details
         * This function evaluates one or more scalar expressions at each quadrature point
         * of every cell belonging to the given mesh entity. It performs a highly optimized,
         * SIMD-friendly numerical quadrature and stores the final weighted sums into the 
         * corresponding output pointers.
         *
         * The integration is performed in three phases:
         * -# **Analysis & Binding** (O(1)): Inspects all targets to determine required
         * shape evaluation flags, geometric dependencies, and FE spaces. Binds each
         * expression to a compact space-index mapping.
         * -# **Reference Cache Build**: For each cell type present in the entity, precomputes
         * and caches reference shape functions, padded strides, and quadrature points. 
         * This cache is strictly read-only during the parallel phase.
         * -# **Parallel Integration**: Iterates over all cells using OpenMP. For each cell,
         * it computes the geometric mapping (Jacobian, determinant, physical derivatives),
         * evaluates all expressions in batch (SoA format), and computes the reduction.
         *
         * @par Affine Cell Optimization
         * For affine cells (linear geometry or explicitly flagged), the Jacobian and its
         * derived quantities (determinant, inverse, physical derivatives) are computed only
         * once per cell and hoisted out of the quadrature loop.
         *
         * @par Thread Safety & Parallelism
         * Each OpenMP thread maintains its own local evaluation context and accumulator array. 
         * Results are safely reduced into the final output pointers inside a critical section 
         * at the end of the entity loop.
         *
         * @par Arena Memory Management
         * Physical derivatives and expression intermediates are allocated via a zero-copy,
         * O(1) per-thread Arena allocator (@ref EvalContext). A strict bookmark/release 
         * mechanism guarantees that temporary AST buffers are reclaimed after each target, 
         * while SIMD-aligned physical derivatives are preserved for the cell's lifetime.
         *
         * @param[in] targets
         * A list of (expression, output pointer) pairs. Each expression must be scalar
         * (i.e., @ref CfemExpression::getNumComponents() == 1). If the output pointer
         * is non-null, its underlying value is **overwritten** (reset to zero before 
         * accumulation). Null output pointers are silently skipped.
         *
         * @param[in] entity
         * The mesh entity over which to integrate (defines the subset of cells to visit
         * via @ref MeshEntity::getCellIds()).
         *
         * @pre All expressions in @p targets must be compatible with the mesh's spatial
         * dimension and geometric order.
         * @pre Output pointers, if non-null, must remain valid for the duration of the call.
         *
         * @note If @p targets is empty or all expressions are null, the function returns
         * immediately with near-zero overhead.
         *
         * @warning This function is **not reentrant**: it must not be called concurrently
         * on the same @ref ScalarFunctionalEvaluator instance, as internal state
         * (e.g., scheme parameters) is shared across threads.
         *
         * @see EvalContext
         * @see MappingInfoBatchView
         * @see ReferenceElement::computeMappingBatch
         */
        void integrate(const std::vector<EvaluationTarget> &targets, const MeshEntity &entity);
        
        /**
         * @brief Numerically integrates a set of scalar expressions over a boundary mesh entity.
         *
         * @details
         * This function evaluates scalar expressions at the quadrature points of boundary facets.
         * Boundary integration in Finite Elements requires a highly optimized **Hybrid Geometry Architecture**:
         * the surface measure (dS) and physical normals are derived from the facet's geometry, while 
         * the inverse Jacobian required for physical derivatives (gradients) must be extracted from the 
         * parent volume cell.
         *
         * The integration is performed in three highly optimized phases:
         * -# **Analysis & Binding** (O(1)): Inspects all targets to determine required shape flags
         * (including geometric normal dependencies) and binds expressions to FE spaces.
         * -# **Hybrid Cache Build**: For each facet type and its parent cell type, precomputes 
         * reference shape functions on the facet, and simultaneously projects these quadrature 
         * points into the parent reference space to cache volume shape functions.
         * -# **Parallel Integration**: Iterates over all boundary facets using OpenMP. It dynamically
         * constructs a zero-copy **Hybrid View** (@ref MappingInfoBatchView) combining the facet's 
         * determinant/normals with the parent's inverse Jacobian, evaluates the expressions in 
         * SIMD-friendly batches, and reduces the results.
         *
         * @par Affine Cell Optimization
         * If the parent cell is affine, both the facet mapping and parent volume mapping (Jacobians)
         * are hoisted out of the quadrature loop and computed only once per facet.
         *
         * @par Arena Memory Management & Hybrid Views
         * Memory allocations for physical derivatives are handled by a zero-allocation, O(1) 
         * thread-local Arena (@ref EvalContext). The engine creates a Zero-Copy Hybrid View that 
         * perfectly isolates the math kernels from the geometric complexity of embedded manifolds.
         * Temporary AST evaluation buffers are efficiently reclaimed via a bookmark/release mechanism.
         *
         * @par Thread Safety & Parallelism
         * Each OpenMP thread utilizes its own local context and accumulator array. The final
         * boundary integrals are safely reduced into the output pointers within a critical section.
         *
         * @param[in] targets
         * A list of (expression, output pointer) pairs. Each expression must be scalar
         * (i.e., @ref CfemExpression::getNumComponents() == 1). If the output pointer
         * is non-null, its underlying value is **overwritten** (reset to zero before 
         * accumulation). Null output pointers are silently skipped.
         *
         * @param[in] boundary
         * The mesh entity representing the surface/boundary over which to integrate. 
         * Must contain valid facets (@ref MeshEntity::getFacetIds()).
         *
         * @param[in] volume
         * The mesh entity representing the ambient volume. This is required to identify 
         * the valid parent cell ("left" or "right" cell) connected to the boundary facet 
         * in order to compute spatial gradients.
         *
         * @pre All boundary facets must be strictly attached to at least one valid parent cell 
         * present in the @p volume entity. Unattached facets are silently ignored.
         * @pre Output pointers, if non-null, must remain valid for the duration of the call.
         *
         * @warning This function is **not reentrant**: it must not be called concurrently
         * on the same @ref ScalarFunctionalEvaluator instance.
         *
         * @see EvalContext
         * @see MappingInfoBatchView
         * @see ScalarFunctionalEvaluator::integrate
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

    public:
        /**
         * @brief Constructs a ScalarFunctionalEvaluator for a specific mesh.
         * @param mesh Shared pointer to the mesh to integrate over.
         * @param params Optional parameters to define quadrature orders per element type.
         */
        explicit ScalarFunctionalEvaluator(std::shared_ptr<const Mesh> mesh,
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
        
        /** @name Flux and Boundary Integral Computations */
        ///@{
        /**
         * @brief Computes the surface integral (flux) of an expression over a boundary.
         * @details This function evaluates surface integrals, typically used for Neumann 
         * boundary conditions, traction forces, or flow rates. It internally utilizes the 
         * highly optimized @ref integrateBoundary engine.
         * 
         * If the symbolic expression requires spatial gradients (e.g., calculating viscous 
         * stresses or heat flux \f$ \nabla u \cdot \mathbf{n} \f$), the parent @p volume_entity 
         * is strictly required to provide the correct inverse Jacobian for the chain rule.
         * 
         * @param expr The symbolic expression to integrate over the surface.
         * @param boundary_entity The boundary domain (facets) over which to integrate.
         * @param volume_entity The ambient volume domain containing the parent cells of the boundary facets.
         * @return The scalar result of the surface integral.
         */
        CfemReal computeFlux(const std::shared_ptr<const sym::CfemExpression> &expr, const MeshEntity& boundary_entity, const MeshEntity& volume_entity);
        
        /**
         * @brief Computes the surface integral over a boundary.
         * @note This overload implicitly uses the entire mesh ("Universe" entity) as the parent volume.
         * @param expr The symbolic expression to integrate.
         * @param boundary_entity The boundary domain.
         * @return The scalar result of the surface integral.
         */
        CfemReal computeFlux(const std::shared_ptr<const sym::CfemExpression> &expr, const MeshEntity& boundary_entity);
        
        /**
         * @brief Computes the surface integral using boundary and volume physical names.
         * @overload
         */
        CfemReal computeFlux(const std::shared_ptr<const sym::CfemExpression> &expr, const std::string& boundaryName, const std::string& volumeName);
        
        /**
         * @brief Computes the surface integral using a boundary physical name.
         * @note Implicitly uses the entire mesh ("Universe") as the parent volume.
         * @overload
         */
        CfemReal computeFlux(const std::shared_ptr<const sym::CfemExpression> &expr, const std::string& boundaryName);
        
        /**
         * @brief Computes the surface integral using boundary and volume physical IDs.
         * @overload
         */
        CfemReal computeFlux(const std::shared_ptr<const sym::CfemExpression> &expr, CfemInt boundaryId, CfemInt volumeId);
        
        /**
         * @brief Computes the surface integral using a boundary physical ID.
         * @note Implicitly uses the entire mesh ("Universe") as the parent volume.
         * @overload
         */
        CfemReal computeFlux(const std::shared_ptr<const sym::CfemExpression> &expr, CfemInt boundaryId);
        ///@}
    };
}