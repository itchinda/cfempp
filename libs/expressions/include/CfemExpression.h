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
#include "EvalContext.h"
#include "Vector3D.h"
#include "Tensor3D.h"
#include "FESpace.h"
#include <set>
#include <span>
#include <memory>
#include <ostream>


namespace cfem::sym
{

    /** *************************************************************************
     * @enum MathOp
     * @brief Defines the mathematical operations supported by the symbolic engine.
     * 
     * Includes standard binary operations, vector/tensor operations, and a 
     * comprehensive suite of unary transcendental and algebraic functions.
     **************************************************************************** */
    enum class MathOp {
        // Binary
        Add, Sub, Mul, Div, Dot, Inner, Cross, Outer, Pow, Max, Min, Atan2,
        // Unary
        Negate, Transpose, Trace, Det, Inv,
        Abs, Sqrt, Cbrt, 
        Sign, Heaviside,
        Sin, Cos, Tan, Cot, Sec, Csc, 
        Asin, Acos, Atan, Acot, Asec, Acsc, 
        Sinh, Cosh, Tanh, Asinh, Acosh, Atanh, 
        Exp, Log, Gamma, Erf, Erfc, Zeta
    };

    /** *************************************************************************
     * @enum ExpressionType
     * @brief Identifies the specific node type within the Abstract Syntax Tree (AST).
     * 
     * Used for fast Runtime Type Information (RTTI), expression folding, 
     * and targeted optimizations without relying on dynamic casting.
     **************************************************************************** */
    enum class ExpressionType {
        Constant,           ///< A literal value or vector of values
        ScalarSymbolic,     ///< A SymEngine/Symbolic spatial or temporal variable (x, y, t)
        Variable,           ///< A generic unknown variable
        Zero,               ///< Specialized zero node for rapid algebraic simplification
        Unary,              ///< Mathematical functions (sin, cos, exp, log, sqrt, abs)
        Binary,             ///< Arithmetic operators (+, -, *, /, power)
        Scaled,             ///< Optimized node for scalar multiplication (k * expression)
        VectorExpression,   ///< A dynamically evaluated vector
        VectorSymbolic,     ///< Container for 3 scalar symbolic variables
        TensorExpression,   ///< Container for 9 scalars (or Voigt notation)
        TensorSymbolic,     ///< Symbolic tensor
        Gradient,           ///< Spatial gradient of an expression: Grad(u)
        MagnitudeSq,        ///< Squared magnitude optimization: |u|^2
        Divergence,         ///< Divergence of a vector field: div(u)
        Lambda,             ///< Custom user-defined C++ lambda functions
        LinearEvaluator,    ///< Optimized linear combinations
        FEField,            ///< A Finite Element Field bounded to a mesh
        Trial,              ///< Trial function for bilinear forms
        Test                ///< Test function for weak formulations
    };

    /** *************************************************************************
     * @struct SpaceRequirements
     * @brief Tracks the Finite Element spaces and shape flags required by an expression.
     * 
     * During the preparation phase, the AST is traversed to collect all `FESpace` 
     * dependencies. This allows the integrator to pre-allocate caches and only 
     * compute the strictly necessary shape functions or derivatives.
     **************************************************************************** */
    struct SpaceRequirements {
        std::vector<std::pair<const FESpace*, ShapeEvaluationFlags>> data;

        /**
         * @brief Registers a dependency on an FE space with specific shape requirements.
         * @param space Pointer to the required FE Space.
         * @param flags Bitmask of required evaluations (e.g., Value | FirstDerivatives).
         * @note If the space is already registered, the flags are merged (bitwise OR).
         */
        void addRequirement(const FESpace* space, ShapeEvaluationFlags flags) {
            for (auto& pair : data) {
                if (pair.first == space) {
                    pair.second |= flags; 
                    return;
                }
            }
            data.push_back({space, flags});
        }
    };

    /** *************************************************************************
     * @class CfemExpression
     * @brief Abstract base class for all nodes in the symbolic Abstract Syntax Tree (AST).
     * 
     * This class provides the unified interface for evaluating mathematical expressions
     * at quadrature points. It supports both single-point scalar evaluation and 
     * high-performance, SIMD-friendly batch evaluation. It also handles automatic 
     * symbolic differentiation via lazy evaluation.
     **************************************************************************** */
    class CfemExpression
    {
    protected:
        
        mutable std::shared_ptr<const CfemExpression> m_gradient = nullptr;
        mutable std::shared_ptr<const CfemExpression> m_timeDeriv = nullptr;
        
        /**
         * @brief Factory method to generate the symbolic spatial gradient of this node.
         * @return A new expression tree representing the gradient.
         */
        virtual std::shared_ptr<const CfemExpression> createGradient() const = 0;
        
        /**
         * @brief Factory method to generate the symbolic time derivative of this node.
         * @return A new expression tree representing the time derivative.
         */
        virtual std::shared_ptr<const CfemExpression> createTimeDerivative() const = 0;
        
        /**
         * @brief Internal hook for AST nodes to register their FE space dependencies.
         */
        virtual void doCollectSpaceRequirements(SpaceRequirements& reqs,  ShapeEvaluationFlags inheritedFlags) const {}

    public:
        virtual ~CfemExpression() = default;

        /**
         * @brief Returns the number of scalar components (1 for scalar, 3 for vector, 9 for tensor).
         */
        virtual CfemInt  getNumComponents() const = 0;
        
        /**
         * @brief Returns the AST node type for RTTI and optimizations.
         */
        virtual ExpressionType getType() const = 0;

        /**
         * @brief Evaluates the expression at a single reference point.
         * @param cellId The global ID of the cell being evaluated.
         * @param xi The reference coordinates of the evaluation point.
         * @param out Span pointing to the pre-allocated output buffer.
         * @param ctx Thread-local evaluation context (Arena allocator).
         */
        virtual void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const = 0;
        
        /**
         * @brief Evaluates the expression for a batch of quadrature points (SIMD-friendly).
         * @param cellId The global ID of the cell being evaluated.
         * @param nQp The number of quadrature points in the batch.
         * @param outValues Span pointing to the pre-allocated SoA output buffer.
         * @param ctx Thread-local batch evaluation context.
         */
        virtual void evaluateBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const = 0;
        
        /**
         * @brief Directly evaluates the spatial gradient at a single point.
         * @note Defaults to lazily building the gradient tree and evaluating it.
         */
        virtual void evaluateGradient(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const {
            auto g = this->grad(); 
            if (g && g.get() != this) {
                g->evaluate(cellId, xi, out, ctx);
                return;
            }
            CFEM_THROW("Gradient evaluation failed.");
        }

        /**
         * @brief Directly evaluates the spatial gradient for a batch of points.
         */
        virtual void evaluateGradientBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext &ctx) const {
            auto g = this->grad(); 
            if (g && g.get() != this) {
                g->evaluateBatch(cellId, nQp, outValues, ctx);
                return;
            }
            CFEM_THROW("Gradient batch evaluation failed.");
        }

        /**
         * @brief Performs pre-computation and tree optimization before the integration loop begins.
         */
        virtual void prepare() const {};

        // --- AST Optimization & Traversal Flags ---
        virtual bool isPointwise() const { return false; }
        virtual bool isBinary() const    { return false; }
        virtual bool isZero() const      { return false; }
        virtual bool isContinuous() const { return true; }
        virtual bool isConstant() const  { return false; }
        virtual bool dependsOnFEField() const { return false; }
        virtual bool dependsOnShape() const { return true; }
        virtual bool dependsOnGeometricNormal() const { return false; }

        /**
         * @brief Binds required FE Spaces to unique indices used by the EvalContext.
         * @param spaceToId Mapping from FESpace pointers to compact integer IDs.
         */
        virtual void bind(const std::map<const FESpace*, int>& spaceToId) const {}

        /**
         * @brief Returns the minimum shape evaluation flags required by this specific node.
         */
        virtual ShapeEvaluationFlags requiredShapeFlags() const {
            return ShapeEvaluationFlags::None;
        }

        /**
         * @brief Traverses the AST to collect all unique FESpace instances required for evaluation.
         */
        virtual void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const {
            // By default, no dependency on an FE space.
        }

        /**
         * @brief Recursively collects the specific shape flags required for each FE Space.
         */
        void collectSpaceRequirements(SpaceRequirements& reqs, 
                                      ShapeEvaluationFlags inheritedFlags = ShapeEvaluationFlags::Value) const 
        {
            doCollectSpaceRequirements(reqs, inheritedFlags);
        }

        /**
         * @brief Computes the partial derivative with respect to a specific variable index.
         */
        virtual std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const {return nullptr; } 

        /**
         * @brief Returns the lazily-evaluated spatial gradient of this expression.
         * @return A cached shared pointer to the gradient AST.
         */
        std::shared_ptr<const CfemExpression> grad() const { 
            if ( !m_gradient ){
                m_gradient = createGradient();
            }
            if ( !m_gradient ){
                CFEM_THROW("Could not differentiate expression: " << toString());
            }
            return m_gradient;
         }

        /**
         * @brief Returns the lazily-evaluated time derivative of this expression.
         * @return A cached shared pointer to the time derivative AST.
         */
        std::shared_ptr<const CfemExpression> timeDerivative() const { 
            if ( !m_timeDeriv ){
                m_timeDeriv = createTimeDerivative();
            }
            if ( !m_timeDeriv ){
                CFEM_THROW("Could not differentiate expression: " << toString());
            }
            return m_timeDeriv;
        }

        /**
         * @brief Returns the divergence of the expression (applicable for vector fields).
         */
        virtual std::shared_ptr<const CfemExpression> div() const { return nullptr; } 

        /**
         * @brief Generates a human-readable string representation of the expression.
         */
        virtual std::string toString() const { return "???"; }

        /**
         * @brief Returns the operator precedence (used for safe string formatting with parentheses).
         */
        virtual CfemInt precedence() const {  return 10; }

    };

    /**
     * @brief Stream operator to easily print the symbolic expression.
     */
    inline std::ostream& operator<< (std::ostream& os, const std::shared_ptr<const CfemExpression>& expr){
        return os << expr->toString();
    }
    
    /**
     * @brief Helper function to convert a MathOp enum to its string symbol (e.g., MathOp::Add -> "+").
     */
    std::string opToSym(MathOp op);


}

#include <set>
#include <string>
#include <stdexcept>

namespace cfem::sym {
    
    /**
     * @class VariableRegistry
     * @brief A global registry to manage symbolic variable names and prevent collisions.
     * 
     * Ensures that standard spatial/temporal coordinates ("x", "y", "z", "t") are 
     * protected, and prevents the user from accidentally defining overlapping variables 
     * which would break symbolic differentiation.
     */
    class VariableRegistry {
        inline static std::set<std::string> reserved_names = {"x", "y", "z", "t"};
        inline static std::set<std::string> active_variables;

    public:
        /**
         * @brief Checks if a variable name is a reserved system coordinate.
         */
        static bool isReserved(const std::string& name) {
            return reserved_names.contains(name);
        }

        /**
         * @brief Registers a new variable name in the system.
         * @param name The requested variable name.
         * @throws CFEM_ERROR If the name is already actively in use by another symbolic object.
         */
        static void registerVariable(const std::string& name) {
            // If it's a reserved coordinate name, we allow it (it's "pre-validated")
            if (isReserved(name)) return;

            // Otherwise, check if it's already used by someone else to avoid collisions
            if (active_variables.contains(name)) {
                CFEM_THROW("Variable name '" << name << "' is already in use.");
            }
            active_variables.insert(name);
        }

        /**
         * @brief Clears all user-defined active variables from the registry.
         */
        static void clear() {
            active_variables.clear();
        }
    };
}