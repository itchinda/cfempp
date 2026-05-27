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
     * @brief 
     * 
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

    enum class ExpressionType {
        Constant,           /// A literal value or vector of values
        ScalarSymbolic,     /// A SymEngine/Symbolic variable (x, y, t)
        Variable,           /// A FEM field (u, v, p)
        Zero,               /// Specialized zero node for optimization
        Unary,              /// sin, cos, exp, log, sqrt, abs
        Binary,             /// +, -, *, /, power
        Scaled,             /// k * expression (the "optimization" node)
        VectorExpression,
        VectorSymbolic,     /// Container for 3 scalars
        TensorExpression,   /// Container for 9 scalars (or Voigt)
        TensorSymbolic,
        Gradient,           /// Grad(u)
        MagnitudeSq,        /// |u|^2
        Divergence,         /// div(u)
        Lambda,             /// Custom user-defined functions
        LinearEvaluator,
        FEField,
        Trial,
        Test
    };

    struct SpaceRequirements {
        std::vector<std::pair<const FESpace*, ShapeEvaluationFlags>> data;

        void addRequirement(const FESpace* space, ShapeEvaluationFlags flags) {
            for (auto& pair : data) {
                if (pair.first == space) {
                    pair.second |= flags; // Fusionne les drapeaux si l'espace existe
                    return;
                }
            }
            data.push_back({space, flags});
        }
    };

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class CfemExpression : public std::enable_shared_from_this<CfemExpression>
    {
    protected:
        
        mutable std::shared_ptr<const CfemExpression> m_gradient = nullptr;
        mutable std::shared_ptr<const CfemExpression> m_timeDeriv = nullptr;
        virtual std::shared_ptr<const CfemExpression> createGradient() const = 0;
        virtual std::shared_ptr<const CfemExpression> createTimeDerivative() const = 0;
        virtual void doCollectSpaceRequirements(SpaceRequirements& reqs,  ShapeEvaluationFlags inheritedFlags) const {}

    public:
        virtual ~CfemExpression() = default;

        virtual CfemInt  getNumComponents() const = 0;
        virtual ExpressionType getType() const = 0;

        virtual void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const = 0;
        virtual void evaluateBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext& ctx) const = 0;
        
        virtual void evaluateGradient(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const {
            auto g = this->grad(); 
            if (g && g.get() != this) {
                g->evaluate(cellId, xi, out, ctx);
                return;
            }
            CFEM_THROW("Gradient evaluation failed.");
        }

        virtual void evaluateGradientBatch(CfemInt cellId, CfemInt nQp, std::span<CfemReal> outValues, EvalContext &ctx) const {
            auto g = this->grad(); 
            if (g && g.get() != this) {
                g->evaluateBatch(cellId, nQp, outValues, ctx);
                return;
            }
            CFEM_THROW("Gradient batch evaluation failed.");
        }

        virtual void prepare() const {};

        virtual bool isPointwise() const { return false; }
        virtual bool isBinary() const    { return false; }
        virtual bool isZero() const      { return false; }
        virtual bool isContinuous() const { return true; }

        virtual bool isConstant() const  { return false; }
        virtual bool dependsOnFEField() const { return false; }
        virtual bool dependsOnShape() const { return true; }
        virtual bool dependsOnGeometricNormal() const { return false; }

        virtual void bind(const std::map<const FESpace*, int>& spaceToId) const {}

        virtual ShapeEvaluationFlags requiredShapeFlags() const {
            return ShapeEvaluationFlags::None;
        }

        virtual void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const {
            // Par défaut, pas de dépendance à un espace FE
        }

        void collectSpaceRequirements(SpaceRequirements& reqs, 
                                      ShapeEvaluationFlags inheritedFlags = ShapeEvaluationFlags::Value) const 
        {
            doCollectSpaceRequirements(reqs, inheritedFlags);
        }

        virtual std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const {return nullptr; } 

        std::shared_ptr<const CfemExpression> grad() const { 
            if ( !m_gradient ){
                m_gradient = createGradient();
            }
            if ( !m_gradient ){
                CFEM_THROW("Could not differentiate expression: " << toString());
            }
            return m_gradient;
         }

        std::shared_ptr<const CfemExpression> timeDerivative() const { 
            if ( !m_timeDeriv ){
                m_timeDeriv = createTimeDerivative();
            }
            if ( !m_timeDeriv ){
                CFEM_THROW("Could not differentiate expression: " << toString());
            }
            return m_timeDeriv;
        }

        virtual std::shared_ptr<const CfemExpression> div() const { return nullptr; } 

        virtual std::string toString() const { return "???"; }

        virtual CfemInt precedence() const {  return 10; }

    };

    inline std::ostream& operator<< (std::ostream& os, const std::shared_ptr<const CfemExpression>& expr){
        return os << expr->toString();
    }
    
    std::string opToSym(MathOp op);


}

#include <set>
#include <string>
#include <stdexcept>

namespace cfem::sym {
    class VariableRegistry {
        inline static std::set<std::string> reserved_names = {"x", "y", "z", "t"};
        inline static std::set<std::string> active_variables;

    public:
        static bool isReserved(const std::string& name) {
            return reserved_names.contains(name);
        }

        static void registerVariable(const std::string& name) {
            // If it's a reserved coordinate name, we allow it (it's "pre-validated")
            if (isReserved(name)) return;

            // Otherwise, check if it's already used by someone else to avoid collisions
            if (active_variables.contains(name)) {
                CFEM_THROW("Variable name '" << name << "' is already in use.");
            }
            active_variables.insert(name);
        }

        static void clear() {
            // If it's a reserved coordinate name, we allow it (it's "pre-validated")
            active_variables.clear();
        }
    };
}

