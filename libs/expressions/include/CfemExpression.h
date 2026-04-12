#pragma once

#include "cfemutils.h"
#include "EvalContext.h"
#include "Vector3D.h"
#include "TensorO2.h"
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
        Add, Sub, Mul, Div, Dot, Inner, Cross, Pow, Max, Min, Atan2,
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

    public:
        virtual ~CfemExpression() = default;

        virtual CfemInt getDimension() const = 0;
        virtual ExpressionType getType() const = 0;

        virtual void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const = 0;
        virtual void evaluateGradient(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const {
            auto g = this->grad(); 
            if (g && g.get() != this) {
                g->evaluate(cellId, xi, out, ctx);
                return;
            }
            CFEM_ERROR << "Gradient evaluation failed.";
        }

        virtual void prepare() const {};

        virtual bool isPointwise() const { return false; }
        virtual bool isBinary() const    { return false; }
        virtual bool isZero() const      { return false; }
        virtual bool isContinuous() const { return true; }

        virtual bool isConstant() const  { return false; }
        virtual bool dependsOnFEField() const { return false; }
        virtual bool dependsOnShape() const { return true; }

        virtual void bind(const std::map<const FESpace*, int>& spaceToId) const {}

        virtual ShapeEvaluationFlags requiredShapeFlags() const {
            return ShapeEvaluationFlags::None;
        }

        virtual void collectRequiredSpaces(std::set<std::shared_ptr<const FESpace>>& spaces) const {
            // Par défaut, pas de dépendance à un espace FE
        }

        virtual std::shared_ptr<const CfemExpression> diff(CfemInt var_idx) const {return nullptr; } 

        std::shared_ptr<const CfemExpression> grad() const { 
            if ( !m_gradient ){
                m_gradient = createGradient();
            }
            if ( !m_gradient ){
                CFEM_ERROR << "Could not differentiate expression: " << toString();
            }
            return m_gradient;
         }

        std::shared_ptr<const CfemExpression> timeDerivative() const { 
            if ( !m_timeDeriv ){
                m_timeDeriv = createTimeDerivative();
            }
            if ( !m_timeDeriv ){
                CFEM_ERROR << "Could not differentiate expression: " << toString();
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
                CFEM_ERROR << "Variable name '" << name << "' is already in use.";
            }
            active_variables.insert(name);
        }
    };
}

