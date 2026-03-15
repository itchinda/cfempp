#pragma once

#include "cfemutils.h"
#include "EvalContext.h"
#include "Vector3D.h"
#include "TensorO2.h"
#include <span>
#include <memory>



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
        Negate, Transpose, Abs, Sqrt, Cbrt, 
        Sign, Heaviside,
        Sin, Cos, Tan, Cot, Sec, Csc, 
        Asin, Acos, Atan, Acot, Asec, Acsc, 
        Sinh, Cosh, Tanh, Asinh, Acosh, Atanh, 
        Exp, Log, Gamma, Erf, Erfc, Zeta
    };

    enum class ExpressionType {
        Constant,           // A literal value or vector of values
        ScalarSymbolic,     // A SymEngine/Symbolic variable (x, y, t)
        Variable,           // A FEM field (u, v, p)
        Zero,               // Specialized zero node for optimization
        Unary,              // sin, cos, exp, log, sqrt, abs
        Binary,             // +, -, *, /, power
        Scaled,             // k * expression (the "optimization" node)
        VectorSymbolic,     // Container for 3 scalars
        TensorSymbolic,     // Container for 9 scalars (or Voigt)
        Gradient,           // Grad(u)
        MagnitudeSq,        // |u|^2
        Divergence,         // div(u)
        Lambda,             // Custom user-defined functions
        FEField
    };

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class CfemExpression : public std::enable_shared_from_this<CfemExpression>
    {
    public:
        virtual ~CfemExpression() = default;

        virtual CfemInt getDimension() const = 0;
        virtual ExpressionType getType() const = 0;

        virtual void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const = 0;

        virtual void prepare() const {};

        virtual bool isPointwise() const { return false; }
        virtual bool isBinary() const { return false; }
        virtual bool isZero() const {return false; }
        virtual bool isConstant() const {return false; }

        virtual std::shared_ptr<const CfemExpression> grad() const { return nullptr; }

        virtual std::string toString() const { return "???"; }
        virtual int precedence() const {  return 10; }

    };

    inline void hash_combine(size_t& seed, size_t v) {
        // Standard hash_combine from Boost/C++ standards proposals
        seed ^= v + 0x9e3779b9 + (seed << 6) + (seed >> 2);
    }
    

    inline std::string opToSym(MathOp op) {
        switch (op) {
            case MathOp::Add: return "+"; break;
            case MathOp::Sub: return "-"; break;
            case MathOp::Mul: return "*"; break;
            case MathOp::Div: return "/"; break;
            case MathOp::Dot: return "dot"; break;
            case MathOp::Inner: return "Inner"; break;
            case MathOp::Cross: return "cross"; break;
            case MathOp::Pow: return "pow"; break;
            case MathOp::Max: return "max"; break; 
            case MathOp::Min: return "min"; break;
            case MathOp::Atan2: return "atan2"; break;
            // Matrix Operations
            case MathOp::Transpose: return "transpose"; break;
            // Basic Math
            case MathOp::Negate: return "-"; break;
            case MathOp::Abs:    return "abs"; break;
            case MathOp::Sqrt:   return "sqrt"; break;
            case MathOp::Cbrt:   return "cbrt"; break;
            // Logical / Step Functions
            case MathOp::Sign:   return "sign"; break;
            case MathOp::Heaviside: "Heaviside"; break;
            // Trigonometry
            case MathOp::Sin: return "sin"; break;
            case MathOp::Cos: return "cos"; break;
            case MathOp::Tan: return "tan"; break;
            case MathOp::Cot: return "cot"; break; // Relies on IEEE 754 infinity for div by 0
            case MathOp::Sec: return "sec"; break;
            case MathOp::Csc: return "csc"; break;
            // Inverse Trigonometry
            case MathOp::Asin: return "asin"; break;
            case MathOp::Acos: return "acos"; break;
            case MathOp::Atan: return "atan"; break;
            case MathOp::Acot: return "acot"; break;
            case MathOp::Asec: return "asec"; break;
            case MathOp::Acsc: return "acsc"; break;
            // Hyperbolic
            case MathOp::Sinh:  return "sinh" ; break;
            case MathOp::Cosh:  return "cosh" ; break;
            case MathOp::Tanh:  return "tanh" ; break;
            case MathOp::Asinh: return "asinh" ; break;
            case MathOp::Acosh: return "acosh" ; break;
            case MathOp::Atanh: return "atanh" ; break;
            // Exponential & Special Functions
            case MathOp::Exp:   return "exp"; break;
            case MathOp::Log:   return "log"; break;
            case MathOp::Gamma: return "Gamma"; break;
            case MathOp::Erf:   return "Erf"; break;
            case MathOp::Erfc:  return "Erfc"; break;
            case MathOp::Zeta:  return "Zeta"; break;

            default:
                CFEM_ERROR("Unsupported math operation.");
                break;
        }
        return "???";
    }


}

