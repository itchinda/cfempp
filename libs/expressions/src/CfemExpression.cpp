#include "CfemExpression.h"

namespace cfem::sym
{
    std::string opToSym(MathOp op) {
        switch (op) {
            case MathOp::Add:   return "+";     break;
            case MathOp::Sub:   return "-";     break;
            case MathOp::Mul:   return "*";     break;
            case MathOp::Div:   return "/";     break;
            case MathOp::Dot:   return "dot";   break;
            case MathOp::Inner: return "Inner"; break;
            case MathOp::Cross: return "cross"; break;
            case MathOp::Pow:   return "pow";   break;
            case MathOp::Max:   return "max";   break; 
            case MathOp::Min:   return "min";   break;
            case MathOp::Atan2: return "atan2"; break;
            // Matrix Operations
            case MathOp::Transpose: return "transpose"; break;
            case MathOp::Trace:     return "trace"; break;
            case MathOp::Det:       return "det"; break;
            case MathOp::Inv:       return "inv"; break;
            // Basic Math
            case MathOp::Negate: return "-"; break;
            case MathOp::Abs:    return "abs"; break;
            case MathOp::Sqrt:   return "sqrt"; break;
            case MathOp::Cbrt:   return "cbrt"; break;
            // Logical / Step Functions
            case MathOp::Sign:   return "sign"; break;
            case MathOp::Heaviside: return "heaviside"; break;
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
                CFEM_ERROR << "Unsupported math operation.";
                break;
        }
        return "???";
    }
}