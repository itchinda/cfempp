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

    /** *************************************************************************
     * @brief 
     * 
     **************************************************************************** */
    class CfemExpression : public std::enable_shared_from_this<CfemExpression>
    {
    public:
        virtual ~CfemExpression() = default;

        virtual CfemInt getDimension() const = 0;

        virtual void evaluate(CfemInt cellId, const Vector3D &xi, std::span<CfemReal> out, EvalContext &ctx) const = 0;

        virtual void prepare() const {};

        virtual bool isPointwise() const { return false; }
        virtual bool isBinary() const { return false; }
        virtual bool isZero() const {return false; }
        virtual bool isConstant() const {return false; }

        virtual std::shared_ptr<const CfemExpression> grad() const { return nullptr; }
    };

}