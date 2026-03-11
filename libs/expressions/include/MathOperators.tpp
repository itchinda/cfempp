#pragma once


#include "BinaryExpression.h"
#include "ConstantExpression.h"
#include "GradExpression.h"
#include "LambdaExpression.h"
#include "MagnitudeSqExpression.h"
#include "UnaryExpression.h"
#include "VariableExpression.h"
#include "symengine/ScalarSymbolicExpression.h"
#include "symengine/SymEngineSymbols.h"
#include "symengine/VectorSymbolicExpression.h"
#include "symengine/TensorSymbolicExpression.h"
#include "ZeroExpression.h"

namespace cfem::sym
{

    // Utilitaire interne pour appliquer une fonction mathématique sur une constante
    template <typename Func>
    inline std::shared_ptr<const CfemExpression> foldUnary(
        std::shared_ptr<const CfemExpression> u, 
        MathOp op, 
        Func func) 
    {
        if (auto cU = std::dynamic_pointer_cast<const ConstantExpression>(u)) {
            const auto& v = cU->getValues();
            std::vector<CfemReal> res(v.size());
            for(size_t i = 0; i < v.size(); ++i) res[i] = func(v[i]);
            return std::make_shared<ConstantExpression>(std::move(res));
        }
        return std::make_shared<UnaryExpression>(u, op);
    }

    // Helper pour plier (fold) deux expressions constantes
    template <typename Func>
    std::shared_ptr<const CfemExpression> foldBinary(
        std::shared_ptr<const CfemExpression> lhs, 
        std::shared_ptr<const CfemExpression> rhs, 
        MathOp op, 
        CfemInt outDim,
        Func func) 
    {
        // On vérifie si lhs et rhs sont des constantes
        if (lhs->isConstant() && rhs->isConstant()) {
            auto cLhs = std::static_pointer_cast<const ConstantExpression>(lhs);
            auto cRhs = std::static_pointer_cast<const ConstantExpression>(rhs);
            
            const auto& vL = cLhs->getValues();
            const auto& vR = cRhs->getValues();
            std::vector<CfemReal> res(outDim);
            
            // On exécute la fonction lambda qui contient la logique mathématique
            func(vL, vR, res);
            
            return std::make_shared<ConstantExpression>(std::move(res));
        }
        
        // Si au moins un des deux n'est pas constant, on crée le noeud binaire
        return std::make_shared<BinaryExpression>(lhs, rhs, op);
    }

    // Fonction utilitaire pour générer un tenseur/vecteur de "1.0"
    inline std::shared_ptr<const CfemExpression> makeOnesLike(std::shared_ptr<const CfemExpression> u) {
        return std::make_shared<ConstantExpression>(std::vector<CfemReal>(u->getDimension(), 1.0));
    }

        inline std::shared_ptr<const CfemExpression> operator-(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u;
        return foldUnary(u, MathOp::Negate, [](CfemReal x) { return -x; });
    }

    inline std::shared_ptr<const CfemExpression> operator-(std::shared_ptr<const CfemExpression> lhs, 
                                                           std::shared_ptr<const CfemExpression> rhs) 
    {
        if (lhs->getDimension() != rhs->getDimension()) {
            CFEM_ERROR("Dimension mismatch in operator-");
        }

        if (rhs->isZero()) return lhs; // x - 0 = x
        if (lhs->isZero()) return -rhs; // 0 - x = -x

        return foldBinary(lhs, rhs, MathOp::Sub, lhs->getDimension(), 
            [](const auto& vL, const auto& vR, auto& res) {
                for(size_t i = 0; i < vL.size(); ++i) res[i] = vL[i] - vR[i];
            });
    }

    inline std::shared_ptr<const CfemExpression> operator+(std::shared_ptr<const CfemExpression> lhs, 
                                                           std::shared_ptr<const CfemExpression> rhs) 
    {
        if (lhs->getDimension() != rhs->getDimension()) {
            CFEM_ERROR("Dimension mismatch in operator+");
        }

        if (rhs->isZero()) return lhs;
        if (lhs->isZero()) return rhs;

        return foldBinary(lhs, rhs, MathOp::Add, lhs->getDimension(), 
            [](const auto& vL, const auto& vR, auto& res) {
                for(size_t i = 0; i < vL.size(); ++i) res[i] = vL[i] + vR[i];
            });
    }

    inline std::shared_ptr<const CfemExpression> operator*(std::shared_ptr<const CfemExpression> lhs, 
                                                        std::shared_ptr<const CfemExpression> rhs) 
    {
        CfemInt dimL = lhs->getDimension();
        CfemInt dimR = rhs->getDimension();

        if (dimL != dimR && dimL != 1 && dimR != 1) {
            CFEM_ERROR("Invalid operator*");
        }

        CfemInt outDim = std::max(dimL, dimR);

        if (lhs->isZero() || rhs->isZero()) {
            if (outDim == 1) return zeros::Zero1D;
            if (outDim == 3) return zeros::Zero3D;
            if (outDim == 9) return zeros::ZeroTensor;
            return std::make_shared<ZeroExpression>(outDim);
        }

        if (lhs->isConstant() && dimL == 1) {
            auto cLhs = std::static_pointer_cast<const ConstantExpression>(lhs);
            if (cLhs->getValues()[0] == 1.0) return rhs;
        }
        if (rhs->isConstant() && dimR == 1) {
            auto cRhs = std::static_pointer_cast<const ConstantExpression>(rhs);
            if (cRhs->getValues()[0] == 1.0) return lhs;
        }

        return foldBinary(lhs, rhs, MathOp::Mul, outDim, 
            [dimL, dimR, outDim](const auto& vL, const auto& vR, auto& res) {
                if (dimL == dimR) { 
                    for(size_t i=0; i < outDim; ++i) res[i] = vL[i] * vR[i];
                } else if (dimL == 1) { 
                    for(size_t i=0; i < outDim; ++i) res[i] = vL[0] * vR[i];
                } else { 
                    for(size_t i=0; i < outDim; ++i) res[i] = vL[i] * vR[0];
                }
            });
    }

    
    inline std::shared_ptr<const CfemExpression> operator/(std::shared_ptr<const CfemExpression> lhs, 
                                                           std::shared_ptr<const CfemExpression> rhs) 
    {
        CfemInt dimL = lhs->getDimension();
        CfemInt dimR = rhs->getDimension();

        if (dimL != dimR && dimR != 1) {
            CFEM_ERROR("Invalid operator/: Divisor must be a scalar or match numerator dimension. Got " 
                    + std::to_string(dimL) + " / " + std::to_string(dimR));
        }

        if (rhs->isZero()) {
            CFEM_ERROR("Division by Zero detected in symbolic expression setup!");
        }

        if (lhs->isZero()) {
            return lhs;
        }

        if (rhs->isConstant() && dimR == 1) {
            auto cRhs = std::static_pointer_cast<const ConstantExpression>(rhs);
            if (cRhs->getValues()[0] == 1.0) {
                return lhs;
            }
        }

        return foldBinary(lhs, rhs, MathOp::Div, dimL, 
            [dimL, dimR](const auto& vL, const auto& vR, auto& res) {
                if (dimL == dimR) { // Division Élément par élément
                    for(size_t i = 0; i < dimL; ++i) {
                        if (vR[i] == 0.0) CFEM_ERROR("Constant folding detected division by zero!");
                        res[i] = vL[i] / vR[i];
                    }
                } else { // Division Multi-dim / Scalaire
                    CfemReal divisor = vR[0];
                    if (divisor == 0.0) CFEM_ERROR("Constant folding detected division by zero!");
                    for(size_t i = 0; i < dimL; ++i) {
                        res[i] = vL[i] / divisor;
                    }
                }
            });
    }

    inline std::shared_ptr<const CfemExpression> transpose(std::shared_ptr<const CfemExpression> u) {
        if (u->getDimension() != 9) return u; 
        if (u->isZero()) return u;

        if (u->isConstant()) {
            auto cU = std::static_pointer_cast<const ConstantExpression>(u);
            const auto& v = cU->getValues();
            return std::make_shared<ConstantExpression>(std::vector<CfemReal>{
                v[0], v[3], v[6],
                v[1], v[4], v[7],
                v[2], v[5], v[8]
            });
        }

        return std::make_shared<UnaryExpression>(u, MathOp::Transpose);
    }

    inline std::shared_ptr<const CfemExpression> dot(std::shared_ptr<const CfemExpression> lhs, 
                                                 std::shared_ptr<const CfemExpression> rhs)
    {
        CfemInt dimL = lhs->getDimension();
        CfemInt dimR = rhs->getDimension();

        // Détermination de la dimension de sortie
        CfemInt outDim = 1; 
        if (dimL == 9 && dimR == 3) outDim = 3;      // Tenseur . Vecteur = Vecteur
        else if (dimL == 3 && dimR == 9) outDim = 3; // Vecteur . Tenseur = Vecteur
        else if (dimL == 9 && dimR == 9) outDim = 9; // Tenseur . Tenseur = Tenseur
        else if (dimL == 3 && dimR == 3) outDim = 1; // Vecteur . Vecteur = Scalaire
        else {
            CFEM_ERROR("Invalid dimensions for dot(): " + std::to_string(dimL) + " . " + std::to_string(dimR));
        }

        // Zero Pruning
        if (lhs->isZero() || rhs->isZero()) {
            if (outDim == 1) return zeros::Zero1D;
            if (outDim == 3) return zeros::Zero3D;
            return zeros::ZeroTensor;
        }

        return foldBinary(lhs, rhs, MathOp::Dot, outDim, 
            [dimL, dimR](const auto& vL, const auto& vR, auto& res) {
                if (dimL == 3 && dimR == 3) { // Scalaire
                    res[0] = vL[0]*vR[0] + vL[1]*vR[1] + vL[2]*vR[2];
                }
                else if (dimL == 9 && dimR == 3) { // Tenseur * Vecteur
                    for (int i = 0; i < 3; ++i)
                        res[i] = vL[i*3+0]*vR[0] + vL[i*3+1]*vR[1] + vL[i*3+2]*vR[2];
                }
                else if (dimL == 3 && dimR == 9) { // Vecteur * Tenseur
                    for (int j = 0; j < 3; ++j)
                        res[j] = vL[0]*vR[0+j] + vL[1]*vR[3+j] + vL[2]*vR[6+j];
                }
                else if (dimL == 9 && dimR == 9) { // MatMul
                    for (int i = 0; i < 3; ++i) {
                        for (int j = 0; j < 3; ++j) {
                            res[i*3+j] = vL[i*3+0]*vR[0*3+j] + vL[i*3+1]*vR[1*3+j] + vL[i*3+2]*vR[2*3+j];
                        }
                    }
                }
            });
    }

    inline std::shared_ptr<const CfemExpression> inner(std::shared_ptr<const CfemExpression> lhs, 
                                                   std::shared_ptr<const CfemExpression> rhs) 
    {
        if (lhs->getDimension() != 9 || rhs->getDimension() != 9) {
            CFEM_ERROR("inner() product (double contraction) requires two dim 9 tensors.");
        }

        if (lhs->isZero() || rhs->isZero()) return zeros::Zero1D;

        return foldBinary(lhs, rhs, MathOp::Inner, 1, 
            [](const auto& vL, const auto& vR, auto& res) {
                CfemReal sum = 0;
                for(int i = 0; i < 9; ++i) sum += vL[i] * vR[i];
                res[0] = sum;
            });
    }

    inline std::shared_ptr<const CfemExpression> cross(std::shared_ptr<const CfemExpression> lhs, 
                                                   std::shared_ptr<const CfemExpression> rhs)
    {
        if (lhs->getDimension() != 3 || rhs->getDimension() != 3) {
            CFEM_ERROR("cross() product is only defined for dim 3 vectors.");
        }

        if (lhs->isZero() || rhs->isZero()) return zeros::Zero3D;

        return foldBinary(lhs, rhs, MathOp::Cross, 3, 
            [](const auto& vL, const auto& vR, auto& res) {
                res[0] = vL[1] * vR[2] - vL[2] * vR[1];
                res[1] = vL[2] * vR[0] - vL[0] * vR[2];
                res[2] = vL[0] * vR[1] - vL[1] * vR[0];
            });
    }

    inline std::shared_ptr<const CfemExpression> atan2(std::shared_ptr<const CfemExpression> y, 
                                                       std::shared_ptr<const CfemExpression> x)
    {
        CfemInt dimY = y->getDimension();
        CfemInt dimX = x->getDimension();

        // Détermination de la dimension (autorise le mix Scalaire/Vecteur)
        if (dimY != dimX && dimY != 1 && dimX != 1) {
            CFEM_ERROR("Invalid dimensions for atan2(): " + std::to_string(dimY) + " , " + std::to_string(dimX));
        }

        CfemInt outDim = std::max(dimY, dimX);

        // On utilise std::atan2 dans le foldBinary
        return foldBinary(y, x, MathOp::Atan2, outDim, 
            [dimY, dimX, outDim](const auto& vY, const auto& vX, auto& res) {
                if (dimY == dimX) { // Élément par élément
                    for(size_t i = 0; i < outDim; ++i) res[i] = std::atan2(vY[i], vX[i]);
                } 
                else if (dimY == 1) { // Scalaire , Vecteur/Tenseur
                    for(size_t i = 0; i < outDim; ++i) res[i] = std::atan2(vY[0], vX[i]);
                } 
                else { // Vecteur/Tenseur , Scalaire
                    for(size_t i = 0; i < outDim; ++i) res[i] = std::atan2(vY[i], vX[0]);
                }
            });
    }

    inline std::shared_ptr<const CfemExpression> pow(std::shared_ptr<const CfemExpression> base, 
                                                     std::shared_ptr<const CfemExpression> exp)
    {
        CfemInt dimBase = base->getDimension();
        CfemInt dimExp  = exp->getDimension();

        if (dimBase != dimExp && dimBase != 1 && dimExp != 1) {
            CFEM_ERROR("Invalid dimensions for pow(): " + std::to_string(dimBase) + " ^ " + std::to_string(dimExp));
        }

        CfemInt outDim = std::max(dimBase, dimExp);

        // Optimisation : base ^ 1 = base
        if (exp->isConstant() && dimExp == 1) {
            auto cExp = std::static_pointer_cast<const ConstantExpression>(exp);
            if (cExp->getValues()[0] == 1.0) {
                return base; 
            }
        }

        // On n'optimise pas agressivement base^0 en ConstantExpression(1.0) 
        // ici pour ne pas devoir instancier un objet Ones3D/OnesTensor, 
        // le foldBinary le fera très bien pour les constantes.

        return foldBinary(base, exp, MathOp::Pow, outDim, 
            [dimBase, dimExp, outDim](const auto& vBase, const auto& vExp, auto& res) {
                if (dimBase == dimExp) { 
                    for(size_t i = 0; i < outDim; ++i) res[i] = std::pow(vBase[i], vExp[i]);
                } 
                else if (dimBase == 1) { 
                    for(size_t i = 0; i < outDim; ++i) res[i] = std::pow(vBase[0], vExp[i]);
                } 
                else { 
                    for(size_t i = 0; i < outDim; ++i) res[i] = std::pow(vBase[i], vExp[0]);
                }
            });
    }

    inline std::shared_ptr<const CfemExpression> abs(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u;
        return foldUnary(u, MathOp::Abs, [](CfemReal x) { return std::abs(x); });
    }

    inline std::shared_ptr<const CfemExpression> sqrt(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u;
        return foldUnary(u, MathOp::Sqrt, [](CfemReal x) { return (x>0) ? std::sqrt(x) : 0.0; });
    }

    inline std::shared_ptr<const CfemExpression> cbrt(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u;
        return foldUnary(u, MathOp::Cbrt, [](CfemReal x) { return std::cbrt(x); });
    }

    inline std::shared_ptr<const CfemExpression> sign(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // sign(0) = 0
        return foldUnary(u, MathOp::Sign, [](CfemReal x) { return (x > 0) ? 1.0 : ((x < 0) ? -1.0 : 0.0); });
    }

    inline std::shared_ptr<const CfemExpression> heaviside(std::shared_ptr<const CfemExpression> u) {
        // H(0) = 1.0 (selon notre convention)
        if (u->isZero()) return makeOnesLike(u);
        return foldUnary(u, MathOp::Heaviside, [](CfemReal x) { return (x >= 0) ? 1.0 : 0.0; });
    }

    inline std::shared_ptr<const CfemExpression> sin(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // sin(0) = 0
        return foldUnary(u, MathOp::Sin, [](CfemReal x) { return std::sin(x); });
    }

    inline std::shared_ptr<const CfemExpression> cos(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return makeOnesLike(u); // cos(0) = 1
        return foldUnary(u, MathOp::Cos, [](CfemReal x) { return std::cos(x); });
    }

    inline std::shared_ptr<const CfemExpression> tan(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // tan(0) = 0
        return foldUnary(u, MathOp::Tan, [](CfemReal x) { return std::tan(x); });
    }

    inline std::shared_ptr<const CfemExpression> cot(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Cot, [](CfemReal x) { return 1.0 / std::tan(x); });
    }

    inline std::shared_ptr<const CfemExpression> sec(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return makeOnesLike(u); // sec(0) = 1/cos(0) = 1
        return foldUnary(u, MathOp::Sec, [](CfemReal x) { return 1.0 / std::cos(x); });
    }

    inline std::shared_ptr<const CfemExpression> csc(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Csc, [](CfemReal x) { return 1.0 / std::sin(x); });
    }

    inline std::shared_ptr<const CfemExpression> asin(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // asin(0) = 0
        return foldUnary(u, MathOp::Asin, [](CfemReal x) { return std::asin(x); });
    }

    inline std::shared_ptr<const CfemExpression> acos(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Acos, [](CfemReal x) { return std::acos(x); });
    }

    inline std::shared_ptr<const CfemExpression> atan(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // atan(0) = 0
        return foldUnary(u, MathOp::Atan, [](CfemReal x) { return std::atan(x); });
    }

    inline std::shared_ptr<const CfemExpression> acot(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Acot, [](CfemReal x) { return M_PI_2 - std::atan(x); });
    }

    inline std::shared_ptr<const CfemExpression> asec(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Asec, [](CfemReal x) { return std::acos(1.0/x); });
    }

    inline std::shared_ptr<const CfemExpression> acsc(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Acsc, [](CfemReal x) { return std::asin(1.0/x); });
    }

    inline std::shared_ptr<const CfemExpression> sinh(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // sinh(0) = 0
        return foldUnary(u, MathOp::Sinh, [](CfemReal x) { return std::sinh(x); });
    }

    inline std::shared_ptr<const CfemExpression> cosh(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return makeOnesLike(u); // cosh(0) = 1
        return foldUnary(u, MathOp::Cosh, [](CfemReal x) { return std::cosh(x); });
    }

    inline std::shared_ptr<const CfemExpression> tanh(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // tanh(0) = 0
        return foldUnary(u, MathOp::Tanh, [](CfemReal x) { return std::tanh(x); });
    }

    inline std::shared_ptr<const CfemExpression> asinh(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // asinh(0) = 0
        return foldUnary(u, MathOp::Asinh, [](CfemReal x) { return std::asinh(x); });
    }

    inline std::shared_ptr<const CfemExpression> acosh(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Acosh, [](CfemReal x) { return std::acosh(x); });
    }

    inline std::shared_ptr<const CfemExpression> atanh(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // atanh(0) = 0
        return foldUnary(u, MathOp::Atanh, [](CfemReal x) { return std::atanh(x); });
    }

    inline std::shared_ptr<const CfemExpression> exp(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return makeOnesLike(u); // exp(0) = 1
        return foldUnary(u, MathOp::Exp, [](CfemReal x) { return std::exp(x); });
    }

    inline std::shared_ptr<const CfemExpression> log(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Log, [](CfemReal x) { return (x>0) ? std::log(x) : -1e20; });
    }

    inline std::shared_ptr<const CfemExpression> erf(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // erf(0) = 0
        return foldUnary(u, MathOp::Erf, [](CfemReal x) { return std::erf(x); });
    }

    inline std::shared_ptr<const CfemExpression> erfc(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return makeOnesLike(u); // erfc(0) = 1
        return foldUnary(u, MathOp::Erfc, [](CfemReal x) { return std::erfc(x); });
    }

    inline std::shared_ptr<const CfemExpression> tgamma(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Gamma, [](CfemReal x) { return std::tgamma(x); });
    }

    inline std::shared_ptr<const CfemExpression> zeta(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Zeta, [](CfemReal x) { return std::riemann_zeta(x); });
    }

    inline std::shared_ptr<const CfemExpression> scalar(const std::function<CfemReal(const Vector3D&)>& f, 
                                                        std::shared_ptr<const CfemExpression> gradient){
        return LambdaExpression::fromScalar(f, gradient);
    }

    inline std::shared_ptr<const CfemExpression> vector(const std::function<Vector3D(const Vector3D&)>& f,
                                                        std::shared_ptr<const CfemExpression> gradient) {
        return LambdaExpression::fromVector(f, gradient);
    }

    inline std::shared_ptr<const CfemExpression> tensor(
                                            const std::function<TensorO2(const Vector3D &)> &f,
                                            std::shared_ptr<const CfemExpression> gradient )
    {
        return LambdaExpression::fromTensor(f, gradient);
    }

    inline std::shared_ptr<const CfemExpression> grad(const std::shared_ptr<const CfemExpression>& e) {
        if (!e) return nullptr;

        if (auto lamb = std::dynamic_pointer_cast<const LambdaExpression>(e)) {
            auto analyticalGrad = lamb->grad();
            if (analyticalGrad) {
                return analyticalGrad;
            }
        }

        if (auto symScalar = std::dynamic_pointer_cast<const ScalarSymbolicExpression>(e)) {
            return std::make_shared<VectorSymbolicExpression>(
                symScalar->diff(x_rcp), symScalar->diff(y_rcp), symScalar->diff(z_rcp)
            );
        }
        
        if (auto vect = std::dynamic_pointer_cast<const VectorSymbolicExpression>(e)) {
                        
            return std::make_shared<TensorSymbolicExpression>(
                vect->getComponent(0)->diff(x_rcp), vect->getComponent(0)->diff(y_rcp), vect->getComponent(0)->diff(z_rcp), 
                vect->getComponent(1)->diff(x_rcp), vect->getComponent(1)->diff(y_rcp), vect->getComponent(1)->diff(z_rcp), 
                vect->getComponent(2)->diff(x_rcp), vect->getComponent(2)->diff(y_rcp), vect->getComponent(2)->diff(z_rcp)
            );

        }

        if (auto tensor = std::dynamic_pointer_cast<const TensorSymbolicExpression>(e)) {
                        
            // return std::make_shared<TensorO3SymbolicExpression>(
            //     tensor->getComponent(0)->diff(x_rcp), tensor->getComponent(0)->diff(y_rcp), tensor->getComponent(0)->diff(z_rcp), 
            //     tensor->getComponent(1)->diff(x_rcp), tensor->getComponent(1)->diff(y_rcp), tensor->getComponent(1)->diff(z_rcp), 
            //     tensor->getComponent(2)->diff(x_rcp), tensor->getComponent(2)->diff(y_rcp), tensor->getComponent(2)->diff(z_rcp),
            //     tensor->getComponent(3)->diff(x_rcp), tensor->getComponent(3)->diff(y_rcp), tensor->getComponent(3)->diff(z_rcp), 
            //     tensor->getComponent(4)->diff(x_rcp), tensor->getComponent(4)->diff(y_rcp), tensor->getComponent(4)->diff(z_rcp), 
            //     tensor->getComponent(5)->diff(x_rcp), tensor->getComponent(5)->diff(y_rcp), tensor->getComponent(5)->diff(z_rcp),
            //     tensor->getComponent(6)->diff(x_rcp), tensor->getComponent(6)->diff(y_rcp), tensor->getComponent(6)->diff(z_rcp), 
            //     tensor->getComponent(7)->diff(x_rcp), tensor->getComponent(7)->diff(y_rcp), tensor->getComponent(7)->diff(z_rcp), 
            //     tensor->getComponent(8)->diff(x_rcp), tensor->getComponent(8)->diff(y_rcp), tensor->getComponent(8)->diff(z_rcp)
            // );
            CFEM_ERROR("TODO: Gradient of a second order tensor currentely not supported. Please implement the class TensorO3SymbolicExpression.");

        }

        if (auto bin = std::dynamic_pointer_cast<const BinaryExpression>(e)) {
            auto op = bin->getOp();
            if (op == MathOp::Add || op == MathOp::Sub) {
                auto gradL = grad(bin->getLHS());
                auto gradR = grad(bin->getRHS());
                
                if (op == MathOp::Add) return gradL + gradR;
                else return gradL - gradR;
            }
            
            // Product Rule: grad(A*B) = grad(A)*B + A*grad(B)
            if (op == MathOp::Mul) {
                auto A = bin->getLHS();
                auto B = bin->getRHS();
                // (grad A * B) + (A * grad B)
                return (grad(A) * B) + (A * grad(B));
            }

            // TODO : Division by a scalar
        }

        // 4. NEW: Magnitude Squared grad(|u|^2) -> 2 * u * grad(u)
        if (auto magSq = std::dynamic_pointer_cast<const MagnitudeSqExpression>(e)) {
            auto inner = magSq->getInner(); 
            std::vector<CfemReal> lval = {2.0};
            auto two = std::make_shared<ConstantExpression>(lval); // Assuming 1D constant helper
            
            // Note: If 'inner' is a vector, 'grad(inner)' is a tensor.
            // The result of this multiplication must be a vector.
            return two * inner * grad(inner);
        }

        // 5. Fallback for Leaf Nodes (FEField, etc.)
        // This creates a GradExpression which will be evaluated numerically
        return std::make_shared<GradExpression>(e);
    }

    inline std::shared_ptr<const CfemExpression> div(const std::shared_ptr<const CfemExpression>& e) {
        if (auto v = std::dynamic_pointer_cast<const VectorSymbolicExpression>(e)) {
            auto scalar_comp = std::static_pointer_cast<const ScalarSymbolicExpression>(v->getComponent(0));
            auto res = v->getComponent(0)->get_expression().diff(x_rcp) + 
                       v->getComponent(1)->get_expression().diff(y_rcp) + 
                       v->getComponent(2)->get_expression().diff(z_rcp);
            return scalar(res);
        }

        if (auto t = std::dynamic_pointer_cast<const TensorSymbolicExpression>(e)) {
            auto r0 = t->getComponent(0)->get_expression().diff(x_rcp) + 
                      t->getComponent(1)->get_expression().diff(y_rcp) + 
                      t->getComponent(2)->get_expression().diff(z_rcp);
            auto r1 = t->getComponent(3)->get_expression().diff(x_rcp) + 
                      t->getComponent(4)->get_expression().diff(y_rcp) + 
                      t->getComponent(5)->get_expression().diff(z_rcp);
            auto r2 = t->getComponent(6)->get_expression().diff(x_rcp) + 
                      t->getComponent(7)->get_expression().diff(y_rcp) + 
                      t->getComponent(8)->get_expression().diff(z_rcp);

            
            // TODO : Etendre aux lambda
            return std::make_shared<VectorSymbolicExpression>(scalar(r0), scalar(r1), scalar(r2));
        }

        

        CFEM_ERROR("div() only supports Symbolic Vector/Tensor types currently.");
        return nullptr;
    }

    inline std::shared_ptr<const CfemExpression> mag_sqr(std::shared_ptr<const CfemExpression> a) {
        return std::make_shared<MagnitudeSqExpression>(std::move(a));
    }

    inline std::shared_ptr<const CfemExpression> scalar(const Expression& ex) {
        return std::make_shared<ScalarSymbolicExpression>(ex, std::vector{x_rcp, y_rcp, z_rcp});
    }

    inline std::shared_ptr<const CfemExpression> vector(const Expression& ex, 
                                                        const Expression& ey, 
                                                        const Expression& ez) {
        return std::make_shared<VectorSymbolicExpression>(scalar(ex), scalar(ey), scalar(ez));
    }


    // SymEngine:
    /**
     * @name Common math functions
     *  @{
     */
    inline Expression sin(const Expression &e) { return SymEngine::sin(e); }
    inline Expression cos(const Expression &e) { return SymEngine::cos(e); }
    inline Expression tan(const Expression &e) { return SymEngine::tan(e); }
    inline Expression cot(const Expression &e) { return SymEngine::cot(e); }
    inline Expression sec(const Expression &e) { return SymEngine::sec(e); }
    inline Expression csc(const Expression &e) { return SymEngine::csc(e); }

    inline Expression asin(const Expression &e) { return SymEngine::asin(e); }
    inline Expression acos(const Expression &e) { return SymEngine::acos(e); }
    inline Expression atan(const Expression &e) { return SymEngine::atan(e); }
    inline Expression acot(const Expression &e) { return SymEngine::acot(e); }
    inline Expression asec(const Expression &e) { return SymEngine::asec(e); }
    inline Expression acsc(const Expression &e) { return SymEngine::acsc(e); }

    inline Expression sinh(const Expression &e) { return SymEngine::sinh(e); }
    inline Expression cosh(const Expression &e) { return SymEngine::cosh(e); }
    inline Expression tanh(const Expression &e) { return SymEngine::tanh(e); }
    inline Expression coth(const Expression &e) { return SymEngine::coth(e); }
    inline Expression asinh(const Expression &e) { return SymEngine::asinh(e); }
    inline Expression acosh(const Expression &e) { return SymEngine::acosh(e); }
    inline Expression atanh(const Expression &e) { return SymEngine::atanh(e); }

    inline Expression exp(const Expression &e) { return SymEngine::exp(e); }
    inline Expression log(const Expression &e) { return SymEngine::log(e); }
    inline Expression ln(const Expression &e) { return SymEngine::log(e); }

    inline Expression sqrt(const Expression &e) { return SymEngine::sqrt(e); }
    inline Expression cbrt(const Expression &e) { return SymEngine::cbrt(e); }
    inline Expression pow(const Expression &b, const Expression &e) { return SymEngine::pow(b, e); }

    inline Expression abs(const Expression &e) { return SymEngine::abs(e); }

    inline Expression max(const Expression &a, const Expression &b) { return SymEngine::max({a, b}); }
    inline Expression min(const Expression &a, const Expression &b) { return SymEngine::min({a, b}); }

    inline Expression gamma(const Expression &e) { return SymEngine::gamma(e); }
    inline Expression erf(const Expression &e) { return SymEngine::erf(e); }
    inline Expression erfc(const Expression &e) { return SymEngine::erfc(e); }
    inline Expression zeta(const Expression &e) { return SymEngine::zeta(e); }
    /** @} */

}