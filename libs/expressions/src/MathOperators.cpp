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

#include "AlgebraicSimplifier.h"
#include "BinaryExpression.h"
#include "CachedExpression.h"
#include "ConstantExpression.h"
#include "CoordinateExpression.h"
#include "GeometricNormalExpression.h"
#include "GradExpression.h"
#include "LambdaExpression.h"
#include "MagnitudeSqExpression.h"
#include "MathOperators.h"
#include "ScaledExpression.h"
#include "symengine/ScalarSymbolicExpression.h"
#include "symengine/SymEngineSymbols.h"
#include "symengine/TensorSymbolicExpression.h"
#include "symengine/VectorSymbolicExpression.h"
#include "TensorExpression.h"
#include "UnaryExpression.h"
#include "VariableExpression.h"
#include "VectorExpression.h"
#include "ZeroExpression.h"

using namespace cfem::sym::utils;

namespace cfem::sym
{
    std::shared_ptr<const CfemExpression> operator-(std::shared_ptr<const CfemExpression> u) 
    {
        if (u->isZero()) return u;

        // 2. Constant folding check
        if (auto cU = std::dynamic_pointer_cast<const ConstantExpression>(u)) {
            const auto& v = cU->getValues();
            std::vector<CfemReal> res(v.size());
            for(size_t i = 0; i < v.size(); ++i) res[i] = -v[i];
            return std::make_shared<ConstantExpression>(std::move(res));
        }

        if (auto sym = std::dynamic_pointer_cast<const ScalarSymbolicExpression>(u)) {
            return std::make_shared<ScalarSymbolicExpression>( -sym->get_expression(), sym->get_vars() );
        }

        if (auto sym = std::dynamic_pointer_cast<const VectorSymbolicExpression>(u)) {
            auto exp_x = std::make_shared<ScalarSymbolicExpression>( -sym->getComponent(0)->get_expression(), sym->getComponent(0)->get_vars() );
            auto exp_y = std::make_shared<ScalarSymbolicExpression>( -sym->getComponent(1)->get_expression(), sym->getComponent(1)->get_vars() );
            auto exp_z = std::make_shared<ScalarSymbolicExpression>( -sym->getComponent(2)->get_expression(), sym->getComponent(2)->get_vars() );
            return std::make_shared<VectorSymbolicExpression>(exp_x, exp_y, exp_z);
        }

        if (auto sym = std::dynamic_pointer_cast<const TensorSymbolicExpression>(u)) {
            auto exp_xx = std::make_shared<ScalarSymbolicExpression>( -sym->getComponent(0)->get_expression(), sym->getComponent(0)->get_vars() );
            auto exp_xy = std::make_shared<ScalarSymbolicExpression>( -sym->getComponent(1)->get_expression(), sym->getComponent(1)->get_vars() );
            auto exp_xz = std::make_shared<ScalarSymbolicExpression>( -sym->getComponent(2)->get_expression(), sym->getComponent(2)->get_vars() );
            auto exp_yx = std::make_shared<ScalarSymbolicExpression>( -sym->getComponent(3)->get_expression(), sym->getComponent(3)->get_vars() );
            auto exp_yy = std::make_shared<ScalarSymbolicExpression>( -sym->getComponent(4)->get_expression(), sym->getComponent(4)->get_vars() );
            auto exp_yz = std::make_shared<ScalarSymbolicExpression>( -sym->getComponent(5)->get_expression(), sym->getComponent(5)->get_vars() );
            auto exp_zx = std::make_shared<ScalarSymbolicExpression>( -sym->getComponent(6)->get_expression(), sym->getComponent(6)->get_vars() );
            auto exp_zy = std::make_shared<ScalarSymbolicExpression>( -sym->getComponent(7)->get_expression(), sym->getComponent(7)->get_vars() );
            auto exp_zz = std::make_shared<ScalarSymbolicExpression>( -sym->getComponent(8)->get_expression(), sym->getComponent(8)->get_vars() );
            return std::make_shared<TensorSymbolicExpression>(exp_xx, exp_xy, exp_xz,
                                                              exp_yx, exp_yy, exp_yz,
                                                              exp_zx, exp_zy, exp_zz);
        }
        
        // 3. Smart Scaling: - (k * u) -> (-k) * u
        if (auto uScaled = std::dynamic_pointer_cast<const ScaledExpression>(u)) {
            return std::make_shared<ScaledExpression>(-uScaled->getCoeff(), uScaled->getBaseShared());
        }

        // 4. Default: Treat -u as (-1.0) * u
        return std::make_shared<ScaledExpression>(-1.0, u);
    }

    std::shared_ptr<const CfemExpression> operator-(std::shared_ptr<const CfemExpression> lhs, 
                                                    std::shared_ptr<const CfemExpression> rhs) 
    {
        CFEM_ASSERT(lhs && rhs, "");
        return lhs + (-rhs); 
    }

    std::shared_ptr<const CfemExpression> operator+(std::shared_ptr<const CfemExpression> lhs, 
                                                    std::shared_ptr<const CfemExpression> rhs)
    {
        if (lhs->getNumComponents() != rhs->getNumComponents())
            CFEM_ERROR << "Dimension mismatch in operator << ";

        if (lhs->isZero()) return rhs;
        if (rhs->isZero()) return lhs;


        if (lhs->getType() == ExpressionType::ScalarSymbolic && rhs->getType() == ExpressionType::ScalarSymbolic){
            auto L = std::dynamic_pointer_cast<const ScalarSymbolicExpression>(lhs);
            auto R = std::dynamic_pointer_cast<const ScalarSymbolicExpression>(rhs);
            auto vars = L->get_vars();
            if (vars != R->get_vars()) CFEM_ERROR << "Symbolic Expressions must depend on the same list of variables";

            return std::make_shared<ScalarSymbolicExpression>( L->get_expression() + R->get_expression(), vars);
        }

        if (lhs->getType() == ExpressionType::VectorSymbolic && rhs->getType() == ExpressionType::VectorSymbolic){
            auto L = std::dynamic_pointer_cast<const VectorSymbolicExpression>(lhs);
            auto R = std::dynamic_pointer_cast<const VectorSymbolicExpression>(rhs);

            auto vars = L->getComponent(0)->get_vars();
            if (vars != R->getComponent(0)->get_vars()) CFEM_ERROR << "Symbolic Expressions must depend on the same list of variables";

            return std::make_shared<VectorSymbolicExpression>( std::dynamic_pointer_cast<const ScalarSymbolicExpression>( L->getComponent(0) + R->getComponent(0) ),
                                                               std::dynamic_pointer_cast<const ScalarSymbolicExpression>( L->getComponent(1) + R->getComponent(1) ),
                                                               std::dynamic_pointer_cast<const ScalarSymbolicExpression>( L->getComponent(2) + R->getComponent(2) ) );
        }

        if (lhs->getType() == ExpressionType::TensorSymbolic && rhs->getType() == ExpressionType::TensorSymbolic){
            auto L = std::dynamic_pointer_cast<const TensorSymbolicExpression>(lhs);
            auto R = std::dynamic_pointer_cast<const TensorSymbolicExpression>(rhs);

            auto vars = L->getComponent(0)->get_vars();
            if (vars != R->getComponent(0)->get_vars()) CFEM_ERROR << "Symbolic Expressions must depend on the same list of variables";

            return std::make_shared<TensorSymbolicExpression>( std::dynamic_pointer_cast<const ScalarSymbolicExpression>( L->getComponent(0) + R->getComponent(0) ),
                                                               std::dynamic_pointer_cast<const ScalarSymbolicExpression>( L->getComponent(1) + R->getComponent(1) ),
                                                               std::dynamic_pointer_cast<const ScalarSymbolicExpression>( L->getComponent(2) + R->getComponent(2) ),
                                                               std::dynamic_pointer_cast<const ScalarSymbolicExpression>( L->getComponent(3) + R->getComponent(3) ),
                                                               std::dynamic_pointer_cast<const ScalarSymbolicExpression>( L->getComponent(4) + R->getComponent(4) ),
                                                               std::dynamic_pointer_cast<const ScalarSymbolicExpression>( L->getComponent(5) + R->getComponent(5) ),
                                                               std::dynamic_pointer_cast<const ScalarSymbolicExpression>( L->getComponent(6) + R->getComponent(6) ),
                                                               std::dynamic_pointer_cast<const ScalarSymbolicExpression>( L->getComponent(7) + R->getComponent(7) ),
                                                               std::dynamic_pointer_cast<const ScalarSymbolicExpression>( L->getComponent(8) + R->getComponent(8) ) );
        }

        return simplifyAdd(lhs, rhs);
    }

    std::shared_ptr<const CfemExpression> operator*(std::shared_ptr<const CfemExpression> lhs, 
                                                    std::shared_ptr<const CfemExpression> rhs) 
    {
        CfemInt dimL = lhs->getNumComponents();
        CfemInt dimR = rhs->getNumComponents();

        CfemInt outDim = std::max(dimL, dimR);

        // Cas Scalaire (L ou R est dimension 1)
        if (dimL == 1 || dimR == 1) {
            outDim = std::max(dimL, dimR);
        }
        // Cas Tenseur * Tenseur (MatMul : 9 * 9 -> 9)
        else if (dimL == 9 && dimR == 9) {
            outDim = 9;
        }
        // Cas Tenseur * Vecteur (9 * 3 -> 3)
        else if (dimL == 9 && dimR == 3) {
            outDim = 3;
        }
        // Cas Vecteur * Tenseur (3 * 9 -> 3) - optionnel selon tes besoins
        // else if (dimL == 3 && dimR == 9) {
        //     outDim = 3;
        // }
        // Cas Hadamard ou Erreur (ex: 3 * 3 -> 3 si tu veux du terme à terme)
        else if (dimL == dimR) {
            outDim = dimL;
        }
        else {
            CFEM_ERROR << "Invalid operator* : Cannot multiply dimension " << dimL << " with " << dimR;
        }

        if (lhs->isZero() || rhs->isZero()) {
            if (outDim == 1) return zeros::Zero1D;
            if (outDim == 3) return zeros::Zero3D;
            if (outDim == 9) return zeros::ZeroTensor;
            return std::make_shared<ZeroExpression>(outDim);
        }

        if (lhs->isConstant() && dimL == 1) {
            CfemReal k = std::static_pointer_cast<const ConstantExpression>(lhs)->getValues()[0];
            if (k == 1.0) return rhs;
            
            // If we are scaling something already scaled: k1 * (k2 * u) -> (k1*k2) * u
            if (auto rScaled = std::dynamic_pointer_cast<const ScaledExpression>(rhs)) {
                return std::make_shared<ScaledExpression>(k * rScaled->getCoeff(), rScaled->getBaseShared());
            }
            return std::make_shared<ScaledExpression>(k, rhs);
        }

        if (rhs->isConstant() && dimR == 1) {
            CfemReal k = std::static_pointer_cast<const ConstantExpression>(rhs)->getValues()[0];
            if (k == 1.0) return lhs;
            if (auto lScaled = std::dynamic_pointer_cast<const ScaledExpression>(lhs)) {
                return std::make_shared<ScaledExpression>(k * lScaled->getCoeff(), lScaled->getBaseShared());
            }
            return std::make_shared<ScaledExpression>(k, lhs);
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

    
    std::shared_ptr<const CfemExpression> operator/(std::shared_ptr<const CfemExpression> lhs, 
                                                    std::shared_ptr<const CfemExpression> rhs) 
    {
        CfemInt dimL = lhs->getNumComponents();
        CfemInt dimR = rhs->getNumComponents();

        if (dimL != dimR && dimR != 1) {
            CFEM_ERROR << "Invalid operator/: Divisor must be a scalar or match numerator dimension. Got " 
                    << std::to_string(dimL) << " / " << std::to_string(dimR);
        }

        if (rhs->isZero()) {
            CFEM_ERROR << "Division by Zero detected in symbolic expression setup!";
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
                        if (vR[i] == 0.0) CFEM_ERROR << "Constant folding detected division by zero!";
                        res[i] = vL[i] / vR[i];
                    }
                } else { // Division Multi-dim / Scalaire
                    CfemReal divisor = vR[0];
                    if (divisor == 0.0) CFEM_ERROR << "Constant folding detected division by zero!";
                    for(size_t i = 0; i < dimL; ++i) {
                        res[i] = vL[i] / divisor;
                    }
                }
            });
    }

    std::shared_ptr<const CfemExpression> identityO2TensorExpr(){
        return std::make_shared<const ConstantExpression>(IdentityTensor3D());
    }


    std::shared_ptr<const CfemExpression> transpose(std::shared_ptr<const CfemExpression> u) {
        if (u->getNumComponents() != 9) CFEM_ERROR << "`transpose()` expects a second order tensor argument."; 
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


    std::shared_ptr<const CfemExpression> trace(std::shared_ptr<const CfemExpression> u){
        if (u->getNumComponents() != 9) CFEM_ERROR << "`trace()` expects a second order tensor argument."; 
        if (u->isZero()) return std::make_shared<ZeroExpression>(0.0);

        if (u->isConstant()) {
            auto cU = std::static_pointer_cast<const ConstantExpression>(u);
            const auto& v = cU->getValues();
            return std::make_shared<ConstantExpression>(v[0] + v[4] + v[8]);
        }

        return std::make_shared<UnaryExpression>(u, MathOp::Trace);
    }


    std::shared_ptr<const CfemExpression> dot(std::shared_ptr<const CfemExpression> lhs, 
                                                     std::shared_ptr<const CfemExpression> rhs)
    {
        CfemInt dimL = lhs->getNumComponents();
        CfemInt dimR = rhs->getNumComponents();

        // Détermination de la dimension de sortie
        CfemInt outDim = 1; 
        if (dimL == 9 && dimR == 3) outDim = 3;      // Tenseur . Vecteur = Vecteur
        else if (dimL == 3 && dimR == 9) outDim = 3; // Vecteur . Tenseur = Vecteur
        else if (dimL == 9 && dimR == 9) outDim = 9; // Tenseur . Tenseur = Tenseur
        else if (dimL == 3 && dimR == 3) outDim = 1; // Vecteur . Vecteur = Scalaire
        else {
            CFEM_ERROR << "Invalid dimensions for dot(): " << std::to_string(dimL) << " . " << std::to_string(dimR);
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
                    for (CfemInt i = 0; i < 3; ++i)
                        res[i] = vL[i*3+0]*vR[0] + vL[i*3+1]*vR[1] + vL[i*3+2]*vR[2];
                }
                else if (dimL == 3 && dimR == 9) { // Vecteur * Tenseur
                    for (CfemInt j = 0; j < 3; ++j)
                        res[j] = vL[0]*vR[0+j] + vL[1]*vR[3+j] + vL[2]*vR[6+j];
                }
                else if (dimL == 9 && dimR == 9) { // MatMul
                    for (CfemInt i = 0; i < 3; ++i) {
                        for (CfemInt j = 0; j < 3; ++j) {
                            res[i*3+j] = vL[i*3+0]*vR[0*3+j] + vL[i*3+1]*vR[1*3+j] + vL[i*3+2]*vR[2*3+j];
                        }
                    }
                }
            });
    }

    std::shared_ptr<const CfemExpression> inner(std::shared_ptr<const CfemExpression> lhs, 
                                                std::shared_ptr<const CfemExpression> rhs) 
    {
        if (lhs->getNumComponents() != 9 || rhs->getNumComponents() != 9) {
            CFEM_ERROR << "inner() product (double contraction) requires two dim 9 tensors.";
        }

        if (lhs->isZero() || rhs->isZero()) return zeros::Zero1D;

        return foldBinary(lhs, rhs, MathOp::Inner, 1, 
            [](const auto& vL, const auto& vR, auto& res) {
                CfemReal sum = 0;
                for(CfemInt i = 0; i < 9; ++i) sum += vL[i] * vR[i];
                res[0] = sum;
            });
    }

    std::shared_ptr<const CfemExpression> cross(std::shared_ptr<const CfemExpression> lhs, 
                                                std::shared_ptr<const CfemExpression> rhs)
    {
        if (lhs->getNumComponents() != 3 || rhs->getNumComponents() != 3) {
            CFEM_ERROR << "cross() product is only defined for dim 3 vectors.";
        }

        if (lhs->isZero() || rhs->isZero()) return zeros::Zero3D;

        return foldBinary(lhs, rhs, MathOp::Cross, 3, 
            [](const auto& vL, const auto& vR, auto& res) {
                res[0] = vL[1] * vR[2] - vL[2] * vR[1];
                res[1] = vL[2] * vR[0] - vL[0] * vR[2];
                res[2] = vL[0] * vR[1] - vL[1] * vR[0];
            });
    }

    std::shared_ptr<const CfemExpression> outer(std::shared_ptr<const CfemExpression> lhs, 
                                                std::shared_ptr<const CfemExpression> rhs)
    {
        if (lhs->getNumComponents() != 3 || rhs->getNumComponents() != 3) {
            CFEM_ERROR << "outer() product is only defined for dim 3 vectors.";
        }

        if (lhs->isZero() || rhs->isZero()) return zeros::ZeroTensor;

        return foldBinary(lhs, rhs, MathOp::Outer, 9, 
            [](const auto& vL, const auto& vR, auto& res) {
                for (CfemInt i = 0; i < 3; ++i) {
                    for (CfemInt j = 0; j < 3; ++j) {
                        res[i*3+j] = vL[i]*vR[j];
                    }
                }
            });
    }

    std::shared_ptr<const CfemExpression> atan2(std::shared_ptr<const CfemExpression> y, 
                                                std::shared_ptr<const CfemExpression> x)
    {
        CfemInt dimY = y->getNumComponents();
        CfemInt dimX = x->getNumComponents();

        // Détermination de la dimension (autorise le mix Scalaire/Vecteur)
        if (dimY != dimX && dimY != 1 && dimX != 1) {
            CFEM_ERROR << "Invalid dimensions for atan2(): " << std::to_string(dimY) << " , " << std::to_string(dimX);
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

    std::shared_ptr<const CfemExpression> pow(std::shared_ptr<const CfemExpression> base, 
                                              std::shared_ptr<const CfemExpression> exp)
    {
        CfemInt dimBase = base->getNumComponents();
        CfemInt dimExp  = exp->getNumComponents();

        if (dimBase != dimExp && dimBase != 1 && dimExp != 1) {
            CFEM_ERROR << "Invalid dimensions for pow(): " << std::to_string(dimBase) << " ^ " << std::to_string(dimExp);
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

    // Fonction utilitaire pour générer un tenseur/vecteur de "1.0"
    std::shared_ptr<const CfemExpression> makeOnesLike(std::shared_ptr<const CfemExpression> u) {
        return std::make_shared<ConstantExpression>(std::vector<CfemReal>(u->getNumComponents(), 1.0));
    }

    
    std::shared_ptr<const CfemExpression> abs(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u;
        return foldUnary(u, MathOp::Abs, [](CfemReal x) { return std::abs(x); });
    }

    std::shared_ptr<const CfemExpression> sqrt(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u;
        return foldUnary(u, MathOp::Sqrt, [](CfemReal x) { return (x>0) ? std::sqrt(x) : 0.0; });
    }

    std::shared_ptr<const CfemExpression> cbrt(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u;
        return foldUnary(u, MathOp::Cbrt, [](CfemReal x) { return std::cbrt(x); });
    }

    std::shared_ptr<const CfemExpression> sign(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // sign(0) = 0
        return foldUnary(u, MathOp::Sign, [](CfemReal x) { return (x > 0) ? 1.0 : ((x < 0) ? -1.0 : 0.0); });
    }

    std::shared_ptr<const CfemExpression> heaviside(std::shared_ptr<const CfemExpression> u) {
        // H(0) = 1.0 (selon notre convention)
        if (u->isZero()) return makeOnesLike(u);
        return foldUnary(u, MathOp::Heaviside, [](CfemReal x) { return (x >= 0) ? 1.0 : 0.0; });
    }

    std::shared_ptr<const CfemExpression> sin(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // sin(0) = 0
        return foldUnary(u, MathOp::Sin, [](CfemReal x) { return std::sin(x); });
    }

    std::shared_ptr<const CfemExpression> cos(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return makeOnesLike(u); // cos(0) = 1
        return foldUnary(u, MathOp::Cos, [](CfemReal x) { return std::cos(x); });
    }

    std::shared_ptr<const CfemExpression> tan(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // tan(0) = 0
        return foldUnary(u, MathOp::Tan, [](CfemReal x) { return std::tan(x); });
    }

    std::shared_ptr<const CfemExpression> cot(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Cot, [](CfemReal x) { return 1.0 / std::tan(x); });
    }

    std::shared_ptr<const CfemExpression> sec(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return makeOnesLike(u); // sec(0) = 1/cos(0) = 1
        return foldUnary(u, MathOp::Sec, [](CfemReal x) { return 1.0 / std::cos(x); });
    }

    std::shared_ptr<const CfemExpression> csc(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Csc, [](CfemReal x) { return 1.0 / std::sin(x); });
    }

    std::shared_ptr<const CfemExpression> asin(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // asin(0) = 0
        return foldUnary(u, MathOp::Asin, [](CfemReal x) { return std::asin(x); });
    }

    std::shared_ptr<const CfemExpression> acos(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Acos, [](CfemReal x) { return std::acos(x); });
    }

    std::shared_ptr<const CfemExpression> atan(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // atan(0) = 0
        return foldUnary(u, MathOp::Atan, [](CfemReal x) { return std::atan(x); });
    }

    std::shared_ptr<const CfemExpression> acot(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Acot, [](CfemReal x) { return M_PI_2 - std::atan(x); });
    }

    std::shared_ptr<const CfemExpression> asec(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Asec, [](CfemReal x) { return std::acos(1.0/x); });
    }

    std::shared_ptr<const CfemExpression> acsc(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Acsc, [](CfemReal x) { return std::asin(1.0/x); });
    }

    std::shared_ptr<const CfemExpression> sinh(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // sinh(0) = 0
        return foldUnary(u, MathOp::Sinh, [](CfemReal x) { return std::sinh(x); });
    }

    std::shared_ptr<const CfemExpression> cosh(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return makeOnesLike(u); // cosh(0) = 1
        return foldUnary(u, MathOp::Cosh, [](CfemReal x) { return std::cosh(x); });
    }

    std::shared_ptr<const CfemExpression> tanh(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // tanh(0) = 0
        return foldUnary(u, MathOp::Tanh, [](CfemReal x) { return std::tanh(x); });
    }

    std::shared_ptr<const CfemExpression> asinh(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // asinh(0) = 0
        return foldUnary(u, MathOp::Asinh, [](CfemReal x) { return std::asinh(x); });
    }

    std::shared_ptr<const CfemExpression> acosh(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Acosh, [](CfemReal x) { return std::acosh(x); });
    }

    std::shared_ptr<const CfemExpression> atanh(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // atanh(0) = 0
        return foldUnary(u, MathOp::Atanh, [](CfemReal x) { return std::atanh(x); });
    }

    std::shared_ptr<const CfemExpression> exp(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return makeOnesLike(u); // exp(0) = 1
        return foldUnary(u, MathOp::Exp, [](CfemReal x) { return std::exp(x); });
    }

    std::shared_ptr<const CfemExpression> log(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Log, [](CfemReal x) { return (x>0) ? std::log(x) : -1e20; });
    }

    std::shared_ptr<const CfemExpression> erf(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return u; // erf(0) = 0
        return foldUnary(u, MathOp::Erf, [](CfemReal x) { return std::erf(x); });
    }

    std::shared_ptr<const CfemExpression> erfc(std::shared_ptr<const CfemExpression> u) {
        if (u->isZero()) return makeOnesLike(u); // erfc(0) = 1
        return foldUnary(u, MathOp::Erfc, [](CfemReal x) { return std::erfc(x); });
    }

    std::shared_ptr<const CfemExpression> tgamma(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Gamma, [](CfemReal x) { return std::tgamma(x); });
    }

    std::shared_ptr<const CfemExpression> zeta(std::shared_ptr<const CfemExpression> u) {
        return foldUnary(u, MathOp::Zeta, [](CfemReal x) { return std::riemann_zeta(x); });
    }

    std::shared_ptr<const CfemExpression> constant(CfemReal c)
    {
        return std::make_shared<const ConstantExpression>(c);
    }

    std::shared_ptr<const CfemExpression> constant(const Vector3D& v)
    {
        return std::make_shared<const ConstantExpression>(v);
    }

    std::shared_ptr<const CfemExpression> constant(const Tensor3D& t)
    {
        return std::make_shared<const ConstantExpression>(t);
    }

    std::shared_ptr<const CfemExpression> vector(std::shared_ptr<const CfemExpression> expr)
    {
        if (expr->getNumComponents() != 3 ) CFEM_ERROR << "expr must be a vector (dimesion 3) expression.";
        return expr;
    }

    std::shared_ptr<const CfemExpression> vector(std::shared_ptr<const CfemExpression> expr_x, 
                                                 std::shared_ptr<const CfemExpression> expr_y,
                                                 std::shared_ptr<const CfemExpression> expr_z)
    {
        if (expr_x->getNumComponents() != 1 || expr_y->getNumComponents() != 1 || expr_z->getNumComponents() != 1 )
            CFEM_ERROR << "expr_x, expr_y and expr_z must be scalar (dimesion 1) expressions.";

        return std::make_shared<const VectorExpression>(expr_x, expr_y, expr_z);
    }

    std::shared_ptr<const CfemExpression> tensor(std::shared_ptr<const CfemExpression> expr)
    {
        if (expr->getNumComponents() != 9 ) CFEM_ERROR << "expr must be a scalar expression.";
        return expr;
    }

    std::shared_ptr<const CfemExpression> tensor(std::shared_ptr<const CfemExpression> expr_xx,
                                                 std::shared_ptr<const CfemExpression> expr_xy,
                                                 std::shared_ptr<const CfemExpression> expr_xz,
                                                 std::shared_ptr<const CfemExpression> expr_yx,
                                                 std::shared_ptr<const CfemExpression> expr_yy,
                                                 std::shared_ptr<const CfemExpression> expr_yz,
                                                 std::shared_ptr<const CfemExpression> expr_zx,
                                                 std::shared_ptr<const CfemExpression> expr_zy,
                                                 std::shared_ptr<const CfemExpression> expr_zz)
    {
        return std::make_shared<const TensorExpression>(std::move(expr_xx), std::move(expr_xy), std::move(expr_xz), 
                                                        std::move(expr_yx), std::move(expr_yy), std::move(expr_yz), 
                                                        std::move(expr_zx), std::move(expr_zy), std::move(expr_zz));
    }

    std::shared_ptr<const CfemExpression> symmetricTensor(std::shared_ptr<const CfemExpression> xx, 
                                                          std::shared_ptr<const CfemExpression> yy, 
                                                          std::shared_ptr<const CfemExpression> zz,
                                                          std::shared_ptr<const CfemExpression> yz, 
                                                          std::shared_ptr<const CfemExpression> xz, 
                                                          std::shared_ptr<const CfemExpression> xy) 
    {
        return std::make_shared<const TensorExpression>(xx, xy, xz, 
                                                        xy, yy, yz, 
                                                        xz, yz, zz);
    }

    /**
     * @brief Creates a diagonal tensor expression: diag(a, b, c)
     */
    std::shared_ptr<const CfemExpression> diag(std::shared_ptr<const CfemExpression> a,
                                               std::shared_ptr<const CfemExpression> b,
                                               std::shared_ptr<const CfemExpression> c)
    {
        return std::make_shared<const TensorExpression>(
            a,    zeros::Zero1D, zeros::Zero1D,
            zeros::Zero1D, b,    zeros::Zero1D,
            zeros::Zero1D, zeros::Zero1D, c
        );
    }

    std::shared_ptr<const CfemExpression> diag(CfemReal a,
                                                      CfemReal b,
                                                      CfemReal c)
    {
        return std::make_shared<const TensorExpression>(
            constant(a),    zeros::Zero1D, zeros::Zero1D,
            zeros::Zero1D, constant(b),    zeros::Zero1D,
            zeros::Zero1D, zeros::Zero1D, constant(c)
        );
    }

    std::shared_ptr<const CfemExpression> identity()
    {
        auto one = std::make_shared<const ConstantExpression>(1.0);
        return diag(one, one, one);
    }

    std::shared_ptr<const CfemExpression> grad(const std::shared_ptr<const CfemExpression>& e) {
        if (!e) return nullptr;
        return e->grad();
    }

    /**
     * @brief Computes the divergence of a vector or tensor expression.
     * - For a Vector: returns a Scalar (trace of the gradient).
     * - For a Tensor: returns a Vector (row-wise divergence).
     */
    std::shared_ptr<const CfemExpression> div(const std::shared_ptr<const CfemExpression>& e) {
        if (!e) return nullptr;

        // Handle Symbolic Vector
        if (auto v = std::dynamic_pointer_cast<const VectorSymbolicExpression>(e)) {
            auto res = v->getComponent(0)->get_expression().diff(x_rcp) + 
                       v->getComponent(1)->get_expression().diff(y_rcp) + 
                       v->getComponent(2)->get_expression().diff(z_rcp);

            return scalar(res);
        }

        // Handle Symbolic Tensor
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

            return vector(scalar(r0), scalar(r1), scalar(r2));
        }

        CfemInt dim = e->getNumComponents();

        // --- Divergence of a Vector (Dim 3) -> Returns a Scalar ---
        if (dim == 3) {
            // div(v) = dvx/dx + dvy/dy + dvz/dz
            auto g = e->grad(); 
            auto dvx_dx = utils::extractComponent(g, 0); // Grad component (0,0)
            auto dvy_dy = utils::extractComponent(g, 4); // Grad component (1,1)
            auto dvz_dz = utils::extractComponent(g, 8); // Grad component (2,2)
            
            return (dvx_dx + dvy_dy) + dvz_dz;
        }

        // --- Divergence of a Tensor (Dim 9) -> Returns a Vector ---
        if (dim == 9) {
            // Row 0 divergence
            auto row0_div = (utils::extractComponent(e->grad(), 0) + // dT00/dx
                             utils::extractComponent(e->grad(), 4) + // dT01/dy
                             utils::extractComponent(e->grad(), 8)); // dT02/dz
                            
            // Row 1 divergence
            auto row1_div = (utils::extractComponent(e->grad(), 9) +  // Wait: Indexing depends on Grad(Tensor) layout
                             utils::extractComponent(e->grad(), 13) + 
                             utils::extractComponent(e->grad(), 17)); 

            // Alternative safe implementation using row extraction:
            auto row0 = vector(utils::extractComponent(e, 0), utils::extractComponent(e, 1), utils::extractComponent(e, 2));
            auto row1 = vector(utils::extractComponent(e, 3), utils::extractComponent(e, 4), utils::extractComponent(e, 5));
            auto row2 = vector(utils::extractComponent(e, 6), utils::extractComponent(e, 7), utils::extractComponent(e, 8));

            return vector(div(row0), div(row1), div(row2));
        }

        CFEM_ERROR << "div() is not defined for expressions with dimension " << dim;
        return nullptr;
    }

    std::shared_ptr<const CfemExpression> mag_sqr(std::shared_ptr<const CfemExpression> a) {
        return std::make_shared<MagnitudeSqExpression>(std::move(a));
    }

    std::shared_ptr<const CfemExpression> scalar(std::shared_ptr<const CfemExpression> expr)
    {
        if (expr->getNumComponents() != 1 ) CFEM_ERROR << "expr must be a scalar (dimension 1) expression.";
        return expr;
    }

    std::shared_ptr<const CfemExpression> scalar(const SymEngine::Expression& ex) 
    {
        return std::make_shared<ScalarSymbolicExpression>(ex, std::vector{x_rcp, y_rcp, z_rcp, t_rcp});
    }

    std::shared_ptr<const CfemExpression> vector(const SymEngine::Expression& ex, 
                                                 const SymEngine::Expression& ey, 
                                                 const SymEngine::Expression& ez) 
    {
        return std::make_shared<VectorSymbolicExpression>(scalar(ex), scalar(ey), scalar(ez));
    }

    std::shared_ptr<const CfemExpression> variable(const std::string name, CfemInt dimIndex)
    {
        return std::make_shared<CoordinateExpression>(name, dimIndex);
    }

    std::shared_ptr<const CfemExpression> geometric_normal() {
        return std::make_shared<GeometricNormalExpression>();
    }

    std::shared_ptr<const CfemExpression>  grad_s(std::shared_ptr<const CfemExpression> expr, 
                                                  std::shared_ptr<const CfemExpression> normal) 
    {
        auto normale = (normal == nullptr) ? geometric_normal() : normal;
        auto P = identity() - outer(normale, normale); // Projector
        
        if (expr->getNumComponents() == 1) return P * grad(expr);
        else return grad(expr) * P;
    }

    SymEngine::Expression sin(const SymEngine::Expression &e) { return SymEngine::sin(e); }
    SymEngine::Expression cos(const SymEngine::Expression &e) { return SymEngine::cos(e); }
    SymEngine::Expression tan(const SymEngine::Expression &e) { return SymEngine::tan(e); }
    SymEngine::Expression cot(const SymEngine::Expression &e) { return SymEngine::cot(e); }
    SymEngine::Expression sec(const SymEngine::Expression &e) { return SymEngine::sec(e); }
    SymEngine::Expression csc(const SymEngine::Expression &e) { return SymEngine::csc(e); }

    SymEngine::Expression asin(const SymEngine::Expression &e) { return SymEngine::asin(e); }
    SymEngine::Expression acos(const SymEngine::Expression &e) { return SymEngine::acos(e); }
    SymEngine::Expression atan(const SymEngine::Expression &e) { return SymEngine::atan(e); }
    SymEngine::Expression acot(const SymEngine::Expression &e) { return SymEngine::acot(e); }
    SymEngine::Expression asec(const SymEngine::Expression &e) { return SymEngine::asec(e); }
    SymEngine::Expression acsc(const SymEngine::Expression &e) { return SymEngine::acsc(e); }

    SymEngine::Expression sinh(const SymEngine::Expression &e) { return SymEngine::sinh(e); }
    SymEngine::Expression cosh(const SymEngine::Expression &e) { return SymEngine::cosh(e); }
    SymEngine::Expression tanh(const SymEngine::Expression &e) { return SymEngine::tanh(e); }
    SymEngine::Expression coth(const SymEngine::Expression &e) { return SymEngine::coth(e); }
    SymEngine::Expression asinh(const SymEngine::Expression &e) { return SymEngine::asinh(e); }
    SymEngine::Expression acosh(const SymEngine::Expression &e) { return SymEngine::acosh(e); }
    SymEngine::Expression atanh(const SymEngine::Expression &e) { return SymEngine::atanh(e); }

    SymEngine::Expression exp(const SymEngine::Expression &e) { return SymEngine::exp(e); }
    SymEngine::Expression log(const SymEngine::Expression &e) { return SymEngine::log(e); }
    SymEngine::Expression ln(const SymEngine::Expression &e) { return SymEngine::log(e); }

    SymEngine::Expression sqrt(const SymEngine::Expression &e) { return SymEngine::sqrt(e); }
    SymEngine::Expression cbrt(const SymEngine::Expression &e) { return SymEngine::cbrt(e); }
    SymEngine::Expression pow(const SymEngine::Expression &b, const SymEngine::Expression &e) { return SymEngine::pow(b, e); }

    SymEngine::Expression abs(const SymEngine::Expression &e) { return SymEngine::abs(e); }

    SymEngine::Expression max(const SymEngine::Expression &a, const SymEngine::Expression &b) { return SymEngine::max({a, b}); }
    SymEngine::Expression min(const SymEngine::Expression &a, const SymEngine::Expression &b) { return SymEngine::min({a, b}); }

    SymEngine::Expression gamma(const SymEngine::Expression &e) { return SymEngine::gamma(e); }
    SymEngine::Expression erf(const SymEngine::Expression &e) { return SymEngine::erf(e); }
    SymEngine::Expression erfc(const SymEngine::Expression &e) { return SymEngine::erfc(e); }
    SymEngine::Expression zeta(const SymEngine::Expression &e) { return SymEngine::zeta(e); }
}

namespace cfem::sym::utils {

    inline std::shared_ptr<const CfemExpression> extractComponent(std::shared_ptr<const CfemExpression> expr, CfemInt idx) {
        if (!expr) return nullptr;

        // Direct Access (The "Leaf" nodes or explicit containers)
        if (auto v = std::dynamic_pointer_cast<const VectorExpression>(expr))         return v->getComponent(idx);
        if (auto v = std::dynamic_pointer_cast<const TensorExpression>(expr))         return v->getComponent(idx);
        if (auto v = std::dynamic_pointer_cast<const VectorSymbolicExpression>(expr)) return v->getComponent(idx);
        if (auto v = std::dynamic_pointer_cast<const TensorSymbolicExpression>(expr)) return v->getComponent(idx);

        // Scaled Expression (Linear)
        if (auto v = std::dynamic_pointer_cast<const ScaledExpression>(expr)) {
            return std::make_shared<const ScaledExpression>(v->getCoeff(), extractComponent(v->getBaseShared(), idx));
        }

        // Unary Expressions
        if (auto v = std::dynamic_pointer_cast<const UnaryExpression>(expr)) {
            MathOp op = v->getOp();
            // Handle Transpose: Swap index based on 3x3 flattening
            if (op == MathOp::Transpose) {
                static const CfemInt t_map[9] = {0, 3, 6, 1, 4, 7, 2, 5, 8};
                return extractComponent(v->getInner(), t_map[idx]);
            }
            // Pointwise (Sin, Cos, Exp, etc.): Apply op to the extracted component
            if (v->getNumComponents() > 1) {
                return std::make_shared<const UnaryExpression>(extractComponent(v->getInner(), idx), op);
            }
            return v; // Already a scalar
        }

        // Binary Expressions
        if (auto b = std::dynamic_pointer_cast<const BinaryExpression>(expr)) {
            MathOp op = b->getOp();
            auto lhs = b->getLHS();
            auto rhs = b->getRHS();

            if (op == MathOp::Add || op == MathOp::Sub) {
                return std::make_shared<const BinaryExpression>(extractComponent(lhs, idx), extractComponent(rhs, idx), op);
            }
            
            if (op == MathOp::Mul) {
                if (lhs->getNumComponents() == 1) return lhs * extractComponent(rhs, idx);
                if (rhs->getNumComponents() == 1) return extractComponent(lhs, idx) * rhs;
            }

            if (op == MathOp::Div && rhs->getNumComponents() == 1) {
                return extractComponent(lhs, idx) / rhs;
            }

            // Cross Product Expansion: u x v
            if (op == MathOp::Cross) {
                if (idx == 0) return extractComponent(lhs, 1)*extractComponent(rhs, 2) - extractComponent(lhs, 2)*extractComponent(rhs, 1);
                if (idx == 1) return extractComponent(lhs, 2)*extractComponent(rhs, 0) - extractComponent(lhs, 0)*extractComponent(rhs, 2);
                if (idx == 2) return extractComponent(lhs, 0)*extractComponent(rhs, 1) - extractComponent(lhs, 1)*extractComponent(rhs, 0);
            }
            
            if (op == MathOp::Dot || op == MathOp::Inner) return b; // Result is already scalar
        }

        // Fallback: Lambda Evaluation (The Safety Net)
        auto eval_func = [expr, idx](const Vector3D& pos, CfemReal t) -> CfemReal {
            CfemReal temp[9] = {0.0};
            EvalContext ctx;
            ctx.singleWS.physPos = pos;
            ctx.currentTime = t;
            expr->evaluate(-1, pos, temp, ctx); 
            return temp[idx];
        };

        std::shared_ptr<const CfemExpression> comp_time_deriv = nullptr;
        auto td = expr->timeDerivative();
        if (td && !td->isZero()) comp_time_deriv = extractComponent(td, idx);

        return LambdaExpression<decltype(eval_func)>::fromScalar(std::move(eval_func), nullptr, comp_time_deriv);
    }
    
}