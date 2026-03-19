#pragma once


#include "BinaryExpression.h"
#include "CachedExpression.h"
#include "ConstantExpression.h"
#include "CoordinateExpression.h"
#include "GradExpression.h"
#include "LambdaExpression.h"
#include "MagnitudeSqExpression.h"
#include "UnaryExpression.h"
#include "VariableExpression.h"
#include "symengine/ScalarSymbolicExpression.h"
#include "symengine/SymEngineSymbols.h"
#include "symengine/VectorSymbolicExpression.h"
#include "symengine/TensorSymbolicExpression.h"
#include "VectorExpression.h"
#include "TensorExpression.h"
#include "ZeroExpression.h"
#include "ScaledExpression.h"
#include "AlgebraicSimplifier.h"
#include "MathOperators.h"


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

    // inline bool isIdentical(const std::shared_ptr<const CfemExpression>& a, 
    //                 const std::shared_ptr<const CfemExpression>& b) {
    //     if (!a || !b) return a == b;
    //     if (a.get() == b.get()) return true; // Fastest check
    //     if (a->getDimension() != b->getDimension()) return false;

    //     auto cA = std::dynamic_pointer_cast<const ConstantExpression>(a);
    //     auto cB = std::dynamic_pointer_cast<const ConstantExpression>(b);
    //     if (cA && cB) {
    //         return cA->getValues() == cB->getValues();
    //     }

    //     // 1. Check for UnaryExpression equality (e.g., sin(u) == sin(u))
    //     auto uA = std::dynamic_pointer_cast<const UnaryExpression>(a);
    //     auto uB = std::dynamic_pointer_cast<const UnaryExpression>(b);
    //     if (uA && uB) {
    //         return (uA->getOp() == uB->getOp() && isIdentical(uA->getInner(), uB->getInner()));
    //     }

    //     // 2. Check for ScalarSymbolicExpression (SymEngine wrapper)
    //     auto sA = std::dynamic_pointer_cast<const ScalarSymbolicExpression>(a);
    //     auto sB = std::dynamic_pointer_cast<const ScalarSymbolicExpression>(b);
    //     if (sA && sB) {
    //         return sA->get_expression() == sB->get_expression(); 
    //     }

    //     auto vA = std::dynamic_pointer_cast<const VectorSymbolicExpression>(a);
    //     auto vB = std::dynamic_pointer_cast<const VectorSymbolicExpression>(b);
    //     if (vA && vB) {
    //         for (int i = 0; i < 3; ++i) {
    //             if (!isIdentical(vA->getComponent(i), vB->getComponent(i))) return false;
    //         }
    //         return true;
    //     }

    //     // 5. Tensor Check (assuming 3x3 or Voigt notation)
    //     auto tA = std::dynamic_pointer_cast<const TensorSymbolicExpression>(a);
    //     auto tB = std::dynamic_pointer_cast<const TensorSymbolicExpression>(b);
    //     if (tA && tB) {
    //         for (size_t i = 0; i < 9; ++i) {
    //             if (!isIdentical(tA->getComponent(i), tB->getComponent(i))) return false;
    //         }
    //         return true;
    //     }

    //     return false; // Default: assume different if not explicitly matched
    // }


    inline std::shared_ptr<const CfemExpression> operator-(std::shared_ptr<const CfemExpression> u) 
    {
        if (u->isZero()) return u;

        // 2. Constant folding check
        if (auto cU = std::dynamic_pointer_cast<const ConstantExpression>(u)) {
            const auto& v = cU->getValues();
            std::vector<CfemReal> res(v.size());
            for(size_t i = 0; i < v.size(); ++i) res[i] = -v[i];
            return std::make_shared<ConstantExpression>(std::move(res));
        }

        // 3. Smart Scaling: - (k * u) -> (-k) * u
        if (auto uScaled = std::dynamic_pointer_cast<const ScaledExpression>(u)) {
            return std::make_shared<ScaledExpression>(-uScaled->getCoeff(), uScaled->getBaseShared());
        }

        // 4. Default: Treat -u as (-1.0) * u
        return std::make_shared<ScaledExpression>(-1.0, u);
    }

    inline std::shared_ptr<const CfemExpression> operator-(std::shared_ptr<const CfemExpression> lhs, 
                                                       std::shared_ptr<const CfemExpression> rhs) 
    {
        return lhs + (-rhs); 
    }



    // inline std::shared_ptr<const CfemExpression> operator+(std::shared_ptr<const CfemExpression> lhs, 
    //                                                    std::shared_ptr<const CfemExpression> rhs) 
    // {
    //     if (lhs->getDimension() != rhs->getDimension()) {
    //         std::cout << "Dim left: " << lhs->getDimension()  << " Dim Right: " << rhs->getDimension() << std::endl;
    //         CFEM_ERROR << "Dimension mismatch in operator << ";
    //     }

    //     if (rhs->isZero()) return lhs;
    //     if (lhs->isZero()) return rhs;

    //     if (lhs.get() == rhs.get()) {
    //         return std::make_shared<const ScaledExpression>(2.0, lhs); 
    //     }

    //     // 2. Unary Structural Check
    //     auto lU = std::dynamic_pointer_cast<const UnaryExpression>(lhs);
    //     auto rU = std::dynamic_pointer_cast<const UnaryExpression>(rhs);

    //     if (lU && rU) {
    //         // If it's the same math operation (e.g. SIN) acting on the exact same object
    //         if (lU->getOp() == rU->getOp() && lU->getInner().get() == rU->getInner().get()) {
    //             return std::make_shared<const ScaledExpression>(2.0, lhs);
    //         }
    //     }

    //     // --- Scale Awareness ---
    //     auto lS = std::dynamic_pointer_cast<const ScaledExpression>(lhs);
    //     auto rS = std::dynamic_pointer_cast<const ScaledExpression>(rhs);

    //     if (lS && rS && lS->getBasePtr() == rS->getBasePtr()) {
    //         return std::make_shared<const ScaledExpression>(lS->getCoeff() + rS->getCoeff(), lS->getBaseShared());
    //     }
    //     if (lS && lS->getBasePtr() == rhs.get()) {
    //         auto coeff = lS->getCoeff()+ 1.0;
    //         if (NumericPolicy::isZero(coeff)) 
    //             return std::make_shared<const ZeroExpression>(lS->getDimension());
    //         return std::make_shared<const ScaledExpression>(coeff, rhs);
    //     }
    //     if (rS && rS->getBasePtr() == lhs.get()) {
    //         auto coeff = rS->getCoeff()+ 1.0;
    //         if (NumericPolicy::isZero(coeff)) 
    //             return std::make_shared<const ZeroExpression>(rS->getDimension());
    //         return std::make_shared<const ScaledExpression>(coeff, lhs);
    //     }

    //     // --- RECURSIVE PEEK LOGIC ---
    //     auto bl = std::dynamic_pointer_cast<const BinaryExpression>(lhs);
    //     auto br = std::dynamic_pointer_cast<const BinaryExpression>(rhs);

    //     bool blIsSum = bl && (bl->getOp() == MathOp::Add || bl->getOp() == MathOp::Sub);
    //     bool brIsSum = br && (br->getOp() == MathOp::Add || br->getOp() == MathOp::Sub);

    //     auto getEffectiveRHS = [](const std::shared_ptr<const BinaryExpression>& b) -> std::shared_ptr<const CfemExpression> {
    //         if (b->getOp() == MathOp::Add) return b->getRHS();
    //         if (b->getOp() == MathOp::Sub) return std::make_shared<const ScaledExpression>(-1.0, b->getRHS());
    //         return nullptr;
    //     };

    //     auto isSimplified = [](const std::shared_ptr<const CfemExpression>& e) {
    //         auto bin = std::dynamic_pointer_cast<const BinaryExpression>(e);
    //         return !bin || (bin->getOp() != MathOp::Add && bin->getOp() != MathOp::Sub);
    //     };

    //     // PRIORITY 1: Try merging entire RHS into LHS's branches
    //     if (blIsSum) {
    //         auto L1 = bl->getLHS();
    //         auto L2 = getEffectiveRHS(bl);

    //         auto tryRight = L2 + rhs;
    //         if (isSimplified(tryRight)) return L1 + tryRight;

    //         auto tryLeft = L1 + rhs;
    //         if (isSimplified(tryLeft)) return tryLeft + L2;
    //     } 
        
    //     // PRIORITY 2: Try merging entire LHS into RHS's branches
    //     if (brIsSum) {
    //         auto R1 = br->getLHS();
    //         auto R2 = getEffectiveRHS(br);

    //         auto tryLeft = lhs + R1;
    //         if (isSimplified(tryLeft)) return tryLeft + R2;

    //         auto tryRight = lhs + R2;
    //         if (isSimplified(tryRight)) return R1 + tryRight;
    //     }

    //     // PRIORITY 3: Only if whole-tree merges failed, try dismantling both
    //     if (blIsSum && brIsSum) {
    //         auto L1 = bl->getLHS();
    //         auto L2 = getEffectiveRHS(bl);
    //         auto R1 = br->getLHS();
    //         auto R2 = getEffectiveRHS(br);

    //         auto tryLL = L1 + R1;
    //         if (isSimplified(tryLL)) return tryLL + (L2 + R2);

    //         auto tryLR = L1 + R2;
    //         if (isSimplified(tryLR)) return tryLR + (L2 + R1);

    //         auto tryRL = L2 + R1;
    //         if (isSimplified(tryRL)) return tryRL + (L1 + R2);

    //         auto tryRR = L2 + R2;
    //         if (isSimplified(tryRR)) return tryRR + (L1 + R1);
    //     }

    //     // Fallback: Just build a standard Binary Add node
    //     return foldBinary(lhs, rhs, MathOp::Add, lhs->getDimension(), 
    //         [](const auto& vL, const auto& vR, auto& res) {
    //             for(size_t i = 0; i < vL.size(); ++i) res[i] = vL[i] + vR[i];
    //         });
    // }

    inline std::shared_ptr<const CfemExpression>
    operator+(std::shared_ptr<const CfemExpression> lhs, std::shared_ptr<const CfemExpression> rhs)
    {
        if (lhs->getDimension() != rhs->getDimension())
            CFEM_ERROR << "Dimension mismatch in operator << ";

        if (lhs->isZero()) return rhs;
        if (rhs->isZero()) return lhs;

        return simplifyAdd(lhs, rhs);
    }

    inline std::shared_ptr<const CfemExpression> operator*(std::shared_ptr<const CfemExpression> lhs, 
                                                        std::shared_ptr<const CfemExpression> rhs) 
    {
        CfemInt dimL = lhs->getDimension();
        CfemInt dimR = rhs->getDimension();

        if (dimL != dimR && dimL != 1 && dimR != 1) {
            CFEM_ERROR << "Invalid operator*";
        }

        CfemInt outDim = std::max(dimL, dimR);

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

    
    inline std::shared_ptr<const CfemExpression> operator/(std::shared_ptr<const CfemExpression> lhs, 
                                                           std::shared_ptr<const CfemExpression> rhs) 
    {
        CfemInt dimL = lhs->getDimension();
        CfemInt dimR = rhs->getDimension();

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

    inline std::shared_ptr<const CfemExpression> identityO2TensorExpr(){
        return std::make_shared<const ConstantExpression>(IdentityTensor());
    }


    inline std::shared_ptr<const CfemExpression> transpose(std::shared_ptr<const CfemExpression> u) {
        if (u->getDimension() != 9) CFEM_ERROR << "`transpose()` expects a second order tensor argument."; 
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


    inline std::shared_ptr<const CfemExpression> trace(std::shared_ptr<const CfemExpression> u){
        if (u->getDimension() != 9) CFEM_ERROR << "`trace()` expects a second order tensor argument."; 
        if (u->isZero()) return std::make_shared<ZeroExpression>(0.0);

        if (u->isConstant()) {
            auto cU = std::static_pointer_cast<const ConstantExpression>(u);
            const auto& v = cU->getValues();
            return std::make_shared<ConstantExpression>(v[0] + v[4] + v[8]);
        }

        return std::make_shared<UnaryExpression>(u, MathOp::Trace);
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
            CFEM_ERROR << "inner() product (double contraction) requires two dim 9 tensors.";
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

    inline std::shared_ptr<const CfemExpression> atan2(std::shared_ptr<const CfemExpression> y, 
                                                       std::shared_ptr<const CfemExpression> x)
    {
        CfemInt dimY = y->getDimension();
        CfemInt dimX = x->getDimension();

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

    inline std::shared_ptr<const CfemExpression> pow(std::shared_ptr<const CfemExpression> base, 
                                                     std::shared_ptr<const CfemExpression> exp)
    {
        CfemInt dimBase = base->getDimension();
        CfemInt dimExp  = exp->getDimension();

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

    inline std::shared_ptr<const CfemExpression> scalar(
                                            const std::function<CfemReal(const Vector3D &, CfemReal)> &f,
                                            std::shared_ptr<const CfemExpression> gradient,
                                            std::shared_ptr<const CfemExpression> timeDeriv)
    {
        return LambdaExpression::fromScalar(f, gradient, timeDeriv);
    }

    inline std::shared_ptr<const CfemExpression> scalar(std::shared_ptr<const CfemExpression> expr)
    {
        if (expr->getDimension() != 1 ) CFEM_ERROR << "expr must be a scalar (dimension 1) expression.";
        return expr;
    }

    inline std::shared_ptr<const CfemExpression> constant(CfemReal c)
    {
        return std::make_shared<const ConstantExpression>(c);
    }

    inline std::shared_ptr<const CfemExpression> vector(const std::function<Vector3D(const Vector3D &)> &f,
                                                        std::shared_ptr<const CfemExpression> gradient)
    {
        return LambdaExpression::fromVector(f, gradient);
    }

    inline std::shared_ptr<const CfemExpression> vector(
                                            const std::function<Vector3D(const Vector3D &, CfemReal)> &f,
                                            std::shared_ptr<const CfemExpression> gradient,
                                            std::shared_ptr<const CfemExpression> timeDeriv)
    {
        return LambdaExpression::fromVector(f, gradient, timeDeriv);
    }

    inline std::shared_ptr<const CfemExpression> vector(std::shared_ptr<const CfemExpression> expr)
    {
        if (expr->getDimension() != 3 ) CFEM_ERROR << "expr must be a vector (dimesion 3) expression.";
        return expr;
    }

    inline std::shared_ptr<const CfemExpression> vector(std::shared_ptr<const CfemExpression> expr_x, 
                                                 std::shared_ptr<const CfemExpression> expr_y,
                                                 std::shared_ptr<const CfemExpression> expr_z)
    {
        if (expr_x->getDimension() != 1 || expr_y->getDimension() != 1 || expr_z->getDimension() != 1 )
            CFEM_ERROR << "expr_x, expr_y and expr_z must be scalar (dimesion 1) expressions.";

        return std::make_shared<const VectorExpression>(expr_x, expr_y, expr_z);
    }

    inline std::shared_ptr<const CfemExpression> tensor(
                                            const std::function<TensorO2(const Vector3D &)> &f,
                                            std::shared_ptr<const CfemExpression> gradient )
    {
        return LambdaExpression::fromTensor(f, gradient);
    }

    inline std::shared_ptr<const CfemExpression> tensor(
                                            const std::function<TensorO2(const Vector3D &, CfemReal)> &f,
                                            std::shared_ptr<const CfemExpression> gradient,
                                            std::shared_ptr<const CfemExpression> timeDeriv)
    {
        return LambdaExpression::fromTensor(f, gradient, timeDeriv);
    }

    inline std::shared_ptr<const CfemExpression> tensor(std::shared_ptr<const CfemExpression> expr)
    {
        if (expr->getDimension() != 9 ) CFEM_ERROR << "expr must be a scalar expression.";
        return expr;
    }

    inline std::shared_ptr<const CfemExpression> tensor(std::shared_ptr<const CfemExpression> expr_xx,
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

    inline std::shared_ptr<const CfemExpression> symmetricTensor(std::shared_ptr<const CfemExpression> xx, 
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
    inline std::shared_ptr<const CfemExpression> diag(std::shared_ptr<const CfemExpression> a,
                                                      std::shared_ptr<const CfemExpression> b,
                                                      std::shared_ptr<const CfemExpression> c)
    {
        return std::make_shared<const TensorExpression>(
            a,    zeros::Zero1D, zeros::Zero1D,
            zeros::Zero1D, b,    zeros::Zero1D,
            zeros::Zero1D, zeros::Zero1D, c
        );
    }

    inline std::shared_ptr<const CfemExpression> diag(CfemReal a,
                                                      CfemReal b,
                                                      CfemReal c)
    {
        return std::make_shared<const TensorExpression>(
            constant(a),    zeros::Zero1D, zeros::Zero1D,
            zeros::Zero1D, constant(b),    zeros::Zero1D,
            zeros::Zero1D, zeros::Zero1D, constant(c)
        );
    }

    /**
     * @brief Standard Identity matrix factory
     */
    inline std::shared_ptr<const CfemExpression> identity()
    {
        auto one = std::make_shared<const ConstantExpression>(1.0);
        return diag(one, one, one);
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
            CFEM_ERROR << "TODO: Gradient of a second order tensor currentely not supported. Please implement the class TensorO3SymbolicExpression.";

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

    // inline std::shared_ptr<const CfemExpression> div(const std::shared_ptr<const CfemExpression>& e) {
    //     if (auto v = std::dynamic_pointer_cast<const VectorSymbolicExpression>(e)) {
    //         auto scalar_comp = std::static_pointer_cast<const ScalarSymbolicExpression>(v->getComponent(0));
    //         auto res = v->getComponent(0)->get_expression().diff(x_rcp) + 
    //                    v->getComponent(1)->get_expression().diff(y_rcp) + 
    //                    v->getComponent(2)->get_expression().diff(z_rcp);
    //         return scalar(res);
    //     }

    //     if (auto t = std::dynamic_pointer_cast<const TensorSymbolicExpression>(e)) {
    //         auto r0 = t->getComponent(0)->get_expression().diff(x_rcp) + 
    //                   t->getComponent(1)->get_expression().diff(y_rcp) + 
    //                   t->getComponent(2)->get_expression().diff(z_rcp);
    //         auto r1 = t->getComponent(3)->get_expression().diff(x_rcp) + 
    //                   t->getComponent(4)->get_expression().diff(y_rcp) + 
    //                   t->getComponent(5)->get_expression().diff(z_rcp);
    //         auto r2 = t->getComponent(6)->get_expression().diff(x_rcp) + 
    //                   t->getComponent(7)->get_expression().diff(y_rcp) + 
    //                   t->getComponent(8)->get_expression().diff(z_rcp);

            
    //         // TODO : Etendre aux lambda
    //         return std::make_shared<VectorSymbolicExpression>(scalar(r0), scalar(r1), scalar(r2));
    //     }

    //     CFEM_ERROR << "div() only supports Symbolic Vector/Tensor types currently.";
    //     return nullptr;
    // }

    inline std::shared_ptr<const CfemExpression> div(const std::shared_ptr<const CfemExpression>& e) {
        // 1. Handle Symbolic Vector
        if (auto v = std::dynamic_pointer_cast<const VectorSymbolicExpression>(e)) {
            auto res = v->getComponent(0)->get_expression().diff(x_rcp) + 
                    v->getComponent(1)->get_expression().diff(y_rcp) + 
                    v->getComponent(2)->get_expression().diff(z_rcp);
            return scalar(res);
        }

        // 2. Handle Symbolic Tensor
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

            return std::make_shared<VectorSymbolicExpression>(scalar(r0), scalar(r1), scalar(r2));
        }

        // 3. Handle General Expressions (Lambdas, Variables, or Custom Containers)
        CfemInt dim = e->getDimension();
        
        // 2. Recursive Lambda to extract the i-th scalar component safely from ANY vector expression
        std::function<std::shared_ptr<const CfemExpression>(std::shared_ptr<const CfemExpression>, int)> extractComponent = 
            [&extractComponent](std::shared_ptr<const CfemExpression> expr, int idx) -> std::shared_ptr<const CfemExpression> 
        {
            if (!expr) return nullptr;

            // --- Transparent Container Types ---
            if (auto v = std::dynamic_pointer_cast<const VectorSymbolicExpression>(expr)) {
                return v->getComponent(idx);
            }
            else if (auto v = std::dynamic_pointer_cast<const VectorExpression>(expr)) {
                return v->getComponent(idx);
            }
            
            // --- Unary and Scaled Modifiers ---
            else if (auto v = std::dynamic_pointer_cast<const ScaledExpression>(expr)) {
                return std::make_shared<ScaledExpression>(v->getCoeff(), extractComponent(v->getBaseShared(), idx));
            }
            else if (auto v = std::dynamic_pointer_cast<const UnaryExpression>(expr)) {
                return std::make_shared<UnaryExpression>(extractComponent(v->getInner(), idx), v->getOp());
            }

            // --- Binary Operations (The Missing Piece) ---
            else if (auto b = std::dynamic_pointer_cast<const BinaryExpression>(expr)) {
                MathOp op = b->getOp();
                
                // Vector +/- Vector
                if (op == MathOp::Add || op == MathOp::Sub) {
                    return std::make_shared<BinaryExpression>(
                        extractComponent(b->getLHS(), idx), 
                        extractComponent(b->getRHS(), idx), op);
                }
                // Scalar * Vector  OR  Vector * Scalar
                else if (op == MathOp::Mul) {
                    if (b->getLHS()->getDimension() == 1) {
                        return std::make_shared<BinaryExpression>(b->getLHS(), extractComponent(b->getRHS(), idx), op);
                    } else if (b->getRHS()->getDimension() == 1) {
                        return std::make_shared<BinaryExpression>(extractComponent(b->getLHS(), idx), b->getRHS(), op);
                    }
                }
                // Vector / Scalar
                else if (op == MathOp::Div) {
                    if (b->getRHS()->getDimension() == 1) {
                        return std::make_shared<BinaryExpression>(extractComponent(b->getLHS(), idx), b->getRHS(), op);
                    }
                }
                // Cross Product: A x B
                else if (op == MathOp::Cross) {
                    // C_idx = A_i1 * B_i2 - A_i2 * B_i1
                    int i1 = (idx + 1) % 3;
                    int i2 = (idx + 2) % 3;
                    auto term1 = std::make_shared<BinaryExpression>(
                        extractComponent(b->getLHS(), i1), extractComponent(b->getRHS(), i2), MathOp::Mul);
                    auto term2 = std::make_shared<BinaryExpression>(
                        extractComponent(b->getLHS(), i2), extractComponent(b->getRHS(), i1), MathOp::Mul);
                    return std::make_shared<BinaryExpression>(term1, term2, MathOp::Sub);
                }
                // If it's a Pow, Max, Min, etc. applied to vectors, it will fall through to the lambda wrapper below.
            }

            // --- Opaque Black Box Fallback (Lambdas, Variables, Unhandled Ops) ---
            // We use the static factory 'fromScalar' from your LambdaExpression class!
            auto eval_func = [expr, idx](const Vector3D& pos, CfemReal t) -> CfemReal {
                CfemReal temp[3] = {0.0, 0.0, 0.0};
                EvalContext ctx;
                ctx.physPos = pos;
                ctx.currentTime = t;
                expr->evaluate(-1, pos, temp, ctx); 
                return temp[idx];
            };

            std::shared_ptr<const CfemExpression> comp_time_deriv = nullptr;
            try {
                if (auto td = expr->timeDerivative()) {
                    comp_time_deriv = extractComponent(td, idx);
                }
            } catch (...) { /* Safe to ignore if no time derivative exists */ }

            return LambdaExpression::fromScalar(eval_func, nullptr, comp_time_deriv);
        };

        // Divergence of a Vector (Dim 3) -> Returns a Scalar
        if (dim == 3) {
            // div(v) = trace(grad(v))
            auto g = e->grad(); // This returns a Tensor-like expression
            // We use the same component extraction logic we built for grad()
            auto dvx_dx = extractComponent(g, 0);
            auto dvy_dy = extractComponent(g, 4);
            auto dvz_dz = extractComponent(g, 8);
            
            // Return a Binary Add tree: (du/dx + dv/dy) + dw/dz
            return (dvx_dx + dvy_dy) + dvz_dz;
        }

        // Divergence of a Tensor (Dim 9) -> Returns a Vector
        // if (dim == 9) {
        //     auto g = e->grad(); // Note: Gradient of a 2nd order tensor is a 3rd order tensor
        //     // If your system doesn't support 3rd order tensors yet, we use the row-wise logic:
        //     // We assume 'e' can be split into rows.
            
        //     auto row0_div = (extractComponent(e, 0)->grad()->getComponent(0) + 
        //                     extractComponent(e, 1)->grad()->getComponent(1) + 
        //                     extractComponent(e, 2)->grad()->getComponent(2));
                            
        //     auto row1_div = (extractComponent(e, 3)->grad()->getComponent(0) + 
        //                     extractComponent(e, 4)->grad()->getComponent(1) + 
        //                     extractComponent(e, 5)->grad()->getComponent(2));
                            
        //     auto row2_div = (extractComponent(e, 6)->grad()->getComponent(0) + 
        //                     extractComponent(e, 7)->grad()->getComponent(1) + 
        //                     extractComponent(e, 8)->grad()->getComponent(2));

        //     return vector(row0_div, row1_div, row2_div);
        // }

        CFEM_ERROR << "div() is not defined for expressions with dimension " << dim;
        return nullptr;
    }

    inline std::shared_ptr<const CfemExpression> mag_sqr(std::shared_ptr<const CfemExpression> a) {
        return std::make_shared<MagnitudeSqExpression>(std::move(a));
    }

    inline std::shared_ptr<const CfemExpression> scalar(const SymEngine::Expression& ex) {
        return std::make_shared<ScalarSymbolicExpression>(ex, std::vector{x_rcp, y_rcp, z_rcp, t_rcp});
    }

    inline std::shared_ptr<const CfemExpression> vector(const SymEngine::Expression& ex, 
                                                        const SymEngine::Expression& ey, 
                                                        const SymEngine::Expression& ez) {
        return std::make_shared<VectorSymbolicExpression>(scalar(ex), scalar(ey), scalar(ez));
    }

    inline std::shared_ptr<const CfemExpression> variable(const std::string name, CfemInt dimIndex){
        return std::make_shared<CoordinateExpression>(name, dimIndex);
    }



    // SymEngine:
    /**
     * @name Common math functions
     *  @{
     */
    inline SymEngine::Expression sin(const SymEngine::Expression &e) { return SymEngine::sin(e); }
    inline SymEngine::Expression cos(const SymEngine::Expression &e) { return SymEngine::cos(e); }
    inline SymEngine::Expression tan(const SymEngine::Expression &e) { return SymEngine::tan(e); }
    inline SymEngine::Expression cot(const SymEngine::Expression &e) { return SymEngine::cot(e); }
    inline SymEngine::Expression sec(const SymEngine::Expression &e) { return SymEngine::sec(e); }
    inline SymEngine::Expression csc(const SymEngine::Expression &e) { return SymEngine::csc(e); }

    inline SymEngine::Expression asin(const SymEngine::Expression &e) { return SymEngine::asin(e); }
    inline SymEngine::Expression acos(const SymEngine::Expression &e) { return SymEngine::acos(e); }
    inline SymEngine::Expression atan(const SymEngine::Expression &e) { return SymEngine::atan(e); }
    inline SymEngine::Expression acot(const SymEngine::Expression &e) { return SymEngine::acot(e); }
    inline SymEngine::Expression asec(const SymEngine::Expression &e) { return SymEngine::asec(e); }
    inline SymEngine::Expression acsc(const SymEngine::Expression &e) { return SymEngine::acsc(e); }

    inline SymEngine::Expression sinh(const SymEngine::Expression &e) { return SymEngine::sinh(e); }
    inline SymEngine::Expression cosh(const SymEngine::Expression &e) { return SymEngine::cosh(e); }
    inline SymEngine::Expression tanh(const SymEngine::Expression &e) { return SymEngine::tanh(e); }
    inline SymEngine::Expression coth(const SymEngine::Expression &e) { return SymEngine::coth(e); }
    inline SymEngine::Expression asinh(const SymEngine::Expression &e) { return SymEngine::asinh(e); }
    inline SymEngine::Expression acosh(const SymEngine::Expression &e) { return SymEngine::acosh(e); }
    inline SymEngine::Expression atanh(const SymEngine::Expression &e) { return SymEngine::atanh(e); }

    inline SymEngine::Expression exp(const SymEngine::Expression &e) { return SymEngine::exp(e); }
    inline SymEngine::Expression log(const SymEngine::Expression &e) { return SymEngine::log(e); }
    inline SymEngine::Expression ln(const SymEngine::Expression &e) { return SymEngine::log(e); }

    inline SymEngine::Expression sqrt(const SymEngine::Expression &e) { return SymEngine::sqrt(e); }
    inline SymEngine::Expression cbrt(const SymEngine::Expression &e) { return SymEngine::cbrt(e); }
    inline SymEngine::Expression pow(const SymEngine::Expression &b, const SymEngine::Expression &e) { return SymEngine::pow(b, e); }

    inline SymEngine::Expression abs(const SymEngine::Expression &e) { return SymEngine::abs(e); }

    inline SymEngine::Expression max(const SymEngine::Expression &a, const SymEngine::Expression &b) { return SymEngine::max({a, b}); }
    inline SymEngine::Expression min(const SymEngine::Expression &a, const SymEngine::Expression &b) { return SymEngine::min({a, b}); }

    inline SymEngine::Expression gamma(const SymEngine::Expression &e) { return SymEngine::gamma(e); }
    inline SymEngine::Expression erf(const SymEngine::Expression &e) { return SymEngine::erf(e); }
    inline SymEngine::Expression erfc(const SymEngine::Expression &e) { return SymEngine::erfc(e); }
    inline SymEngine::Expression zeta(const SymEngine::Expression &e) { return SymEngine::zeta(e); }
    /** @} */

}