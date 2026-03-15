#include "AlgebraicSimplifier.h"
#include <vector>
#include <memory>
#include "CfemExpression.h"
#include "ScaledExpression.h"
#include "BinaryExpression.h"
#include "ConstantExpression.h"
#include "UnaryExpression.h"
#include "symengine/ScalarSymbolicExpression.h"
#include "symengine/VectorSymbolicExpression.h"
#include "symengine/TensorSymbolicExpression.h"

namespace cfem::sym
{

    void collectTerms(
        const std::shared_ptr<const CfemExpression> &expr,
        double coeff,
        std::vector<Term> &terms)
    {
        if (auto s = std::dynamic_pointer_cast<const ScaledExpression>(expr))
        {
            collectTerms(s->getBaseShared(), coeff * s->getCoeff(), terms);
            return;
        }

        if (auto b = std::dynamic_pointer_cast<const BinaryExpression>(expr))
        {
            if (b->getOp() == MathOp::Add)
            {
                collectTerms(b->getLHS(), coeff, terms);
                collectTerms(b->getRHS(), coeff, terms);
                return;
            }

            if (b->getOp() == MathOp::Sub)
            {
                collectTerms(b->getLHS(), coeff, terms);
                collectTerms(b->getRHS(), -coeff, terms);
                return;
            }
        }

        terms.push_back({coeff, expr});
    }

    void removeZeroTerms(std::vector<Term> &terms)
    {
        terms.erase(
            std::remove_if(
                terms.begin(),
                terms.end(),
                [](const Term &t)
                {
                    return NumericPolicy::isZero(t.coeff);
                }),
            terms.end());
    }

    bool isIdentical(const std::shared_ptr<const CfemExpression> &a,
                     const std::shared_ptr<const CfemExpression> &b)
    {
        if (!a || !b)
            return a == b;
        if (a.get() == b.get())
            return true; // Fastest check
        if (a->getDimension() != b->getDimension())
            return false;

        auto cA = std::dynamic_pointer_cast<const ConstantExpression>(a);
        auto cB = std::dynamic_pointer_cast<const ConstantExpression>(b);
        if (cA && cB)
        {
            return cA->getValues() == cB->getValues();
        }

        // 1. Check for UnaryExpression equality (e.g., sin(u) == sin(u))
        auto uA = std::dynamic_pointer_cast<const UnaryExpression>(a);
        auto uB = std::dynamic_pointer_cast<const UnaryExpression>(b);
        if (uA && uB)
        {
            return (uA->getOp() == uB->getOp() && isIdentical(uA->getInner(), uB->getInner()));
        }

        // 2. Check for ScalarSymbolicExpression (SymEngine wrapper)
        auto sA = std::dynamic_pointer_cast<const ScalarSymbolicExpression>(a);
        auto sB = std::dynamic_pointer_cast<const ScalarSymbolicExpression>(b);
        if (sA && sB)
        {
            return sA->get_expression() == sB->get_expression();
        }

        auto vA = std::dynamic_pointer_cast<const VectorSymbolicExpression>(a);
        auto vB = std::dynamic_pointer_cast<const VectorSymbolicExpression>(b);
        if (vA && vB)
        {
            for (int i = 0; i < 3; ++i)
            {
                if (!isIdentical(vA->getComponent(i), vB->getComponent(i)))
                    return false;
            }
            return true;
        }

        // 5. Tensor Check (assuming 3x3 or Voigt notation)
        auto tA = std::dynamic_pointer_cast<const TensorSymbolicExpression>(a);
        auto tB = std::dynamic_pointer_cast<const TensorSymbolicExpression>(b);
        if (tA && tB)
        {
            for (size_t i = 0; i < 9; ++i)
            {
                if (!isIdentical(tA->getComponent(i), tB->getComponent(i)))
                    return false;
            }
            return true;
        }

        return false; // Default: assume different if not explicitly matched
    }

    std::vector<Term> combineTerms(const std::vector<Term> &terms)
    {
        std::vector<Term> result;

        for (auto &t : terms)
        {
            bool merged = false;

            for (auto &r : result)
            {
                if (isIdentical(r.base, t.base))
                {
                    r.coeff += t.coeff;
                    merged = true;
                    break;
                }
            }

            if (!merged)
                result.push_back(t);
        }

        return result;
    }

    std::shared_ptr<const CfemExpression> rebuildSum(
        const std::vector<Term> &terms)
    {
        if (terms.empty())
            return std::make_shared<ZeroExpression>(1);

        std::shared_ptr<const CfemExpression> result;

        for (auto &t : terms)
        {
            std::shared_ptr<const CfemExpression> term;

            if (NumericPolicy::areEqual(t.coeff, 1.0))
                term = t.base;
            else
                term = std::make_shared<ScaledExpression>(t.coeff, t.base);

            if (!result)
                result = term;
            else
                result = std::make_shared<BinaryExpression>(result, term, MathOp::Add);
        }

        return result;
    }

    std::shared_ptr<const CfemExpression> simplifyAdd( std::shared_ptr<const CfemExpression> lhs,
                                                       std::shared_ptr<const CfemExpression> rhs)
    {
        std::vector<Term> terms;

        collectTerms(lhs, 1.0, terms);
        collectTerms(rhs, 1.0, terms);

        auto combined = combineTerms(terms);

        removeZeroTerms(combined);

        return rebuildSum(combined);
    }

}