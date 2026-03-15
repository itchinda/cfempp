#pragma once

#include <vector>
#include <memory>
#include "CfemExpression.h"

namespace cfem::sym
{

    struct Term
    {
        double coeff;
        std::shared_ptr<const CfemExpression> base;
    };

    void collectTerms(
        const std::shared_ptr<const CfemExpression> &expr,
        double coeff,
        std::vector<Term> &terms);

    void removeZeroTerms(std::vector<Term> &terms);

    bool isIdentical(const std::shared_ptr<const CfemExpression> &a,
                     const std::shared_ptr<const CfemExpression> &b);

    std::vector<Term> combineTerms(const std::vector<Term> &terms);

    std::shared_ptr<const CfemExpression> rebuildSum(const std::vector<Term> &terms);

    std::shared_ptr<const CfemExpression> simplifyAdd(std::shared_ptr<const CfemExpression> lhs,
                                                      std::shared_ptr<const CfemExpression> rhs);

}