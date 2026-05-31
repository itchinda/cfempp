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

#include <memory>
#include <vector>

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