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

#include "ContractException.h"

namespace cfem
{

// ============================================================================
// ContractException
// ============================================================================

ContractException::ContractException(const char* file,
                                     CfemInt line,
                                     const char* expr,
                                     const char* kind,
                                     const char* msg)
: std::logic_error("Contract violation")
{
    std::ostringstream os;

    os << "\n========================================\n"
       << "CFEM CONTRACT VIOLATION\n"
       << "----------------------------------------\n"
       << "Type         : " << kind << "\n"
       << "Location     : " << file << ":" << line << "\n"
       << "Expression   : " << expr << "\n";

    if (msg) os << "Info         : "  << msg << "\n";
    os << "========================================\n";

    m_report = os.str();
}

const char* ContractException::what() const noexcept {
    return m_report.c_str();
}

// ============================================================================
// Derived Exceptions
// ============================================================================

AssertionException::AssertionException(const char*  file,
                                       CfemInt line,
                                       const char*  expr,
                                       const char* msg)
    : ContractException(file, line, expr,"Assertion", msg)
{}


PreconditionException::PreconditionException(const char* file,
                                             CfemInt line,
                                             const char* expr,
                                             const char* msg)
    : ContractException(file, line, expr, "Precondition", msg)
{}


PostconditionException::PostconditionException(const char* file,
                                               CfemInt line,
                                               const char* expr,
                                               const char* msg)
    : ContractException(file, line, expr,"Postcondition", msg)
{}


InvariantException::InvariantException(const char*  file,
                                       CfemInt line,
                                       const char*  expr,
                                       const char* msg)
    : ContractException(file, line, expr, "Invariant", msg)
{}

} // namespace cfem
