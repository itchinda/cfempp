#include "ContractException.h"

namespace cfem
{

// ============================================================================
// ContractException
// ============================================================================

ContractException::ContractException(const char* file,
                                     int line,
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
                                       int line,
                                       const char*  expr,
                                       const char* msg)
    : ContractException(file, line, expr,"Assertion", msg)
{}


PreconditionException::PreconditionException(const char* file,
                                             int line,
                                             const char* expr,
                                             const char* msg)
    : ContractException(file, line, expr, "Precondition", msg)
{}


PostconditionException::PostconditionException(const char* file,
                                               int line,
                                               const char* expr,
                                               const char* msg)
    : ContractException(file, line, expr,"Postcondition", msg)
{}


InvariantException::InvariantException(const char*  file,
                                       int line,
                                       const char*  expr,
                                       const char* msg)
    : ContractException(file, line, expr, "Invariant", msg)
{}

} // namespace cfem
