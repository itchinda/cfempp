#ifndef CFEM_CONTRACT_EXCEPTION_H
#define CFEM_CONTRACT_EXCEPTION_H

#include <stdexcept>
#include <string>
#include <sstream>

namespace cfem
{

/**
 * @file ContractException.h
 * @brief Contract-based exception hierarchy for error handling.
 *
 * This file defines a hierarchy of exception classes used to implement
 * Design by Contract (DbC) in the CFEM framework.
 *
 * The purpose of this system is to:
 *  - Enforce preconditions, postconditions, and invariants
 *  - Detect logical programming errors
 *  - Provide detailed diagnostic information
 *  - Improve robustness and maintainability
 *
 * In release mode (NDEBUG defined), all contract checks are disabled.
 *
 * This design is inspired by:
 *   Horstmann, "Mastering Object-Oriented Design in C++"
 *
 * @author Ismael Tchinda Ngueyong
 * @version 2.0
 */


/**
 * @class ContractException
 * @brief Base class for all contract violations.
 *
 * This class stores detailed information about a violated contract,
 * including:
 *  - Source file name
 *  - Line number
 *  - Failed logical expression
 *  - Descriptive error message
 *
 * It derives from std::logic_error since contract violations represent
 * programming errors rather than runtime failures.
 */
class ContractException : public std::logic_error
{
public:

    /**
     * @brief Constructs a ContractException.
     *
     * @param file       Source file where the violation occurred.
     * @param line       Line number of the violation.
     * @param expression Logical expression that failed.
     * @param kind       Type of contract (Precondition, Postcondition, etc.)
     */
    ContractException(const char* file,
                      int line,
                      const char* expr,
                      const char* kind,
                      const char* msg = nullptr);

    /**
     * @brief Return formatted diagnostic message.
     */
    const char* what() const noexcept override;

private:
    std::string m_report;
};



/**
 * @class AssertionException
 * @brief Exception for assertion violations.
 *
 * Thrown when a CFEM_ASSERT condition fails.
 *
 * Assertions are intended to detect internal programming errors
 * that should never occur in correct code.
 */
class AssertionException : public ContractException
{
public:

    /**
     * @brief Constructs an AssertionException.
     *
     * @param file Source file where the violation occurred
     * @param line Line number where the violation occurred
     * @param expr The string representation of the failed logical expression
     * @param msg  Optional custom message providing error clarification
     */
    AssertionException(const char*  file,
                       int line,
                       const char*  expr,
                       const char* msg = nullptr);
};



/**
 * @class PreconditionException
 * @brief Exception for precondition violations.
 *
 * Thrown when a CFEM_PRECONDITION condition fails.
 *
 * Preconditions specify requirements that must be satisfied
 * before calling a function.
 */
class PreconditionException : public ContractException
{
public:

    /**
     * @brief Constructs a PreconditionException.
     * 
     * @param file Source file where the violation occurred
     * @param line Line number where the violation occurred
     * @param expr The string representation of the failed logical expression
     * @param msg  Optional custom message providing error clarification
     */
    PreconditionException(const char*  file,
                          int line,
                          const char*  expr,
                          const char* msg = nullptr);
};



/**
 * @class PostconditionException
 * @brief Exception for postcondition violations.
 *
 * Thrown when a CFEM_POSTCONDITION condition fails.
 *
 * Postconditions specify requirements that must be satisfied
 * after a function completes.
 */
class PostconditionException : public ContractException
{
public:

    /**
     * @brief Constructs a PostconditionException. 
     *
     * @param file Source file where the violation occurred
     * @param line Line number where the violation occurred
     * @param expr The string representation of the failed logical expression
     * @param msg  Optional custom message providing error clarification
     */
    PostconditionException(const char*  file,
                           int line,
                           const char*  expr,
                           const char* msg = nullptr);
};



/**
 * @class InvariantException
 * @brief Exception for invariant violations.
 *
 * Thrown when a CFEM_INVARIANT condition fails.
 *
 * Invariants define conditions that must always hold
 * for an object during its lifetime.
 */
class InvariantException : public ContractException
{
public:

    /**
     * @brief Constructs an AssertionException.
     *
     * @param file Source file where the violation occurred
     * @param line Line number where the violation occurred
     * @param expr The string representation of the failed logical expression
     * @param msg  Optional custom message providing error clarification
     */
    InvariantException(const char*  file,
                       int line,
                       const char*  expr,
                       const char* msg = nullptr);
};

} // namespace cfem



// ============================================================================
// Contract Verification Macros
// ============================================================================


#ifndef NDEBUG   // Debug mode enabled



// Helper: select macro by arg count
#define CFEM_GET_MACRO(_1,_2,NAME,...) NAME


/* ================= ASSERT ================= */

#define CFEM_ASSERT_1(expr) \
    do { \
        if (!(expr)) \
            throw cfem::AssertionException(__FILE__, __LINE__, #expr); \
    } while(0)


#define CFEM_ASSERT_2(expr, msg) \
    do { \
        if (!(expr)) \
            throw cfem::AssertionException(__FILE__, __LINE__, #expr, msg); \
    } while(0)

/**
 * @def CFEM_ASSERT
 * @brief Verifies an internal consistency condition.
 */
#define CFEM_ASSERT(...) \
    CFEM_GET_MACRO(__VA_ARGS__, CFEM_ASSERT_2, CFEM_ASSERT_1)(__VA_ARGS__)


/* ================= PRECONDITION ================= */

#define CFEM_PRE_1(expr) \
    do { \
        if (!(expr)) \
            throw cfem::PreconditionException(__FILE__, __LINE__, #expr); \
    } while(0)


#define CFEM_PRE_2(expr, msg) \
    do { \
        if (!(expr)) \
            throw cfem::PreconditionException(__FILE__, __LINE__, #expr, msg); \
    } while(0)

/**
 * @def CFEM_PRECONDITION
 * @brief Verifies a function precondition.
 */
#define CFEM_PRECONDITION(...) \
    CFEM_GET_MACRO(__VA_ARGS__, CFEM_PRE_2, CFEM_PRE_1)(__VA_ARGS__)


/* ================= POSTCONDITION ================= */

#define CFEM_POST_1(expr) \
    do { \
        if (!(expr)) \
            throw cfem::PostconditionException(__FILE__, __LINE__, #expr); \
    } while(0)


#define CFEM_POST_2(expr, msg) \
    do { \
        if (!(expr)) \
            throw cfem::PostconditionException(__FILE__, __LINE__, #expr, msg); \
    } while(0)

/**
 * @def CFEM_POSTCONDITION
 * @brief Verifies a function postcondition.
 */
#define CFEM_POSTCONDITION(...) \
    CFEM_GET_MACRO(__VA_ARGS__, CFEM_POST_2, CFEM_POST_1)(__VA_ARGS__)


/* ================= INVARIANT ================= */

#define CFEM_INV_1(expr) \
    do { \
        if (!(expr)) \
            throw cfem::InvariantException(__FILE__, __LINE__, #expr); \
    } while(0)


#define CFEM_INV_2(expr, msg) \
    do { \
        if (!(expr)) \
            throw cfem::InvariantException(__FILE__, __LINE__, #expr, msg); \
    } while(0)

/**
 * @def CFEM_INVARIANT
 * @brief Verifies an object invariant.
 */

#define CFEM_INVARIANT(...) \
    CFEM_GET_MACRO(__VA_ARGS__, CFEM_INV_2, CFEM_INV_1)(__VA_ARGS__)


/**
 * @def CFEM_CHECK_INVARIANTS
 * @brief Calls the invariant verification function.
 */
#define CFEM_CHECK_INVARIANTS() \
    checkInvariant()


#else   // Release mode


#define CFEM_ASSERT(...)         ((void)0)
#define CFEM_PRECONDITION(...)  ((void)0)
#define CFEM_POSTCONDITION(...) ((void)0)
#define CFEM_INVARIANT(...)     ((void)0)
#define CFEM_CHECK_INVARIANTS()  ((void)0)


#endif // NDEBUG


#ifdef NDEBUG
    #define CFEM_ERROR(msg) \
        throw std::runtime_error(std::string("[CFEM ERROR] ") + msg)

    #define CFEM_WARNING(msg) \
        std::cerr << "[CFEM WARNING] " << msg << std::endl
#else
    #define CFEM_ERROR(msg) \
        throw std::runtime_error(std::string("[CFEM ERROR] ") + msg + \
        " File: " + __FILE__ + ":" + std::to_string(__LINE__))

    #define CFEM_WARNING(msg) \
        std::cerr << "[CFEM WARNING] " << msg << \
        " (File: " << __FILE__ << ":" << __LINE__ << ")" << std::endl
#endif


#endif // CFEM_CONTRACT_EXCEPTION_H
