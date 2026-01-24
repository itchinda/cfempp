/** **************************************************************************************
 * @file PetscKSPSolver.h
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2026-01-23
 * 
 * @copyright Copyright (c) 2026
 * 
 *************************************************************************************** */

#ifndef PETSCKSPSOLVER_H
#define PETSCKSPSOLVER_H

#include <petscksp.h>
#include <string>
#include "PetscMatrix.h"
#include "PetscVector.h"

/** **************************************************************************************
 * @class PetscKSPSolver
 * @brief RAII wrapper for PETSc KSP linear solvers.
 *
 * Supports optional prefix-based configuration through the PETSc
 * options database, enabling runtime configuration via command-line
 * or programmatic injection.
 *
 * Example:
 * @code
 * PetscKSPSolver solver("inner");
 * solver.addOption("ksp_type", "gmres");
 * solver.addOption("pc_type", "ilu");
 * solver.setFromOptions();
 * solver.setOperators(A);
 * solver.solve(b, x);
 * @endcode
 *************************************************************************************** */

class PetscKSPSolver {
public:
    PetscKSPSolver(const std::string& prefix = "");
    ~PetscKSPSolver();

    PetscKSPSolver(const PetscKSPSolver&) = delete;
    PetscKSPSolver& operator=(const PetscKSPSolver&) = delete;

    PetscKSPSolver(PetscKSPSolver&& other) noexcept;
    PetscKSPSolver& operator=(PetscKSPSolver&& other) noexcept;

    void addOption(const std::string& name, const std::string& value = "");
    void setOperators(PetscMatrix& Amat, PetscMatrix& Pmat);
    void setOperators(PetscMatrix& A);

    void solve(PetscVector& b, PetscVector& x);

    PetscInt getIterationCount() const;
    KSP getRaw() { return m_ksp; }

    KSPConvergedReason getConvergedReason() const;
    bool hasConverged() const;
    PetscReal getResidualNorm() const;

private:
    KSP m_ksp {nullptr};
    std::string m_prefix;
};

#endif