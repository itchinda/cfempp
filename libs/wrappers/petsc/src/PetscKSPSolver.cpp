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

/** **************************************************************************************
 * @file PetscKSPSolver.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2026-01-23
 * 
 * @copyright Copyright (c) 2026
 * 
 *************************************************************************************** */

#include "PetscKSPSolver.h"
#include "internal/petscUtils.h"
#include "cfemutils.h"

#include <iostream>

/** **************************************************************************************
 * @brief Construct a new Solver with an optional PETSc prefix.
 * @param prefix Optional string to group options (e.g., "outer_")
 *************************************************************************************** */
PetscKSPSolver::PetscKSPSolver(const std::string& prefix) : m_prefix(prefix) {
    cfem::Environment::ensureInitialized();
    PETSC_CHECK(KSPCreate(PETSC_COMM_WORLD, &m_ksp));
    
    if (!m_prefix.empty()) {
        // Ensure prefix ends with underscore for PETSc convention
        if (m_prefix.back() != '_') m_prefix += '_';
        PETSC_CHECK(KSPSetOptionsPrefix(m_ksp, m_prefix.c_str()));
    }
}

/** **************************************************************************************
 * @brief Destroy the Petsc KSP Solver object
 * 
 *************************************************************************************** */
PetscKSPSolver::~PetscKSPSolver() {
    if (m_ksp) {
        PetscErrorCode ierr = KSPDestroy(&m_ksp);
        (void)ierr; // ignore errors in destructor
    }
}

/** **************************************************************************************
 * @brief Move constructor
 * 
 * @param other 
 *************************************************************************************** */

PetscKSPSolver::PetscKSPSolver(PetscKSPSolver&& other) noexcept: m_ksp(other.m_ksp) {
    other.m_ksp = nullptr;
}

/** **************************************************************************************
 * @brief Move assignment
 * 
 * @param other 
 * @return PetscVector& 
 *************************************************************************************** */

PetscKSPSolver& PetscKSPSolver::operator=(PetscKSPSolver&& other) noexcept {
    if (this != &other) {
        if (m_ksp) PETSC_CHECK(KSPDestroy(&m_ksp));
        m_ksp = other.m_ksp;
        other.m_ksp = nullptr;
    }
    return *this;
}

/** **************************************************************************************
 * @brief Manually inject an option into the PETSc database for this solver.
 * @param name Option name (e.g., "ksp_type")
 * @param value Option value (e.g., "gmres")
 *************************************************************************************** */
void PetscKSPSolver::addOption(const std::string& name, const std::string& value) {
    std::string full_name = (name[0] == '-') ? name : "-" + name;
    
    // If a prefix is set and not present in the name, inject it
    if (!m_prefix.empty()) {
        if (full_name.rfind("-" + m_prefix, 0) != 0) {
            full_name.insert(1, m_prefix);
    }
}

    PETSC_CHECK(PetscOptionsSetValue(NULL, full_name.c_str(), value.c_str()));
}

/** **************************************************************************************
 * @brief Set the operators. Amat is used for the residual, Pmat for preconditioning.
 *************************************************************************************** */
void PetscKSPSolver::setOperators(PetscMatrix& Amat, PetscMatrix& Pmat) {
    PETSC_CHECK(KSPSetOperators(m_ksp, Amat.getRaw(), Pmat.getRaw()));
}

/** **************************************************************************************
 * @brief Standard case where the system matrix is also the preconditioning matrix.
 *************************************************************************************** */
void PetscKSPSolver::setOperators(PetscMatrix& A) {
    PETSC_CHECK(KSPSetOperators(m_ksp, A.getRaw(), A.getRaw()));
}

/** **************************************************************************************
 * @brief Solve the linear system Ax = b.
 * 
 * @param b 
 * @param x 
 *************************************************************************************** */
void PetscKSPSolver::solve(PetscVector& b, PetscVector& x) {
    PETSC_CHECK(KSPSetFromOptions(m_ksp));
    PETSC_CHECK(KSPSolve(m_ksp, b.getRaw(), x.getRaw()));
}

/** **************************************************************************************
 * @brief Return the number of iterations performed in the last solve.
 *
 * @note Valid only after solve().
 *
 * @return PetscInt Number of KSP iterations.
 *************************************************************************************** */
PetscInt PetscKSPSolver::getIterationCount() const {
    PetscInt its;
    PETSC_CHECK(KSPGetIterationNumber(m_ksp, &its));
    return its;
}

/** **************************************************************************************
 * @brief Return the reason why the last KSP solve terminated.
 *
 * This method queries PETSc for the convergence status of the most recent
 * call to solve(). The returned value indicates whether the solver:
 * - converged successfully,
 * - diverged,
 * - or stopped for another reason (e.g., iteration limit).
 *
 * PETSc convention:
 * - reason > 0 : solver converged
 * - reason = 0 : solver is still iterating (unlikely after solve())
 * - reason < 0 : solver diverged
 *
 * @note This method must be called after solve().
 *
 * @return KSPConvergedReason Enumeration describing the convergence outcome.
 *************************************************************************************** */
KSPConvergedReason PetscKSPSolver::getConvergedReason() const {
    KSPConvergedReason reason;
    PETSC_CHECK(KSPGetConvergedReason(m_ksp, &reason));
    return reason;
}

/**
 * @brief Check whether the last solve converged successfully.
 *
 * @return true if the solver converged (converged reason > 0).
 * @return false otherwise.
 */
bool PetscKSPSolver::hasConverged() const{
    return getConvergedReason() > 0;
}


/** **************************************************************************************
 * @brief Return the final residual norm from the last KSP solve.
 *
 * This method returns the norm of the residual associated with the
 * most recent call to solve(). The exact definition of the residual
 * (preconditioned or unpreconditioned) depends on the KSP type and
 * solver configuration.
 *
 * @note This value is only valid after a call to solve().
 * @note The norm returned here is the internally monitored norm used
 *       by PETSc for convergence testing.
 *
 * @return PetscReal Final residual norm.
 *************************************************************************************** */
PetscReal PetscKSPSolver::getResidualNorm() const {
    PetscReal rnorm;
    PETSC_CHECK(KSPGetResidualNorm(m_ksp, &rnorm));
    return rnorm;
}