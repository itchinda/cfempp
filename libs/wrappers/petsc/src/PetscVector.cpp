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
 * @file PetscVector.cpp
 * @author your name (you@domain.com)
 * @brief 
 * @version 0.1
 * @date 2026-01-22
 * 
 * @copyright Copyright (c) 2026
 * 
 *************************************************************************************** */

#include "cfemutils.h"
#include "internal/petscUtils.h"
#include "PetscVector.h"

/** **************************************************************************************
 * @brief Construct a distributed PETSc vector.
 *
 * Collective on PETSC_COMM_WORLD.
 * @param local_size Number of local elements on this MPI rank
 * @param global_size Total number of elements across all MPI ranks
 *************************************************************************************** */
PetscVector::PetscVector(CfemInt local_size, CfemInt global_size) {
    cfem::Environment::ensureInitialized();
    PETSC_CHECK(VecCreate(PETSC_COMM_WORLD, &m_Vec));
    PETSC_CHECK(VecSetSizes(m_Vec, local_size, global_size));
    PETSC_CHECK(VecSetFromOptions(m_Vec));
    PETSC_CHECK(VecSetUp(m_Vec));
}

/** **************************************************************************************
 * @brief Destroy the PetscVector object and free PETSc resources.
 *        Ensures that VecDestroy is called to prevent memory leaks in the 
 *        underlying PETSc C library.
 *************************************************************************************** */
PetscVector::~PetscVector() {
    if (m_Vec) {
        PetscErrorCode ierr = VecDestroy(&m_Vec);
        (void)ierr; // ignore errors in destructor
    }
}

/**
 * @brief  Move constructor
 * 
 * @param other 
 */
PetscVector::PetscVector(PetscVector&& other) noexcept: m_Vec(other.m_Vec) {
    other.m_Vec = nullptr;
}

/**
 * @brief  Move assignment
 * 
 * @param other 
 * @return PetscVector& 
 */
PetscVector& PetscVector::operator=(PetscVector&& other) noexcept {
    if (this != &other) {
        if (m_Vec) PETSC_CHECK(VecDestroy(&m_Vec));
        m_Vec = other.m_Vec;
        other.m_Vec = nullptr;
    }
    return *this;
}

/** **************************************************************************************
 * @brief Set a single value in the vector
 * 
 * @param i Global index
 * @param val Value to set
 * @param mode INSERT_VALUES or ADD_VALUES
 *************************************************************************************** */
void PetscVector::setValue(CfemInt i, PetscScalar val, InsertMode mode) {
    PETSC_CHECK(VecSetValue(m_Vec, i, val, mode));
}

/** **************************************************************************************
 * @brief Assemble the vector.
 *
 * Collective on PETSC_COMM_WORLD.
 * 
 *************************************************************************************** */
void PetscVector::assemble() {
    PETSC_CHECK(VecAssemblyBegin(m_Vec));
    PETSC_CHECK(VecAssemblyEnd(m_Vec));
}

/** **************************************************************************************
 * @brief Set all entries to a constant value
 * 
 * @param val 
 *************************************************************************************** */
void PetscVector::setAll(PetscScalar val) {
    PETSC_CHECK(VecSet(m_Vec, val));
}

/** **************************************************************************************
 * @brief Compute the vector norm
 * 
 * @param type Norm type (default NORM_2)
 * @return PetscReal Norm value
 *************************************************************************************** */
PetscReal PetscVector::norm(NormType type) const {
    PetscReal val;
    PETSC_CHECK(VecNorm(m_Vec, type, &val));
    return val;
}

/** **************************************************************************************
 * @brief Perform the operation y = alpha * x + y.
 * 
 * Both vectors must be compatible (same global size and layout).
 *
 * @param alpha Scaling factor.
 * @param x Vector to be added.
 *************************************************************************************** */
void PetscVector::axpy(PetscScalar alpha, const PetscVector& x) {
    PETSC_CHECK(VecAXPY(m_Vec, alpha, x.getRaw()));
}

/** **************************************************************************************
 * @brief Return the number of local elements on this MPI rank.
 * 
 *************************************************************************************** */

PetscInt PetscVector::localSize() const {
    PetscInt n;
    PETSC_CHECK(VecGetLocalSize(m_Vec, &n));
    return n;
}

/** **************************************************************************************
 * @brief Return the full vector size (The total number of elements across all MPI ranks).
 * 
 *************************************************************************************** */

PetscInt PetscVector::globalSize() const {
    PetscInt n;
    PETSC_CHECK(VecGetSize(m_Vec, &n));
    return n;
}