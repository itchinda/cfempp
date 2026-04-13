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

/**
 * @file PetscVector.h
 * @author Ismael Tchinda Ngueyong
 * @brief RAII wrapper around PETSc Vec
 * @version 0.1
 * @date 2026-01-23
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#ifndef PETSCVECTOR_H
#define PETSCVECTOR_H

#include <petscvec.h>


/** **************************************************************************************
 * @class PetscVector
 * @brief RAII wrapper around PETSc Vec.
 *
 * This class manages the lifetime of a distributed PETSc vector.
 * It follows PETSc ownership semantics and must be used collectively
 * on PETSC_COMM_WORLD.
 *
 * Copying is disabled to prevent double-free; move semantics are supported.
 *
 * Typical usage:
 * @code
 * PetscVector x(PETSC_DECIDE, n);
 * x.setAll(0.0);
 * x.assemble();
 * @endcode
 *************************************************************************************** */

class PetscVector {
public:
    PetscVector(int local_size, int global_size);
    ~PetscVector();

    PetscVector(PetscVector&& other) noexcept;
    PetscVector& operator=(PetscVector&& other) noexcept;


    PetscVector(const PetscVector&) = delete;
    PetscVector& operator=(const PetscVector&) = delete;

    Vec&       getRaw()       { return m_Vec; }
    const Vec& getRaw() const { return m_Vec; }
    
    void setValue(int i, PetscScalar val, InsertMode mode = INSERT_VALUES);
    void assemble();

    // Helpful math wrappers
    void setAll(PetscScalar val);
    PetscReal norm(NormType type = NORM_2) const;
    void axpy(PetscScalar alpha, const PetscVector& x); // y = alpha*x + y

    PetscInt globalSize() const;
    PetscInt localSize() const;
private:
    Vec m_Vec {nullptr};
};

#endif