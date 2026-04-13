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

#ifndef CFEM_PETSC_UTILS_H
#define CFEM_PETSC_UTILS_H

#include <petscsys.h>

/**
 * @brief Check PETSc error codes and abort collectively on failure.
 *
 * This macro is safe to use inside:
 * - constructors
 * - destructors
 * - void functions
 *
 * It matches PETSc's fail-fast philosophy for MPI programs.
 */
#define PETSC_CHECK(call)                                     \
    do {                                                      \
        PetscErrorCode ierr = (call);                         \
        if (PetscUnlikely(ierr))                              \
            CHKERRABORT(PETSC_COMM_WORLD, ierr);              \
    } while (0)

#endif // CFEM_PETSC_UTILS_H
