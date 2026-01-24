#ifndef PETSC_UTILS_H
#define PETSC_UTILS_H

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

#endif // ODEPP_PETSC_UTILS_H
