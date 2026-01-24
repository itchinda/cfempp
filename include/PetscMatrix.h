/**
 * @file PetscMatrix.h
 * @author Ismael Tchinda Ngueyong
 * @brief RAII wrapper around PETSc Mat
 * @version 0.1
 * @date 2026-01-23
 * 
 * @copyright Copyright (c) 2026
 * 
 */

#ifndef PETSCMATRIX_H
#define PETSCMATRIX_H

#include <petscmat.h>
#include <utility>

static constexpr PetscInt DEFAULT_D_NZ = 10;
static constexpr PetscInt DEFAULT_O_NZ = 10;

/** **************************************************************************************
 * @class PetscMatrix
 * @brief RAII wrapper around PETSc Mat.
 *
 * This class manages the lifetime of a distributed PETSc Matrix.
 * It follows PETSc ownership semantics and must be used collectively
 * on PETSC_COMM_WORLD.
 *
 * Copying is disabled to prevent double-free; move semantics are supported.
 *
 * Typical usage:
 * @code
 * PetscInt n {10};
 * PetscMatrix A(PETSC_DECIDE, n);
 * A.preallocate(3,3)
 * PetscInt rstart, rend;
 * MatGetOwnershipRange(A.getRaw(), &rstart, &rend);
 * for (PetscInt i = rstart; i < rend; ++i)
 *      A.setValue(i, i, 2.0);
 *      if (i>0) A.setValue(i, i-1, -1.0);
 *      if (i<n-1) A.setValue(i, i+1, -1.0);
 * A.assemble();
 * @endcode
 *************************************************************************************** */

class PetscMatrix {
public:
    explicit PetscMatrix(int local_rows, int global_rows);
    ~PetscMatrix();

    PetscMatrix(const PetscMatrix&) = delete;
    PetscMatrix& operator=(const PetscMatrix&) = delete;

    PetscMatrix(PetscMatrix&& other) noexcept;
    PetscMatrix& operator=(PetscMatrix&& other) noexcept;

    Mat&       getRaw();
    const Mat& getRaw() const;
    void setValue(int row, int col, PetscScalar value, InsertMode mode = INSERT_VALUES);
    void assemble();

    /// Preallocate AIJ storage
    void preallocate(PetscInt d_nz = DEFAULT_D_NZ,
                     PetscInt o_nz = DEFAULT_O_NZ);

    std::pair<PetscInt, PetscInt> getGlobalSize() const;
    std::pair<PetscInt, PetscInt> getLocalSize() const;

    /// Return ownership range [rstart, rend)
    std::pair<PetscInt, PetscInt> getOwnershipRange() const;

private:
    Mat m_Matrix {nullptr};
};

#endif