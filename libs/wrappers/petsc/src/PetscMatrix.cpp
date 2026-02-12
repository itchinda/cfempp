/** **************************************************************************************
 * @file PetscMatrix.cpp
 * @author itchinda
 * @brief Implementation of the PetscMatrix RAII wrapper for PETSc Mat objects.
 * @version 0.1
 * @date 2026-01-22
 * @copyright Copyright (c) 2026
 *************************************************************************************** */

#include "PetscMatrix.h"
#include "internal/petscUtils.h"

#include "cfemutils.h"

/** **************************************************************************************
 * @brief Construct a new PetscMatrix and allocate memory.
 *        This constructor initializes a square sparse matrix distributed across 
 *        the PETSC_COMM_WORLD communicator.
 * @param local_rows Number of rows to be stored on the local processor.
 * @param global_rows Total number of rows across all processors.
 *************************************************************************************** */

 PetscMatrix::PetscMatrix(int local_rows, int global_rows) {
    cfem::Environment::ensureInitialized();
    CFEM_PRECONDITION(global_rows>=0);
    PETSC_CHECK(MatCreate(PETSC_COMM_WORLD, &m_Matrix));
    PETSC_CHECK(MatSetSizes(m_Matrix, local_rows, local_rows, global_rows, global_rows));
    PETSC_CHECK(MatSetFromOptions(m_Matrix));
}

/** **************************************************************************************
 * @brief Destroy the matrix (never throws)
 *************************************************************************************** */

 PetscMatrix::~PetscMatrix() {
    if (m_Matrix) {
        PetscErrorCode ierr = MatDestroy(&m_Matrix);
        (void)ierr; // ignore errors in destructor
    }
}

/** **************************************************************************************
 * @brief Move constructor
 * 
 * @param other 
 *************************************************************************************** */

PetscMatrix::PetscMatrix(PetscMatrix&& other) noexcept: m_Matrix(other.m_Matrix) {
    other.m_Matrix = nullptr;
}

/** **************************************************************************************
 * @brief Move assignment operator
 * 
 * @param other 
 * @return PetscMatrix& 
 *************************************************************************************** */

PetscMatrix& PetscMatrix::operator=(PetscMatrix&& other) noexcept {
    if (this != &other) {
        if (m_Matrix) PETSC_CHECK(MatDestroy(&m_Matrix));
        m_Matrix = other.m_Matrix;
        other.m_Matrix = nullptr;
    }
    return *this;
}

/** **************************************************************************************
 * @brief Returns a reference to the underlying PETSc Mat object.
 * @return Mat& The raw PETSc Matrix handle.
 *************************************************************************************** */
Mat& PetscMatrix::getRaw() { 
    return m_Matrix; 
}

/** **************************************************************************************
 * @brief Returns a const reference to the underlying PETSc Mat object.
 * @return const Mat& The raw PETSc Matrix handle.
 *************************************************************************************** */
const Mat& PetscMatrix::getRaw() const { 
    return m_Matrix; 
}

/**
 * @brief Preallocate matrix storage (AIJ)
 *
 * @param d_nz Number of nonzeros per row on-diagonal
 * @param o_nz Number of nonzeros per row off-diagonal
 */
void PetscMatrix::preallocate(PetscInt d_nz, PetscInt o_nz) {
    PETSC_CHECK(MatMPIAIJSetPreallocation(
        m_Matrix, d_nz, nullptr, o_nz, nullptr));
    PETSC_CHECK(MatSeqAIJSetPreallocation(
        m_Matrix, d_nz, nullptr));

    PETSC_CHECK(MatSetUp(m_Matrix));
}

/** **************************************************************************************
 * @brief Sets a single value in the matrix.
 *        Note: Values are cached locally. You must call assemble() before the 
 *        matrix can be used in computations.
 * @param row The global row index.
 * @param col The global column index.
 * @param value The scalar value to insert.
 * @param mode Either INSERT_VALUES or ADD_VALUES.
 *************************************************************************************** */
void PetscMatrix::setValue(int row, int col, PetscScalar value, InsertMode mode) {
    PETSC_CHECK(MatSetValue(m_Matrix, row, col, value, mode));
}

/** **************************************************************************************
 * @brief Synchronizes the matrix across all processors.
 *        This is a collective operation that performs the necessary MPI 
 *        communication to move entries to their home processors and finalizes 
 *        the sparse matrix structure.
 *************************************************************************************** */
void PetscMatrix::assemble() {
    PETSC_CHECK(MatAssemblyBegin(m_Matrix, MAT_FINAL_ASSEMBLY));
    PETSC_CHECK(MatAssemblyEnd(m_Matrix, MAT_FINAL_ASSEMBLY));
    PETSC_CHECK(MatSetOption(m_Matrix, MAT_NEW_NONZERO_LOCATIONS, PETSC_FALSE));
}

/** **************************************************************************************
 * @brief Return the number of local rows and columns on this MPI rank.
 * 
 *************************************************************************************** */

std::pair<PetscInt, PetscInt> PetscMatrix::getLocalSize() const {
    PetscInt n_rows;
    PetscInt n_cols;
    PETSC_CHECK(MatGetLocalSize(m_Matrix, &n_rows, &n_cols));
    return {n_rows, n_cols};
}

/** **************************************************************************************
 * @brief Return the number of global rows and columns across all MPI ranks in a std::pair.
 * 
 *************************************************************************************** */

std::pair<PetscInt, PetscInt> PetscMatrix::getGlobalSize() const {
    PetscInt n_rows;
    PetscInt n_cols;
    PETSC_CHECK(MatGetSize(m_Matrix, &n_rows, &n_cols));
    return {n_rows, n_cols};
}

/** **************************************************************************************
 * @brief Ownership range [rstart, rend)
 *************************************************************************************** */

std::pair<PetscInt, PetscInt> PetscMatrix::getOwnershipRange() const {
    PetscInt rstart, rend;
    PETSC_CHECK(MatGetOwnershipRange(m_Matrix, &rstart, &rend));
    return {rstart, rend};
}