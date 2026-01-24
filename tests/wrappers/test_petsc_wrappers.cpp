/** **************************************************************************************
 * @file test_petsc_solver.cpp
 * @author Ismael Tchinda Ngueyong
 * @brief Unit tests for PetscVector, PetscMatrix, and PetscKSPSolver wrappers.
 * @version 0.1
 * @date 2026-01-23
 * 
 * @copyright Copyright (c) 2026
 *************************************************************************************** */

#include <gtest/gtest.h>
#include "PetscKSPSolver.h"
#include <cmath>

//-------------------------------------------------------------
// Linear solver correctness: Identity matrix
//-------------------------------------------------------------
TEST(SolverTest, IdentitySolutionCorrectness) {
    const PetscInt n = 10;
    PetscMatrix A(PETSC_DECIDE, n);
    A.preallocate(1, 0);

    auto [rstart, rend] = A.getOwnershipRange();
    for (PetscInt i = rstart; i < rend; ++i)
        A.setValue(i, i, 1.0);
    A.assemble();

    PetscVector b(PETSC_DECIDE, n);
    PetscVector x(PETSC_DECIDE, n);
    b.setAll(5.0);
    b.assemble();
    x.setAll(0.0);
    x.assemble();

    PetscKSPSolver solver;
    solver.setOperators(A);
    solver.solve(b, x);

    EXPECT_TRUE(solver.hasConverged());

    // Compute residual r = Ax - b
    PetscVector r(PETSC_DECIDE, n);
    r.setAll(0.0);
    r.assemble();
    MatMult(A.getRaw(), x.getRaw(), r.getRaw());
    r.axpy(-1.0, b);

    EXPECT_LT(r.norm(NORM_2), 1e-10);
}

//-------------------------------------------------------------
// 1D Laplacian test
//-------------------------------------------------------------
TEST(SolverTest, Laplacian1D) {
    const PetscInt n = 10;
    PetscMatrix A(PETSC_DECIDE, n);
    A.preallocate(3, 3);

    auto [rstart, rend] = A.getOwnershipRange();
    for (PetscInt i = rstart; i < rend; ++i) {
        A.setValue(i, i, 2.0);
        if (i > 0) A.setValue(i, i - 1, -1.0);
        if (i < n - 1) A.setValue(i, i + 1, -1.0);
    }
    A.assemble();

    PetscVector b(PETSC_DECIDE, n);
    PetscVector x(PETSC_DECIDE, n);
    b.setAll(1.0);
    b.assemble();
    x.setAll(0.0);
    x.assemble();

    PetscKSPSolver solver;
    solver.setOperators(A);
    solver.solve(b, x);

    EXPECT_TRUE(solver.hasConverged());

    PetscVector r(PETSC_DECIDE, n);
    r.setAll(0.0);
    r.assemble();
    MatMult(A.getRaw(), x.getRaw(), r.getRaw());
    r.axpy(-1.0, b);

    EXPECT_LT(r.norm(NORM_2), 1e-10);
}

//-------------------------------------------------------------
// Vec axpy test
//-------------------------------------------------------------
TEST(VectorTest, AxpyWorks) {
    PetscVector x(PETSC_DECIDE, 5);
    PetscVector y(PETSC_DECIDE, 5);

    x.setAll(2.0);
    y.setAll(1.0);
    x.assemble();
    y.assemble();

    y.axpy(3.0, x); // y = 3*x + y = 7

    EXPECT_NEAR(y.norm(NORM_2), std::sqrt(5 * 49.0), 1e-12);
}

//-------------------------------------------------------------
// Diagonal matrix multiplication test
//-------------------------------------------------------------
TEST(MatrixTest, DiagonalMatMult) {
    const PetscInt n = 4;
    PetscMatrix A(PETSC_DECIDE, n);
    A.preallocate(1, 0);

    auto [rstart, rend] = A.getOwnershipRange();
    for (PetscInt i = rstart; i < rend; ++i)
        A.setValue(i, i, i + 1);
    A.assemble();

    PetscVector x(PETSC_DECIDE, n);
    PetscVector y(PETSC_DECIDE, n);

    x.setAll(1.0);
    y.setAll(0.0);
    x.assemble();
    y.assemble();

    MatMult(A.getRaw(), x.getRaw(), y.getRaw());

    EXPECT_NEAR(y.norm(NORM_2), std::sqrt(1 + 4 + 9 + 16), 1e-12);
}
