#include <petsc.h>
#include "PetscMatrix.h"
#include "PetscVector.h"
#include "PetscKSPSolver.h"
#include <iostream>

int main(int argc, char **argv) {
    
    PetscCall(PetscInitialize(&argc, &argv, NULL, NULL));

    {
        const int n = 100000; // Global size of the system

        // 1. Setup the Matrix (Laplacian-like)
        PetscMatrix A(PETSC_DECIDE, n);
        
        std::pair<PetscInt, PetscInt> ownership_range {A.getOwnershipRange()};
        for (PetscInt i = ownership_range.first; i < ownership_range.second; i++) {
            A.setValue(i, i, 2.0); // Diagonal
            if (i > 0)     A.setValue(i, i - 1, -1.0); // Left neighbor
            if (i < n - 1) A.setValue(i, i + 1, -1.0); // Right neighbor
        }
        A.assemble();

        PetscVector b(PETSC_DECIDE, n);
        PetscVector x(PETSC_DECIDE, n);

        b.setAll(1.0); // Solve A*x = [1, 1, 1...]
        b.assemble();
        x.setAll(0.0); // Initial guess
        x.assemble();

        // 2. Petsc solver with options
        PetscKSPSolver solver("heat_");
        solver.addOption("ksp_type", "preonly");
        solver.addOption("pc_type", "lu");
        solver.addOption("ksp_view");
        solver.addOption("ksp_monitor");

        // 3. Finalize and Solve
        solver.setOperators(A);
        solver.solve(b, x);

        // 5. Check results
        PetscInt its = solver.getIterationCount();
        PetscReal final_norm = x.norm(NORM_2);

        PetscPrintf(PETSC_COMM_WORLD, "Solve finished in %d iterations.\n", its);
        PetscPrintf(PETSC_COMM_WORLD, "L2 Norm of solution: %g\n", (double)final_norm);

        // View the solution vector if the system is small
        if (n <= 20) {
            VecView(x.getRaw(), PETSC_VIEWER_STDOUT_WORLD);
        }
    } 

    PetscCall(PetscFinalize());
    return 0;
}