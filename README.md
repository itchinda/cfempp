[![CFEM++ Project CI](https://github.com/itchinda/cfempp/actions/workflows/ci.yml/badge.svg)](https://github.com/itchinda/cfempp/actions/workflows/ci.yml)

# CFEM++: High-Performance Finite Element Solver

**CFEM++** is a research-grade C++ library designed for large-scale computational mechanics. It bridges the gap between complex geometric preprocessing in **Gmsh** and high-performance linear algebra via the **PETSc** framework.

Designed to handle **10^7+ Degrees of Freedom (DoFs)**, the project focuses on topological robustness, memory-efficient assembly, and high-order element support.

---

## Key Features

### 1. Robust Mesh Topology
* **Gmsh 4.1 Integration:** Native support for MSH ASCII formats with advanced vertex unique-ification and tag interning.
* **Dual-Layer Entity System:**
    * **Universe (ID -1):** Mandatory computational domain definition.
    * **SYS_BOUNDARY (ID -2):** Automated exterior facet detection to ensure 100% Boundary Condition (BC) coverage.
* **Corner-Based Mapping:** A specialized algorithm for quadratic elements (TET10, HEX27) that identifies physical entities by intersecting geometric vertices (corners), effectively solving the "mismatched diagonal" problem in high-order meshes.


### 2. Symbolic Expression Support
* **SymEngine Integration:** Define analytical solutions and source terms using natural mathematical syntax.
* **Just-In-Time (JIT) Compilation:** Leveraging **LLVM** to compile symbolic expressions into machine code at runtime for near-native evaluation speeds during assembly.
* **Automatic Differentiation:** High-performance `grad()`, `div()`, and `curl()` operators on finite element fields.


### 3. High-Performance Computing (HPC)
Leveraging **PETSc** for its linear algebra core, CFEM++ provides:
* **Distributed Computing:** MPI-parallelized operations.
* **Advanced Solvers:** Optimized for **Preconditioned Conjugate Gradient (PCG)** with **Algebraic Multigrid (GAMG)** preconditioning.
* **Scalability:** Capable of processing a large number of DOFs.



### 4. Comprehensive Element Library
| Dimension | Linear Elements | Quadratic Elements |
| :--- | :--- | :--- |
| **1D** | LINE2 | LINE3 |
| **2D** | TRI3, QUAD4 | TRI6, QUAD9 |
| **3D** | TET4, PRISM6, HEX8 | TET10, PRISM18, HEX27 |

---

## 🛠 Installation & Environment Setup

CFEM++ includes HPC-grade setup scripts that automate the configuration of PETSc, SymEngine, and system dependencies across Linux and macOS (Apple Silicon & Intel).

### Option A: Standard Users (Public Setup)
For researchers and users who want to safely configure the environment or build missing dependencies automatically.

**1. Source the public script:**
```bash
source scripts/setup_env_pub.sh release
```

**2. Missing dependencies?** Use the bootstrap flag to automatically build PETSc and SymEngine in your home directory without requiring root privileges:
```bash
source scripts/setup_env_pub.sh release --bootstrap
```

### Option B: Core Developers (Dev Setup)
For contributors. This script checks tools, automatically calls system package managers (apt, dnf, pacman, brew), and generates the CMakeUserPresets.

**1. Run the developer bootstrap script:**
```bash
./scripts/setup_env_dev.sh release
```

**2. (Optional) For heavy debugging**, use the debug mode:
```bash
./scripts/setup_env_dev.sh debug
```

### Build the Project
CFEM++ generates isolated build directories based on the chosen mode.

```bash
cd build_release
make -j$(nproc)
```

### Run Tests
```bash
# Execute unit tests across 4 MPI ranks
mpiexec -n 4 ./tests/unit_tests
```

## C++ API Example: Field Projection
The following example demonstrates how to generate a tetrahedral mesh, interpolate a symbolic JIT-compiled velocity field, calculate numerical errors (L2/H1 norms), and export the results for ParaView.

```c++
#include "cfemutils.h"
#include "Mesh.h"
#include "FEField.h"
#include "MathOperators.h"
#include "ScalarFunctionalEvaluator.h"
#include "LagrangeSpace.h"
#include "VTUExporter.h"

using namespace cfem;
using namespace cfem::sym;
using namespace cfem::sym::symengine;

int main(int argc, char **argv) {
    // Initialize CFEM++. This handles MPI, PETSc internal state
    cfem::Environment::initialize(argc, argv);

    // Generate Mesh
    auto mesh = Mesh::block(MeshOrder::Linear, 
                            {10, 10, 10}, 
                            Mesh::BlockBounds({0,0,0}, {1,1,1}), 
                            Mesh::BlockElementType::Tetrahedron);

    // Setup P2 Vector Space and Field
    auto V = LagrangeSpace::create(mesh, Order::P2, FieldType::Vector3D);
    auto u = FEField::create("u", V);
        
    // Symbolic Expression
    auto analytical_u = vector(sin(2 * CFEM_PI * x), cos(2 * CFEM_PI * y), sin(CFEM_PI * x * y * z));
    u->interpolate(analytical_u);

    // Error Analysis (Clean syntax, no explicit scheme params needed)
    CfemReal l2Norm_U, h1SemiNorm_U;
    IntegrationSchemeParams schemes_params{}; // Use default number of integration points
    ScalarFunctionalEvaluator integrator(mesh, schemes_params);
    integrator.L2Norm(u - analytical_u, &l2Norm_U, &h1SemiNorm_U);    
    CFEM_INFO << "L2-norm error on u: " << l2Norm_U;
    CFEM_INFO << "H1-semi norm error on u: " << h1SemiNorm_U;
        
    // Binary Export
    VTUExporter vtu_exporter(*mesh);
    vtu_exporter.addExpression("u", u);
    vtu_exporter.setFormat(VTUExporter::Format::BINARY);
    vtu_exporter.write("simulation_results");

    return 0; 
}
```
**From Preprocessing to Visualization**: Our `GmshParser` handles non-trivial manifolds with ease. Shown below is a 3D heart model imported into CFEM++ and visualized in ParaView. The seamless integration ensures that physical tags and boundary conditions remain consistent from the initial `.msh` file to the final visual results.

<p align="center">
  <img src="tests/grenier/heart.png" width="40%">
</p>
<p align="center">
  <em>Figure 1: Complex cardiac geometry imported from Gmsh and visualized in ParaView.</em>
</p>
```c++
// Future API (Expected)
auto u = Trial(V, "u");
auto v = Test(V, "v");
auto kappa = constant(1.0);

auto problem = VariationalProblem(mesh);
problem.addTerm( Integrate(kappa * grad(u) * grad(v), Volume) );

// Solver configuration
NewtonSolver solver; 
solver.setLinearSolver(PetscSolver); // ex: preconfigured PETSc solver (from PetscSolver)
solver.setTolerance(1e-8);

// Solve
solver.solve(problem);
```


## 📝 Metadata & License
- **Lead Developer**: itchinda
- **Version**: 0.1.0 (High-Order HPC Release)
- **Status**: Active Research
- **License**: Copyright (c) 2026 CFEM++ Team. All rights reserved.
