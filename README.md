[![CFEM++ Project CI](https://github.com/itchinda/cfempp/actions/workflows/ci.yml/badge.svg)](https://github.com/itchinda/cfempp/actions/workflows/ci.yml)

# CFEM++: High-Performance Finite Element Solver

**CFEM++** is a research-grade C++ library designed for large-scale computational mechanics. It bridges the gap between complex geometric preprocessing in **Gmsh** and high-performance linear algebra via the **PETSc** framework.

Designed to handle **10^7+ Degrees of Freedom (DoFs)**, the project focuses on topological robustness, memory-efficient assembly, and high-order element support.

---

## 🚀 Key Features

### 1. Robust Mesh Topology
* **Gmsh 4.1 Integration:** Native support for MSH ASCII formats with advanced vertex unique-ification and tag interning.
* **Dual-Layer Entity System:**
    * **Universe (ID -1):** Mandatory computational domain definition.
    * **SYS_BOUNDARY (ID -2):** Automated exterior facet detection to ensure 100% Boundary Condition (BC) coverage.
* **Corner-Based Mapping:** A specialized algorithm for quadratic elements (TET10, HEX27) that identifies physical entities by intersecting geometric vertices (corners), effectively solving the "mismatched diagonal" problem in high-order meshes.


### 2. Advanced Symbolic Engine (JIT Accelerated)
* **SymEngine Integration:** Define analytical solutions and source terms using natural mathematical syntax.
* **Just-In-Time (JIT) Compilation:** Leveraging **LLVM** to compile symbolic expressions into machine code at runtime for near-native evaluation speeds during assembly.
* **Automatic Differentiation:** High-performance `grad()`, `div()`, and `curl()` operators on finite element fields.


### 3. High-Performance Computing (HPC)
Leveraging **PETSc** for its linear algebra core, CFEM++ provides:
* **Distributed Computing:** MPI-parallelized operations using `Vec` and `Mat` (AIJ format).
* **Advanced Solvers:** Optimized for **Preconditioned Conjugate Gradient (PCG)** with **Algebraic Multigrid (GAMG)** preconditioning.
* **Scalability:** Capable of processing ~8M nodes and ~6M cells in under 30 seconds on standard workstation hardware (32GB RAM).



### 4. Comprehensive Element Library
| Dimension | Linear Elements | Quadratic Elements |
| :--- | :--- | :--- |
| **1D** | LINE2 | LINE3 |
| **2D** | TRI3, QUAD4 | TRI6, QUAD9 |
| **3D** | TET4, PRISM6, HEX8 | TET10, PRISM18, HEX27 |

---

## 🛠 Installation & Automated Setup

CFEM++ includes a **Smart Setup Script** that automates the installation of PETSc, SymEngine, and system dependencies.

### 1. Initialize Environment
Run the setup script from the project root. This will verify your compiler, install required `apt` packages, and configure dependencies.

```bash
# Recommended: Standard Release mode
source scripts/setup_env.sh release

# For development: Debug mode with heavy error-checking
source scripts/setup_env.sh debug
         
```

### 2. Build Project
```bash
cd build
make -j$(nproc)
```

### 3. Run Tests
```bash
# Execute unit tests across 4 MPI ranks
mpiexec -n 4 ./tests/unit_tests
```

## 💻 C++ API Example: Field Projection
The following example demonstrates how to generate a tetrahedral mesh, project a symbolic JIT-compiled velocity field, calculate numerical errors (L2/H1 norms), and export the results for ParaView.

```c++
#include "cfemutils.h"
#include "Mesh.h"
#include "FEField.h"
#include "MathOperators.h"
#include "DomainIntegrator.h"
#include "VTUExporter.h"

using namespace cfem;
using namespace cfem::sym;

int main(int argc, char **argv) {
    // 1. One call to rule them all: Handles MPI, PETSc, and CFEM internal state
    cfem::Environment::initialize(argc, argv);

    // 2. Generate Mesh
    auto mesh = Mesh::block(MeshOrder::Linear, 
                            {10, 10, 10}, 
                            Mesh::BlockBounds({0,0,0}, {1,1,1}), 
                            Mesh::BlockElementType::Tetrahedron);

    // 3. Setup P2 Vector Space and Field
    auto velocitySpace = LagrangeSpace::create(mesh, Order::P2, FieldType::Vector3D);
    auto velocity = FEField::create("U", velocitySpace);
        
    // 4. JIT-compiled Symbolic Definition
    auto analytical_u = vector(sin(2 * pi * x), cos(2 * pi * y), sin(pi * x * y * z));

    // 5. High-level Projection
    velocity->project(analytical_u);

    // 6. Error Analysis (Clean syntax, no explicit scheme params needed)
    CfemReal l2Norm_U, h1SemiNorm_U;
    DomainIntegrator::L2Norm(mesh, analytical_u - velocity, &l2Norm_U, &h1SemiNorm_U);
    
    CFEM_INFO << "L2-norm error on U: " << l2Norm_U;
    CFEM_INFO << "H1-semi norm error on U: " << h1SemiNorm_U;
        
    // 7. High-performance Binary Export
    VTUExporter vtu_exporter(*mesh);
    vtu_exporter.addExpression("U", velocity);
    vtu_exporter.setFormat(VTUExporter::Format::BINARY);
    vtu_exporter.write("simulation_results");

    return 0; 
    // cfem::Environment handles finalization automatically on scope exit or via atexit
}
```
**From Preprocessing to Visualization**: Our robust `GmshParser` handles non-trivial manifolds with ease. Shown below is a 3D heart model imported into CFEM++ and visualized in ParaView. The seamless integration ensures that physical tags and boundary conditions remain consistent from the initial `.msh` file to the final visual results.

<div align="center">

![Figure 1: Complex cardiac geometry imported from Gmsh and visualized in ParaView.](tests/grenier/heart.png)

*Figure 1: Complex cardiac geometry imported from Gmsh and visualized in ParaView.*

</div>

<p align="center">
  <img src="tests/grenier/heart.png" width="40%">
</p>
<p align="center">
  <em>Figure 1: Complex cardiac geometry imported from Gmsh and visualized in ParaView.</em>
</p>


## 📝 Metadata & License
- **Lead Developer**: itchinda
- **Version**: 0.1.0 (High-Order HPC Release)
- **Status**: Active Research
- **License**: Copyright (c) 2026 CFEM++ Team. All rights reserved.
