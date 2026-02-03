[![ODEPP Project CI](https://github.com/itchinda/odepp/actions/workflows/ci.yml/badge.svg)](https://github.com/itchinda/odepp/actions/workflows/ci.yml)

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



### 2. High-Performance Computing (HPC)
Leveraging **PETSc** for its linear algebra core, CFEM++ provides:
* **Distributed Computing:** MPI-parallelized operations using `Vec` and `Mat` (AIJ format).
* **Advanced Solvers:** Optimized for **Preconditioned Conjugate Gradient (PCG)** with **Algebraic Multigrid (GAMG)** preconditioning.
* **Scalability:** Capable of processing ~8M nodes and ~6M cells in under 30 seconds on standard workstation hardware (32GB RAM).



### 3. Comprehensive Element Library
| Dimension | Linear Elements | Quadratic Elements |
| :--- | :--- | :--- |
| **1D** | LINE2 | LINE3 |
| **2D** | TRI3, QUAD4 | TRI6, QUAD9 |
| **3D** | TET4, PRISM6, HEX8 | TET10, PRISM18, HEX27 |

---

## 🛠 Installation & Dependencies

Ensure you have the following installed:
* **C++17 Compiler** (GCC, Clang, or MSVC)
* **CMake** (3.15+)
* **PETSc** (configured with MPI and HYPRE/MUMPS recommended)
* **Gmsh** (for mesh generation)

```bash
git clone [https://github.com/username/cfempp.git](https://github.com/username/cfempp.git)
cd cfempp
mkdir build && cd build
cmake ..
make -j
```

```c++
#include "Mesh.h"
#include "GmshParser.h"
#include "VtkExporter.h"

int main() {
    Mesh mesh;
    GmshParser parser;
    VTKExporter exporter(mesh);

    // Load MSH 4.1 file and build internal topology
    if (parser.load(mesh, "gmsh_file.msh")) {
        
        // Display diagnostic report (Memory, Bounding Box, Entities)
        mesh.info(); 

        // Export to vtk format for paraview visualization
        exporter.setExportQuality(true);
        exporter.write("hex.vtk");

        // Next steps: Finite Element Function Spaces / Fields / Assembly / Solve

    }
    return 0;
}

```

## 📊 Performance Benchmark

- **Hardware**: Intel Core Ultra 7 (165U), 32GB RAM.
- **Workload**: 7.9M nodes, 5.8M cells (3D Quadratic).
- **Load Time**: ~29 seconds (ASCII Parsing + Topological Reconstruction).

## 📝 Metadata & License
- **Lead Developer**: itchinda
- **Version**: 0.1.0 (High-Order HPC Release)
- **Status**: Active Research
- **License**: Copyright (c) 2026 CFEM++ Team. All rights reserved.
