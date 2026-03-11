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



### 3. Comprehensive Element Library
| Dimension | Linear Elements | Quadratic Elements |
| :--- | :--- | :--- |
| **1D** | LINE2 | LINE3 |
| **2D** | TRI3, QUAD4 | TRI6, QUAD9 |
| **3D** | TET4, PRISM6, HEX8 | TET10, PRISM18, HEX27 |

---

## 🛠 Installation & Dependencies

CFEM++ requires a specific stack for symbolic JIT and HPC linear algebra.

### System Prerequisites
```bash
# Ubuntu/Debian example
sudo apt-get update
sudo apt-get install -y build-essential gfortran cmake python3 m4 wget \
        libgmp-dev libmpfr-dev libmpc-dev libzstd-dev llvm-dev \
        libedit-dev libcurl4-openssl-dev libomp-dev
sudo apt-get install libpetsc-dev libsymengine-dev libgmp-dev libmpfr-dev libomp-dev
```


Ensure you have the following installed:
* **C++20 Compiler** (GCC, Clang, or MSVC)
* **CMake** (3.15+)
* **PETSc** (configured with MPI and HYPRE/MUMPS recommended)
* **SymEngine** (Symbolic math and JIT compilation)
* **OpenMP** (Parallel multi-threading support, usually bundled with GCC/Clang)


```bash
git clone [https://github.com/itchinda/cfempp.git](https://github.com/itchinda/cfempp.git)
cd cfempp
mkdir build && cd build

# Point to your PETSc and SymEngine installations if they are in non-standard paths
cmake .. -DCMAKE_BUILD_TYPE=Release \
         -DPETSC_DIR=$PETSC_DIR \
         -DPETSC_ARCH=$PETSC_ARCH
         
make -j$(nproc)
```

```c++
#include "cfemutils.h"
#include "Vector3D.h"
#include "Mesh.h"
#include "LagrangeSpace.h"
#include "FEField.h"
#include "SymbolicExpression.h"
#include "DomainIntegrator.h"
#include "VTKExporter.h"

#include <iostream>
#include <cmath>
#include <memory>

using namespace cfem;
using namespace cfem::sym;
static char help[] = "CFEM Example: Field projection on a quadratic mesh.\n\n";

int main(int argc, char **argv)
{
    cfem::Environment::initialize(argc, argv, help);  // Internally handle Petsc initialization and finalization

    auto mesh = Mesh::block(
        MeshOrder::Linear,
        Mesh::BlockSubdivisions{10, 10, 10},
        Mesh::BlockBounds(
            Vector3D{0., 0., 0.},
            Vector3D{1., 1., 1.}),
        Mesh::BlockElementType::Edge);

    auto pressureSpace = LagrangeSpace::create(mesh, Order::P1, FieldType::Scalar);
    auto velocitySpace = LagrangeSpace::create(mesh, Order::P2, FieldType::Vector3D);

    auto pressure = FEField::create("P", pressureSpace);
    auto velocity = FEField::create("U", velocitySpace);
        
    auto analytical_p = scalar(x * x + y * y + z * z + x * y + x * z + z * y + x + y + z + 1); // SymEngine expressions
    auto analytical_u = vector(sin(2 * pi * x) * sin(2 * pi * y) * sin(2 * pi * z),
                               cos(2 * pi * x) * cos(2 * pi * y) * cos(2 * pi * z),
                               cos(2 * pi * x) * sin(2 * pi * y) * cos(2 * pi * z));

    // Nodal Projection
    pressure->project(*analytical_p);
    velocity->project(*analytical_u);

    // Ability to perform mathematical differentiation grad(), div()
    auto fe_gradp = grad(pressure);

    // L2- and H1 semi- norms
    IntegrationSchemeParams schemes_params{};
    CfemReal l2Norm_P, l2Norm_U, h1SemiNorm_P, h1SemiNorm_U;
    DomainIntegrator::L2Norm(pressureSpace, analytical_p - pressure, &l2Norm_P, &h1SemiNorm_P, schemes_params);
    DomainIntegrator::L2Norm(velocitySpace, analytical_u - velocity, &l2Norm_U, &h1SemiNorm_U, schemes_params);
    std::cout << "L2-norm error on P: " << l2Norm_P << " H1-Seminorm error on P: " << h1SemiNorm_P << std::endl;
    std::cout << "L2-norm error on U: " << l2Norm_U << " H1-Seminorm error on U: " << h1SemiNorm_U << std::endl;
        
    // Exportation to  VTU/VTK files (ASCII and BINARY formats supported)
    VTKExporter vtk_exporter(*mesh);
    vtk_exporter.setExportedEntities( {"xmin", "xmax", "ymin", "ymax", "zmin", "zmax"});
    vtk_exporter.exportCellQuality(true);
    vtk_exporter.addExpression("P", pressure);
    vtk_exporter.addExpression("P_exact", analytical_p);
    vtk_exporter.addExpression("U", velocity);
    vtk_exporter.addExpression("U_exact", analytical_u);
    vtk_exporter.addExpression("GradP", fe_gradp);
    vtk_exporter.setFormat(VTKExporter::Format::ASCII);
    vtk_exporter.write("cube");

    return 0;
}

```
**From Preprocessing to Visualization**: Our robust `GmshParser` handles non-trivial manifolds with ease. Shown below is a 3D heart model imported into CFEM++ and visualized in ParaView. The seamless integration ensures that physical tags and boundary conditions remain consistent from the initial `.msh` file to the final visual results.

<p align="center">
  <img src="tests/grenier/heart.png" width="600px" alt="Human Heart Mesh in ParaView">
  <br>
  <i>Figure 1: Complex cardiac geometry imported from Gmsh and visualized in ParaView.</i>
</p>

## 📊 Performance Benchmark

- **Hardware**: Intel Core Ultra 7 (165U), 32GB RAM.
- **Workload**: 7.9M nodes, 5.8M cells (3D Quadratic).
- **Load Time**: ~29 seconds (ASCII Parsing + Topological Reconstruction).

## 📝 Metadata & License
- **Lead Developer**: itchinda
- **Version**: 0.1.0 (High-Order HPC Release)
- **Status**: Active Research
- **License**: Copyright (c) 2026 CFEM++ Team. All rights reserved.
