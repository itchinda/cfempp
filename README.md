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
#include "cfemutils.h"
#include "Mesh.h"
#include "LagrangeSpace.h"
#include "FEField.h"
#include "VTKExporter.h"

using namespace cfem;

int main(int argc, char **argv)
{
    cfem::Environment::initialize(argc, argv);  // This internally call PetscInitialize and PetscFinalize when the program is completed

    auto mesh = Mesh::block(
        MeshOrder::Quadratic,
        Mesh::BlockSubdivisions{10, 10, 10}, // Mesh::BlockSubdivisions{5*(2<<i), 5*(2<<i), 5*(2<<i)},
        Mesh::BlockBounds(
            Vector3D{-1., -1., -1.},
            Vector3D{1., 1., 1.}
            ),
        Mesh::BlockElementType::Tetrahedron);

    // Lagrange spaces creation
    auto p1Space = LagrangeSpace::create(mesh, Order::P1, FieldType::Scalar);   // Pressure (Linear)
    auto p2Space = LagrangeSpace::create(mesh, Order::P2, FieldType::Vector3D); // Velocity (Quadratic)

    // fields creation
    auto pressure = FEField::create("P", p1Space);
    auto velocity = FEField::create("U", p2Space);

    const CfemReal PI = cfem::pi;

    auto analytical_p = [&](const Vector3D &X) -> CfemReal{
        return X.x * X.x + X.y * X.y + X.z * X.z + X.x * X.y + X.x * X.z + X.z * X.y + X.x + X.y + X.z + 1;
    };

    auto analytical_u = [&](const Vector3D &X) -> Vector3D{
        return {
            std::sin(2 * PI * X.x) * std::sin(2 * PI * X.y) * std::sin(2 * PI * X.z),
            std::cos(2 * PI * X.x) * std::cos(2 * PI * X.y) * std::cos(2 * PI * X.z),
            std::cos(2 * PI * X.x) * std::sin(2 * PI * X.y) * std::cos(2 * PI * X.z)};
    };

    // Field projection
    pressure->project(analytical_p);
    velocity->project(analytical_u);

    // Quadrature rule for l2-norm calculation
    IntegrationSchemeParams schemes_params{};
    schemes_params.NumIntgPtTriElem = 16;

    CfemReal lNorm_l2_p{FEField::L2NormDistance(*pressure, analytical_p, schemes_params)};
    CfemReal lNorm_l2_u{FEField::L2NormDistance(*velocity, analytical_u, schemes_params)};
    std::cout << "L2 error on the pressure: " << lNorm_l2_p  <<", L2 error on the velocity: " << lNorm_l2_u  << std::endl;

    // Exportation
    VTKExporter vtk_exporter(*mesh);

    vtk_exporter.addField(*pressure);
    vtk_exporter.addField(*velocity);

    vtk_exporter.write("block.vtk");

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
