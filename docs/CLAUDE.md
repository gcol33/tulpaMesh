# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working
with code in this repository.

## What is tulpaMesh

R package for constrained Delaunay triangulation mesh generation and FEM
matrix assembly for spatial SPDE models. Hard dependency of the tulpa
Bayesian hierarchical modeling engine. Built on the CDT header-only C++
library (vendored in `src/cdt/`, MPL-2.0, documented in
`inst/COPYRIGHTS`).

## Build & Development Commands

``` bash
# Load package (compile C++ and source R)
"C:/Program Files/R/R-4.5.2/bin/Rscript.exe" -e 'devtools::load_all()'

# Run tests (141 tests)
"C:/Program Files/R/R-4.5.2/bin/Rscript.exe" -e 'devtools::test()'

# Run single test file
"C:/Program Files/R/R-4.5.2/bin/Rscript.exe" -e 'devtools::test(filter = "sphere")'

# Full R CMD check
"C:/Program Files/R/R-4.5.2/bin/Rscript.exe" -e 'devtools::check(args = "--no-manual")'

# After editing C++ files: regenerate Rcpp exports then rebuild docs
"C:/Program Files/R/R-4.5.2/bin/Rscript.exe" -e 'Rcpp::compileAttributes(); devtools::document()'
```

## Architecture

    User coordinates + boundary
            │
            ▼
    tulpa_mesh()  [R/mesh.R]          Formula parsing, convex hull extension,
            │                          boundary constraints, grid refinement,
            │                          deduplication, Ruppert refinement
            ├─► cpp_triangulate()      CDT constrained Delaunay
            └─► cpp_ruppert_refine()   Iterative circumcenter insertion
            │
            ▼
    tulpa_mesh S3 object               vertices, triangles, edges (all 1-based)
            │
            ▼
    fem_matrices()  [R/mesh.R]        Sparse matrix assembly
            │
            ├─► cpp_fem_matrices()     2D P1 FEM: mass (C) + stiffness (G)
            ├─► cpp_fem_matrices_3d()  3D surface FEM (spherical meshes)
            ├─► cpp_projection_matrix()   Grid-accelerated barycentric (2D)
            └─► cpp_projection_matrix_3d() 3D surface projection
            │
            ▼
    List(C, G, A, ...)                dgCMatrix sparse matrices

### R source files

| File | Purpose |
|----|----|
| `R/mesh.R` | [`tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh.md), [`fem_matrices()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices.md), sf/boundary helpers |
| `R/sphere.R` | [`tulpa_mesh_sphere()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh_sphere.md), [`print.tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh.md), icosahedral subdivision |
| `R/diagnostics.R` | [`mesh_quality()`](https://github.com/gcol33/tulpaMesh/reference/mesh_quality.md), [`mesh_summary()`](https://github.com/gcol33/tulpaMesh/reference/mesh_summary.md), [`plot.tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/plot.tulpa_mesh.md) |
| `R/barrier.R` | [`barrier_triangles()`](https://github.com/gcol33/tulpaMesh/reference/barrier_triangles.md), point-in-polygon helper |
| `R/convert.R` | [`as_tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/as_tulpa_mesh.md) for fmesher/INLA mesh conversion |
| `R/refine.R` | [`refine_mesh()`](https://github.com/gcol33/tulpaMesh/reference/refine_mesh.md) adaptive refinement from error indicators |

### C++ layer (`src/mesh.cpp`)

Six `// [[Rcpp::export]]` functions: - `cpp_triangulate()` — CDT
constrained Delaunay, erases outer/holes - `cpp_ruppert_refine()` —
Ruppert algorithm with min angle + max area + segment encroachment -
`cpp_fem_matrices()` — 2D consistent mass + gradient stiffness (triplet
format) - `cpp_fem_matrices_3d()` — 3D surface FEM via metric tensor on
tangent plane - `cpp_projection_matrix()` — grid-accelerated barycentric
(O(n_obs × avg_cell) vs O(n_obs × n_tri)) - `cpp_projection_matrix_3d()`
— 3D barycentric via Cramer’s rule

Shared helpers: `BaryResult barycentric()`, `TriGrid` uniform spatial
index, `circumcenter()`, `min_angle_from_lengths()`.

### Index convention

R input is 1-based → converted to 0-based for CDT → output converted
back to 1-based. FEM triplets returned 0-based from C++ and converted in
R before
[`Matrix::sparseMatrix()`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html).

### Numerical tolerances in C++

- Degenerate triangle/barycentric: `< 1e-15`
- Point-inside boundary (2D): `>= -1e-8`
- Point-inside boundary (3D): `>= -1e-6` (sphere curvature)
- Zero-weight skip in projection: `< 1e-15`

## Testing

141 tests across 8 test files (testthat edition 3):

| File | What it tests |
|----|----|
| `test-mesh.R` | Basic triangulation, formula interface, FEM properties, print |
| `test-fmesher-comparison.R` | Cross-validation against fmesher (C, G, A, C0, va, ta) |
| `test-diagnostics.R` | mesh_quality, mesh_summary, plot.tulpa_mesh |
| `test-sf-boundary.R` | POLYGON, MULTIPOLYGON, holes, sfg objects |
| `test-convert.R` | as_tulpa_mesh from fm_mesh_2d |
| `test-ruppert.R` | Ruppert refinement quality, Steiner points, max_area, FEM on refined |
| `test-barrier.R` | barrier_triangles, stiffness zeroing, sf barriers |
| `test-sphere.R` | Icosahedral subdivision, radii, Euler char, 3D FEM, sphere area ≈ 4π |
| `test-refine.R` | Adaptive refinement, fraction/threshold, min_area guard |

## Key Design Decisions

- **Ruppert refinement**: available via `min_angle` parameter. Iterative
  circumcenter insertion with segment encroachment handling.
- **Consistent mass matrix**: standard for SPDE. Lumped mass (diagonal
  C0) available via `lumped = TRUE`.
- **Grid-accelerated projection**: uniform grid spatial index replaces
  brute-force. Fallback to full scan for edge cases.
- **Barrier models**:
  [`barrier_triangles()`](https://github.com/gcol33/tulpaMesh/reference/barrier_triangles.md)
  identifies barrier regions; `fem_matrices(barrier=)` zeros stiffness
  contributions.
- **Spherical meshes**: icosahedral subdivision with 3D FEM via metric
  tensor. Vertices on exact sphere surface.
- **Triplet accumulation**: C++ returns (i, j, x) vectors; R assembles
  into `dgCMatrix` with `repr = "C"`.

## Ecosystem

tulpaMesh → tulpa (SPDE Q-builder, Laplace solver) → tulpaOcc (occupancy
models with spatial fields). Positioned as CRAN-publishable alternative
to fmesher/INLA mesh tooling.
