# tulpaMesh Development Plan

**Date:** 2026-03-25 **Status:** Scaffolded, compiles, 22 tests passing
**Role:** Hard dependency of tulpa for SPDE spatial fields. Also usable
standalone for any spatial triangulation task.

------------------------------------------------------------------------

## What Exists Now

- CDT header-only C++ library vendored in `src/cdt/` (MIT, Artem
  Amirkhanov)
- `cpp_triangulate()`: constrained Delaunay triangulation from points +
  optional edge constraints
- `cpp_fem_matrices()`: linear P1 FEM assembly (mass C, stiffness G)
  from mesh triangles
- `cpp_projection_matrix()`: barycentric interpolation matrix A mapping
  mesh nodes → observation locations
- [`tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh.md):
  R API with formula interface, convex hull extension, boundary
  constraints, max_edge refinement, cutoff deduplication
- [`fem_matrices()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices.md):
  R wrapper returning sparse `dgCMatrix` objects (via Matrix package)
- 22 tests passing

## Architecture

    User coordinates + boundary
            ↓
    tulpaMesh::tulpa_mesh()     ← R API
            ↓ calls
    cpp_triangulate()           ← CDT C++ (constrained Delaunay)
            ↓ returns
    tulpa_mesh object           ← vertices, triangles, edges
            ↓
    tulpaMesh::fem_matrices()   ← R API
            ↓ calls
    cpp_fem_matrices()          ← C++ FEM assembly
    cpp_projection_matrix()     ← C++ barycentric projection
            ↓ returns
    Sparse C, G, A matrices    ← dgCMatrix (Matrix package)
            ↓
    tulpa SPDE Q-builder        ← lives in tulpa, not here

------------------------------------------------------------------------

## Remaining Work

### P0: Must-do before tulpa integration

#### 1. Mesh quality controls

**Files:** `src/mesh.cpp`, `R/mesh.R` **What:** The current mesh has no
refinement quality guarantees. Need: - Minimum angle constraint
(Delaunay refinement à la Ruppert/Chew). CDT doesn’t do this natively —
it’s a constrained Delaunay triangulator, not a mesh generator with
quality guarantees. Options: - Post-process: insert Steiner points at
circumcenters of bad triangles (angle \< threshold), re-triangulate.
This is the Ruppert algorithm and can be implemented on top of CDT. -
Or: accept that for SPDE the mesh quality from CDT + max_edge refinement
is good enough (it usually is — INLA’s fmesher doesn’t guarantee minimum
angles either, it just does CDT + refinement). - **Decision needed:**
Implement Ruppert refinement or accept CDT + max_edge as sufficient? For
v0.1 I’d say accept CDT quality and add Ruppert later if users hit
issues.

#### 2. Mesh quality diagnostics

**Files:** new `R/diagnostics.R` **What:** Functions to assess mesh
quality: -
[`mesh_quality()`](https://github.com/gcol33/tulpaMesh/reference/mesh_quality.md):
returns per-triangle min angle, max angle, aspect ratio, area -
[`mesh_summary()`](https://github.com/gcol33/tulpaMesh/reference/mesh_summary.md):
prints min/median/max of quality metrics -
[`plot.tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/plot.tulpa_mesh.md):
basic mesh visualization (vertices + edges + triangles) - These help
users decide if their mesh needs refinement

#### 3. Boundary handling from sf polygons

**Files:** `R/mesh.R` (the `sf_to_boundary()` internal already exists)
**What:** Currently `sf_to_boundary()` is minimal — it extracts
coordinates from an sf polygon. Needs: - Handle MULTIPOLYGON (multiple
connected components) - Handle holes (inner boundaries as additional
constraint edges) - Test with real sf objects (country boundaries, study
areas) - Add `sf` to Suggests in DESCRIPTION (already there)

#### 4. Matrix package dependency

**Files:** DESCRIPTION **What:**
[`fem_matrices()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices.md)
uses
[`Matrix::sparseMatrix()`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html)
but Matrix isn’t in Imports. Need to add it:

    Imports:
        Rcpp (>= 1.0.12),
        Matrix

Currently works because Matrix is loaded by other packages, but will
fail in a clean session.

#### 5. Wire into tulpa as hard dependency

**Files:** tulpa’s DESCRIPTION, tulpa’s `R/spatial.R` **What:** - Add
`tulpaMesh` to tulpa’s Imports - Create `spatial_spde()` in tulpa
that: 1. Takes coords (or formula + data) + optional boundary 2. Calls
[`tulpaMesh::tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh.md)
to build mesh 3. Calls
[`tulpaMesh::fem_matrices()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices.md)
to get C, G, A 4. Stores these for the SPDE Q-builder - Create the SPDE
Laplace mode function in tulpa’s `laplace_core.cpp`: - Latent field w
lives on mesh nodes (n_mesh elements) - η = Xβ + A·w (projection from
mesh to observations) - Q(κ) = κ⁴C + 2κ²G + GC⁻¹G (classical SPDE
precision) - Priors on w: -0.5 w’Qw - This is a new variant of
`laplace_newton_solve` where the scatter callback maps through A

### P1: Classical SPDE in tulpa

#### 6. Classical SPDE Q-builder in tulpa

**Files:** tulpa `src/spde_qbuilder.h`, `src/spde_qbuilder.cpp`
**What:** Given sparse C, G matrices and hyperparameters (range ρ,
variance σ²): - κ = √(8·ν) / ρ (for ν = 1, α = 2 in 2D) - τ² = 1 / (4π ·
κ² · σ²) - Q = τ² · (κ⁴·C + 2κ²·G + G·C⁻¹·G) - The C⁻¹ is diagonal
(lumped mass matrix) so G·C⁻¹·G is sparse - **Key property:** sparsity
pattern of Q is fixed (same as C + G pattern). Only values change with
θ. This means CHOLMOD `analyze()` once, `factorize()` per grid point. -
Store C, G as `cholmod_sparse*` for direct CHOLMOD operations -
Implement as
`void build_spde_precision(double kappa, double tau, cholmod_sparse* C, cholmod_sparse* G, cholmod_sparse* Q)`
— fills Q values in-place

#### 7. SPDE Laplace mode function in tulpa

**Files:** tulpa `src/laplace_core.cpp` **What:** New
`laplace_mode_spde()` variant:

``` cpp
LaplaceResult laplace_mode_spde(
    y, n_trials, X,
    A_sparse,          // projection matrix (n_obs × n_mesh), cholmod_sparse
    C_sparse, G_sparse, // FEM matrices, cholmod_sparse
    kappa, tau,         // SPDE hyperparameters
    family, phi, max_iter, tol, n_threads,
    x_init, shared_solver
);
```

- Latent field dimension = n_mesh (not n_obs)
- η_i = X_i·β + Σ_j A\_{ij}·w_j
- Scatter: gradient and Hessian contributions mapped through A’
  (transpose of projection)
- Prior: -0.5·w’·Q(κ,τ)·w, gradient = -Q·w, Hessian = Q
- The Hessian H = A’·diag(h_i)·A + Q — this is sparse!
- CHOLMOD handles the full sparse solve

#### 8. Nested Laplace for SPDE (2D grid over range and variance)

**Files:** tulpa `src/laplace_core.cpp`, `R/nested_laplace.R` **What:**
`cpp_nested_laplace_spde()` — C++ grid loop over (ρ, σ²) pairs: - CCD
grid with k=2 hyperparameters → 9 grid points - At each point: rebuild Q
from (ρ, σ²), run inner Laplace, warm-start from previous - CHOLMOD:
analyze once (Q sparsity is the same for all (ρ, σ²)), factorize per
grid point - This is the “INLA replacement” for continuous spatial
fields

#### 9. R-level SPDE API in tulpa

**Files:** tulpa `R/spatial.R` **What:**

``` r

spatial_spde <- function(coords, data = NULL, mesh = NULL,
                         boundary = NULL, max_edge = NULL,
                         nu = 1, prior_range = NULL, prior_sigma = NULL) {
  # If mesh not provided, build one
  if (is.null(mesh)) {
    mesh <- tulpaMesh::tulpa_mesh(coords, data = data,
                                   boundary = boundary, max_edge = max_edge)
  }
  # Compute FEM matrices
  obs_coords <- resolve_coords(coords, data)
  fem <- tulpaMesh::fem_matrices(mesh, obs_coords = obs_coords)

  structure(list(
    type = "spde", mesh = mesh,
    C = fem$C, G = fem$G, A = fem$A,
    n_mesh = fem$n_mesh, nu = nu,
    prior_range = prior_range, prior_sigma = prior_sigma
  ), class = "tulpa_spatial")
}
```

### P2: Rational SPDE (the differentiator)

#### 10. Rational SPDE Q-builder

**Files:** tulpa `src/spde_qbuilder.h` **What:** Bolin, Simas & Xiong
(2023) rational approximation: - For fractional α = ν + d/2 (where d=2
in 2D), the precision operator (κ² - Δ)^α can’t be represented exactly
with the classical SPDE stencil - Rational approximation: P_m(L) ≈ L^α
where L = κ²C + G and P_m is a degree-m rational polynomial - Q = P_m’ ·
C · P_m - The polynomial coefficients come from a best rational
approximation (Padé or Chebyshev) - **Advantage:** arbitrary ν \> 0, not
just ν ∈ {1, 2}. Users can estimate ν from data. - **Implementation:**
Pre-compute rational polynomial coefficients for given m (typically m=2
or m=3 is sufficient). Build Q as a product of sparse matrices. -
Reference: `rSPDE` package by David Bolin (R implementation exists, we’d
reimplement the Q-builder in C++ for CHOLMOD integration)

#### 11. Smoothness parameter estimation

**Files:** tulpa `R/nested_laplace.R` **What:** With rational SPDE, ν
becomes a hyperparameter that can be estimated: - 3D grid over (range,
variance, smoothness) — CCD with k=3 → 15 grid points - Or: fix ν and do
2D grid, then profile over ν - This is what INLA can’t do natively —
tulpa would be the first CRAN package to offer nested Laplace with
estimated Matérn smoothness

### P3: Advanced mesh features (tulpaMesh)

#### 12. Mesh on spheres

**What:** For global-scale species distribution models, need spherical
meshes. CDT works on 2D planes; for spheres, either: - Project to plane
(stereographic), triangulate, project back — works for small regions -
Use icosahedral subdivision — standard for global meshes - Or use
STRIPACK/similar spherical Delaunay algorithm

#### 13. Mesh refinement (Ruppert algorithm)

**What:** Quality-guaranteed mesh generation. Insert circumcenters of
triangles with min angle \< threshold. Converges to a mesh where all
angles ≥ ~20.7°.

#### 14. Barrier models

**What:** Physical barriers (coastlines, rivers) that the spatial field
shouldn’t cross. Implemented by: - Setting the range to near-zero inside
barrier triangles - Or: removing barrier triangles from the FEM
assembly - Reference: Bakka et al. (2019) barrier model

#### 15. Adaptive mesh refinement

**What:** Refine mesh where data density is high, coarsen where sparse.
Uses error indicators from the SPDE solution to guide refinement.
Iterative: solve → refine → re-solve.

------------------------------------------------------------------------

## Dependency Chain

    tulpaMesh (mesh + FEM)
        ↓ hard dep
    tulpa (SPDE Q-builder + CHOLMOD + nested Laplace)
        ↓ hard dep
    tulpaOcc (occupancy models with SPDE spatial fields)

User sees:

``` r

library(tulpaOcc)
fit <- occ(~ elevation, ~ 1, data = survey_data,
           spatial = spatial_spde(coords = ~ x + y, max_edge = c(0.5, 2)))
```

Under the hood: 1. `spatial_spde()` calls
[`tulpaMesh::tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh.md)
→ mesh 2. `spatial_spde()` calls
[`tulpaMesh::fem_matrices()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices.md)
→ C, G, A 3. `occ()` calls tulpa’s `cpp_laplace_fit_spde()` or
`cpp_nested_laplace_spde()` → fit 4. CHOLMOD handles all sparse solves

------------------------------------------------------------------------

## Competitive Position

| Feature | INLA | tulpa + tulpaMesh | spOccupancy | TMB |
|----|----|----|----|----|
| Mesh generation | fmesher (CRAN) | tulpaMesh (CRAN) | No | No |
| Classical SPDE | Yes | P1 (planned) | No | Manual |
| Rational SPDE (fractional ν) | No | P2 (planned) | No | Via rSPDE |
| Estimated smoothness | No (fixed ν) | P2 (planned) | No | Via tmbstan |
| Nested Laplace | Yes | Yes (built) | No | No |
| CHOLMOD solver | PARDISO (proprietary) | CHOLMOD (MIT, ships with R) | No | Eigen |
| CRAN | No | Yes | Yes | Yes |
| Occupancy models | Via inlaOcc | tulpaOcc (native) | Native | Manual |
| GPU | No | Planned (NNGP batching) | No | No |

**Key differentiator:** tulpa would be the only CRAN-available engine
offering rational SPDE + nested Laplace + estimated Matérn smoothness +
occupancy models in a single install.

------------------------------------------------------------------------

## Open Questions

1.  **CDT license in CRAN context:** CDT is MPL-2.0. The vendored
    headers are used only at compile time (header-only). Need to verify
    MPL-2.0 is compatible with MIT (tulpaMesh) for CRAN. MPL-2.0 allows
    “larger work” under different license as long as MPL files retain
    their license. Should be fine — document in `inst/COPYRIGHTS` or
    `LICENSE.note`.

2.  **Consistent vs lumped mass matrix:** Currently using consistent
    mass (area/6 diagonal, area/12 off-diagonal). The SPDE Q-builder
    needs C⁻¹ — lumped mass (diagonal) makes this trivial. Consistent
    mass is more accurate but requires solving a sparse system for
    C⁻¹·G. Decision: use lumped mass for the Q-builder (standard
    practice in INLA), keep consistent mass available for other uses.

3.  **Projection matrix point location:** Current
    `cpp_projection_matrix` does brute-force triangle search (O(n_obs ×
    n_tri)). For large meshes, need a spatial index (KD-tree or triangle
    walk). CDT’s KDTree.h is already vendored — use it for point
    location. Not blocking for v0.1 but needed for n_obs \> 10K.

4.  **fmesher compatibility:** Should tulpaMesh be able to import
    fmesher meshes? Would let users switch from INLA without re-meshing.
    Low priority but good for adoption. A simple
    [`as_tulpa_mesh.fm_mesh_2d()`](https://github.com/gcol33/tulpaMesh/reference/as_tulpa_mesh.md)
    converter would suffice.

5.  **Matrix package version:**
    [`fem_matrices()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices.md)
    uses
    [`Matrix::sparseMatrix()`](https://rdrr.io/pkg/Matrix/man/sparseMatrix.html).
    Need `Matrix (>= 1.5-0)` for the `repr` argument. Add version
    constraint to Imports.

------------------------------------------------------------------------

## Session Log (2026-03-25)

- Scaffolded tulpaMesh package with CDT vendored headers (11 files)
- Implemented `cpp_triangulate()`, `cpp_fem_matrices()`,
  `cpp_projection_matrix()` in C++
- Implemented
  [`tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh.md)
  and
  [`fem_matrices()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices.md)
  R API
- 22 tests passing
- Smoke test: 100 points → 113 vertices, 211 triangles, FEM matrices 94%
  sparse
- Fixed CDT API (`V2d(x,y)` not `V2d::make(x,y)`)

### Also done in tulpa this session:

- Feature 9: Consolidated Laplace code (3100 → 900 lines)
- Feature 1: CHOLMOD sparse solver via Matrix package
- Feature 5: Hot-start (warm-starting across grid points)
- Feature 3a: R-level nested Laplace prototype
- Feature 3b: C++ nested Laplace grid loop with shared CHOLMOD +
  warm-start chain
- Fixed 3 pre-existing hmc_sampler.cpp build errors
- Total: 225 tests passing in tulpa
