# Compute FEM Matrices from a Mesh

Computes the finite element mass (C) and stiffness (G) matrices from a
triangular mesh, plus the projection matrix (A) that maps mesh vertices
to observation locations.

## Usage

``` r
fem_matrices(
  mesh,
  obs_coords = NULL,
  barrier = NULL,
  parallel = FALSE,
  lumped = FALSE
)
```

## Arguments

- mesh:

  A `tulpa_mesh` object.

- obs_coords:

  Observation coordinates (N x 2 matrix). If NULL, the projection matrix
  A is the identity (observations at mesh nodes).

- barrier:

  Optional logical vector of length `n_triangles`. Triangles marked
  `TRUE` are treated as physical barriers (coastlines, rivers): their
  stiffness contributions are zeroed so the spatial field cannot smooth
  across them. Based on Bakka et al. (2019).

- parallel:

  Logical. If `TRUE`, uses parallel FEM assembly via RcppParallel
  (thread-local triplet accumulation). Beneficial for meshes with \>50K
  triangles. Default `FALSE`.

- lumped:

  Logical. If `TRUE`, returns a diagonal lumped mass matrix C0 (vertex
  areas) in addition to the consistent mass matrix C. The lumped mass
  inverse is trivial and needed for the SPDE Q-builder. Default `FALSE`.

## Value

A list with sparse matrices (dgCMatrix class from Matrix package):

- `C`: consistent mass matrix (n_vertices x n_vertices)

- `G`: stiffness matrix (n_vertices x n_vertices)

- `A`: projection matrix (n_obs x n_vertices)

- `n_mesh`: number of mesh vertices

- `C0`: (only if `lumped = TRUE`) diagonal lumped mass matrix

- `va`: (only if `lumped = TRUE`) vertex areas (numeric vector)

- `ta`: (only if `lumped = TRUE`) triangle areas (numeric vector)
