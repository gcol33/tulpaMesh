# Non-stationary FEM Matrices with Spatially Varying Parameters

Computes weighted FEM matrices for non-stationary SPDE models where the
range and variance parameters vary spatially. The weights are per-vertex
kappa(s) and tau(s) values, interpolated within each triangle using the
element-average approximation.

## Usage

``` r
fem_matrices_nonstationary(mesh, kappa, tau)
```

## Arguments

- mesh:

  A `tulpa_mesh` object (2D only).

- kappa:

  Numeric vector of length `n_vertices` giving the spatial scale
  parameter kappa(s) = sqrt(8\*nu) / range(s).

- tau:

  Numeric vector of length `n_vertices` giving the precision scaling
  tau(s).

## Value

A list with sparse matrices:

- `Ck`: kappa²-weighted mass matrix

- `Gk`: kappa²-weighted stiffness matrix

- `Ct`: tau²-weighted mass matrix

- `C`: unweighted consistent mass matrix

- `G`: unweighted stiffness matrix

- `C0`: lumped (diagonal) mass matrix

- `n_mesh`: number of mesh vertices
