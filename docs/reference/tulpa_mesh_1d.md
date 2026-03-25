# Create a 1D Mesh for Temporal SPDE Models

Generates a 1D mesh (interval partition) for temporal components of
spatio-temporal SPDE models. Returns 1D FEM matrices (tridiagonal mass
and stiffness) that can be combined with 2D spatial FEM matrices via
Kronecker products in tulpa's space-time Q-builder.

## Usage

``` r
tulpa_mesh_1d(knots, boundary = NULL, n_extend = 3L)

# S3 method for class 'tulpa_mesh_1d'
print(x, ...)
```

## Arguments

- knots:

  Numeric vector of mesh knot locations (e.g., time points). Will be
  sorted and deduplicated.

- boundary:

  Two-element numeric vector `c(lower, upper)` defining the mesh domain.
  Defaults to `range(knots)`.

- n_extend:

  Number of extra knots to add beyond each boundary, spaced at the
  median knot interval. Default 3.

- x:

  A `tulpa_mesh_1d` object.

- ...:

  Additional arguments (ignored).

## Value

A `tulpa_mesh_1d` object with components:

- `knots`: sorted numeric vector of mesh knot locations

- `n`: number of knots

- `C`: consistent mass matrix (tridiagonal, n x n sparse)

- `G`: stiffness matrix (tridiagonal, n x n sparse)

- `C0`: lumped (diagonal) mass matrix

The `tulpa_mesh_1d` object `x`, returned invisibly.
