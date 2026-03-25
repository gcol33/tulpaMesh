# Compute P2 (Quadratic) FEM Matrices

Computes finite element matrices using 6-node quadratic triangular
elements. Each triangle edge gets a midpoint node, giving 6 basis
functions per triangle instead of 3. Quadratic elements provide better
approximation accuracy with fewer mesh nodes.

## Usage

``` r
fem_matrices_p2(mesh)
```

## Arguments

- mesh:

  A `tulpa_mesh` object (2D only).

## Value

A list with:

- `C`: consistent mass matrix (n_total x n_total sparse)

- `G`: stiffness matrix (n_total x n_total sparse)

- `n_mesh`: total number of nodes (vertices + midpoints)

- `n_vertices`: number of original mesh vertices

- `n_midpoints`: number of added midpoint nodes

- `vertices`: N x 2 matrix of all node coordinates

- `triangles6`: M x 6 connectivity matrix (columns: v0, v1, v2, mid01,
  mid12, mid20)
