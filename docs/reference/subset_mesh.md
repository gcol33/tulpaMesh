# Extract a Submesh from Triangle Indices

Creates a new `tulpa_mesh` containing only the specified triangles, with
vertex indices remapped.

## Usage

``` r
subset_mesh(mesh, triangles)
```

## Arguments

- mesh:

  A `tulpa_mesh` object.

- triangles:

  Integer vector of triangle indices to keep, or a logical vector of
  length `n_triangles`.

## Value

A new `tulpa_mesh` object containing only the selected triangles and
their vertices.
