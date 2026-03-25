# Subdivide a Mesh

Splits each triangle into 4 sub-triangles by inserting edge midpoints,
then re-triangulates. Useful for multi-resolution workflows where a
coarser mesh needs uniform refinement.

## Usage

``` r
subdivide_mesh(mesh)
```

## Arguments

- mesh:

  A `tulpa_mesh` object (2D only).

## Value

A new `tulpa_mesh` object with approximately 4x as many triangles.
