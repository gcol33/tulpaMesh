# Identify Barrier Triangles

Determines which mesh triangles fall inside barrier regions (e.g.,
coastlines, rivers, lakes). A triangle is marked as a barrier if its
centroid falls inside any barrier polygon.

## Usage

``` r
barrier_triangles(mesh, barriers)
```

## Arguments

- mesh:

  A `tulpa_mesh` object.

- barriers:

  An sf/sfc object defining barrier regions, or a list of N x 2
  coordinate matrices (each defining a closed polygon).

## Value

A logical vector of length `n_triangles`. `TRUE` indicates the triangle
is inside a barrier region.
