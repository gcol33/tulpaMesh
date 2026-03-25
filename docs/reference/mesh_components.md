# Identify Disconnected Mesh Components

Finds connected components of a mesh via triangle adjacency (two
triangles are connected if they share an edge).

## Usage

``` r
mesh_components(mesh)
```

## Arguments

- mesh:

  A `tulpa_mesh` object.

## Value

An integer vector of length `n_triangles` giving the component ID (1, 2,
...) for each triangle.
