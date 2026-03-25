# Create a Metric Graph Mesh

Builds a 1D FEM mesh along the edges of a spatial network (roads,
rivers, coastlines). Each edge is discretized into segments with P1
linear elements. Junction nodes where edges meet are shared. Based on
the Whittle-Matern SPDE formulation on metric graphs (Bolin, Simas &
Wallin, 2024).

## Usage

``` r
tulpa_mesh_graph(edges, max_edge = NULL, snap_tolerance = 1e-08)

# S3 method for class 'tulpa_mesh_graph'
print(x, ...)
```

## Arguments

- edges:

  A list of N x 2 coordinate matrices, each defining a polyline edge of
  the network. Or an sf object with LINESTRING geometries.

- max_edge:

  Maximum segment length along edges. Edges longer than this are
  subdivided. Default `NULL` (use vertices as-is).

- snap_tolerance:

  Distance below which endpoints are snapped to the same junction node.
  Default `1e-8`.

- x:

  A `tulpa_mesh_graph` object.

- ...:

  Additional arguments (ignored).

## Value

A `tulpa_mesh_graph` object with components:

- `vertices`: N x 2 matrix of node coordinates

- `segments`: M x 2 integer matrix of segment connectivity (1-based)

- `n_vertices`, `n_segments`: counts

- `C`: consistent mass matrix (tridiagonal-block sparse)

- `G`: stiffness matrix (tridiagonal-block sparse)

- `C0`: lumped (diagonal) mass matrix

- `degree`: integer vector of node degrees (junction = degree \> 2)

The `tulpa_mesh_graph` object `x`, returned invisibly.
