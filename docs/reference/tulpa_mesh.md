# Create a Triangular Mesh for SPDE Spatial Models

Generates a constrained Delaunay triangulation from point coordinates,
optionally with boundary constraints. The resulting mesh can be used
directly with tulpa's SPDE spatial fields.

## Usage

``` r
tulpa_mesh(
  coords,
  data = NULL,
  boundary = NULL,
  max_edge = NULL,
  cutoff = 0,
  extend = 0.1,
  min_angle = NULL,
  max_area = NULL,
  max_steiner = 10000L
)

# S3 method for class 'tulpa_mesh'
print(x, ...)
```

## Arguments

- coords:

  A matrix or data.frame with columns x and y, or a formula like
  `~ x + y` evaluated in `data`.

- data:

  Optional data.frame for formula evaluation.

- boundary:

  Optional boundary specification: a matrix of boundary vertex
  coordinates (N x 2), an sf polygon, or NULL for convex hull.

- max_edge:

  Maximum edge length. A single value or a vector of two values
  `c(inner, outer)` where `inner` controls the study area and `outer`
  controls the extension region.

- cutoff:

  Minimum distance between mesh vertices. Points closer than this are
  merged. Default 0 (no merging).

- extend:

  Numeric extension factor beyond the boundary. Default 0.1 (10% of
  domain diameter). Set to 0 for no extension.

- min_angle:

  Minimum angle (degrees) for Ruppert refinement. If specified, Steiner
  points are inserted at circumcenters of triangles with angles below
  this threshold. Theoretical maximum is ~20.7 degrees; values up to 30
  usually work. Default `NULL` (no refinement).

- max_area:

  Maximum triangle area for refinement. Triangles larger than this are
  refined regardless of angle quality. Default `NULL` (no area
  constraint).

- max_steiner:

  Maximum number of Steiner points to insert during Ruppert refinement.
  Default 10000.

- x:

  A `tulpa_mesh` object.

- ...:

  Additional arguments (ignored).

## Value

A `tulpa_mesh` object with components:

- `vertices`: N x 2 matrix of vertex coordinates

- `triangles`: M x 3 integer matrix of vertex indices (1-based)

- `edges`: K x 2 integer matrix of edge vertex indices (1-based)

- `n_vertices`: number of mesh vertices

- `n_triangles`: number of triangles

- `n_edges`: number of edges

- `n_input_points`: number of original input points

The `tulpa_mesh` object `x`, returned invisibly.

## Examples

``` r
# Simple mesh from random points
set.seed(42)
coords <- cbind(runif(50), runif(50))
mesh <- tulpa_mesh(coords)
print(mesh)

# Mesh with Ruppert refinement (min angle 25 degrees)
mesh_refined <- tulpa_mesh(coords, min_angle = 25)
```
