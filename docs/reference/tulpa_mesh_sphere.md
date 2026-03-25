# Create a Triangular Mesh on the Sphere

Generates a geodesic mesh on the unit sphere by recursive subdivision of
an icosahedron. Optionally refines around data locations using
stereographic projection and local Delaunay triangulation.

## Usage

``` r
tulpa_mesh_sphere(subdivisions = 3L, coords = NULL, radius = 1)
```

## Arguments

- subdivisions:

  Number of recursive icosahedral subdivisions. Each level quadruples
  the triangle count: 0 = 20 triangles, 1 = 80, 2 = 320, 3 = 1280, 4 =
  5120, 5 = 20480. Default 3.

- coords:

  Optional matrix of lon/lat coordinates (degrees) to insert into the
  mesh. Points are projected to the sphere and added as additional
  vertices via local re-triangulation.

- radius:

  Sphere radius. Default 1 (unit sphere). For Earth, use 6371 (km).

## Value

A `tulpa_mesh` object with components:

- `vertices`: N x 3 matrix of (x, y, z) Cartesian coordinates on the
  sphere surface

- `triangles`: M x 3 integer matrix of vertex indices (1-based)

- `edges`: K x 2 integer matrix of edge vertex indices (1-based)

- `lonlat`: N x 2 matrix of (longitude, latitude) in degrees

- `n_vertices`, `n_triangles`, `n_edges`: counts

- `n_input_points`: number of original input points

- `sphere`: list with `radius` and `subdivisions`
