# Spherical and Temporal Meshes

## Spherical meshes

For global-scale species distribution models or climate fields, you need
a mesh on the sphere.
[`tulpa_mesh_sphere()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh_sphere.md)
generates one by recursive subdivision of an icosahedron.

``` r

globe <- tulpa_mesh_sphere(subdivisions = 3)
globe
#> tulpa_mesh (sphere, radius = 1 ):
#>   Vertices:   642 
#>   Triangles:  1280 
#>   Edges:      1920
```

Each subdivision level quadruples the triangle count:

| Level | Vertices | Triangles |
|-------|----------|-----------|
| 0     | 12       | 20        |
| 1     | 42       | 80        |
| 2     | 162      | 320       |
| 3     | 642      | 1280      |
| 4     | 2562     | 5120      |

All vertices lie exactly on the sphere surface:

``` r

radii <- sqrt(rowSums(globe$vertices^2))
range(radii)
#> [1] 1 1
```

### FEM on the sphere

The 3D FEM assembly uses the metric tensor on each triangle’s tangent
plane:

``` r

fem <- fem_matrices(globe, lumped = TRUE)

# Total surface area should approximate 4*pi for unit sphere
cat("Total area:", sum(fem$va), "\n")
#> Total area: 12.50649
cat("4*pi:      ", 4 * pi, "\n")
#> 4*pi:       12.56637
cat("Error:     ", abs(sum(fem$va) - 4 * pi) / (4 * pi) * 100, "%\n")
#> Error:      0.476493 %
```

### Projection with lon/lat coordinates

Pass observation coordinates as (longitude, latitude) in degrees:

``` r

obs <- cbind(lon = c(0, 90, -45, 170), lat = c(0, 45, -30, 60))
fem <- fem_matrices(globe, obs_coords = obs)
dim(fem$A)
#> [1]   4 642
range(Matrix::rowSums(fem$A))  # row sums = 1
#> [1] 1 1
```

### Earth-radius meshes

``` r

earth <- tulpa_mesh_sphere(subdivisions = 3, radius = 6371)
sqrt(sum(earth$vertices[1, ]^2))  # 6371 km
#> [1] 6371
```

### Euler characteristic

A closed sphere satisfies V - E + T = 2:

``` r

globe$n_vertices - globe$n_edges + globe$n_triangles
#> [1] 2
```

## 1D temporal meshes

For spatio-temporal SPDE models, the temporal dimension needs its own
mesh.
[`tulpa_mesh_1d()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh_1d.md)
produces knots with tridiagonal FEM matrices that combine with the
spatial mesh via Kronecker products.

``` r

# Monthly observations over 5 years
times <- seq(2020, 2025, by = 1/12)
m1d <- tulpa_mesh_1d(times)
m1d
#> tulpa_mesh_1d:
#>   Knots:     67 
#>   Range:    [2019.75, 2025.25]
#>   Intervals: 66
```

The extension knots (controlled by `n_extend`) prevent boundary effects:

``` r

range(times)       # data range
#> [1] 2020 2025
range(m1d$knots)   # mesh range (extended)
#> [1] 2019.75 2025.25
```

### 1D FEM properties

Mass and stiffness matrices are tridiagonal and sparse:

``` r

# No extension for cleaner demonstration
m <- tulpa_mesh_1d(seq(0, 1, by = 0.1), n_extend = 0)

# Symmetric, positive definite mass
Matrix::isSymmetric(m$C)
#> [1] TRUE
all(Matrix::diag(m$C) > 0)
#> [1] TRUE

# Stiffness row sums = 0
max(abs(Matrix::rowSums(m$G)))
#> [1] 1.776357e-15

# Total mass = domain length
sum(m$C)
#> [1] 1
```

### Irregular time spacing

The 1D mesh handles irregular intervals naturally:

``` r

# Denser sampling in summer, sparser in winter
times_irr <- c(
  seq(2020, 2020.25, by = 1/52),      # weekly Jan-Mar
  seq(2020.25, 2020.75, by = 1/365),   # daily Apr-Sep
  seq(2020.75, 2021, by = 1/52)        # weekly Oct-Dec
)
m_irr <- tulpa_mesh_1d(times_irr, n_extend = 0)
cat(m_irr$n, "knots from", length(times_irr), "unique time points\n")
#> 210 knots from 211 unique time points
```

## Metric graph meshes

For SPDE models along roads, rivers, or coastlines,
[`tulpa_mesh_graph()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh_graph.md)
builds a 1D FEM mesh along network edges.

``` r

# Simple river network: main channel + tributary
edges <- list(
  cbind(x = seq(0, 10, by = 0.5), y = rep(0, 21)),     # main channel
  cbind(x = c(5, 5, 5, 5), y = c(0, 2, 4, 6))          # tributary
)

g <- tulpa_mesh_graph(edges)
g
#> tulpa_mesh_graph:
#>   Vertices:  24 
#>   Segments:  23 
#>   Junctions: 1   Endpoints: 3
```

### Junction detection

Nodes where edges meet are automatically detected:

``` r

cat("Junctions (degree > 2):", sum(g$degree > 2), "\n")
#> Junctions (degree > 2): 1
cat("Endpoints (degree = 1):", sum(g$degree == 1), "\n")
#> Endpoints (degree = 1): 3
```

### FEM on graphs

The graph FEM matrices are structurally identical to the 1D case but
assembled across the full network:

``` r

Matrix::isSymmetric(g$C)
#> [1] TRUE
max(abs(Matrix::rowSums(g$G)))  # row sums ~ 0
#> [1] 0
```

### Subdividing long edges

``` r

g_fine <- tulpa_mesh_graph(edges, max_edge = 0.3)
cat("Coarse:", g$n_vertices, "vertices\n")
#> Coarse: 24 vertices
cat("Fine:  ", g_fine$n_vertices, "vertices\n")
#> Fine:   62 vertices
```

### sf LINESTRING input

``` r

library(sf)
#> Linking to GEOS 3.13.1, GDAL 3.11.4, PROJ 9.7.0; sf_use_s2() is TRUE
line1 <- st_linestring(cbind(c(0, 5, 10), c(0, 3, 0)))
line2 <- st_linestring(cbind(c(5, 5), c(3, 8)))
g_sf <- tulpa_mesh_graph(st_sfc(line1, line2), max_edge = 1)
g_sf
#> tulpa_mesh_graph:
#>   Vertices:  18 
#>   Segments:  17 
#>   Junctions: 1   Endpoints: 3
```

## Non-stationary FEM

For fields where the range or variance changes across space,
[`fem_matrices_nonstationary()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices_nonstationary.md)
computes weighted FEM matrices in C++:

``` r

set.seed(42)
mesh <- tulpa_mesh(cbind(runif(50), runif(50)), max_edge = 0.15)
n <- mesh$n_vertices

# Range decreases from left to right
kappa <- sqrt(8) / (2 - mesh$vertices[, 1])  # shorter range on the right
tau <- rep(1, n)

ns <- fem_matrices_nonstationary(mesh, kappa, tau)
names(ns)
#> [1] "Ck"     "Gk"     "Ct"     "C"      "G"      "C0"     "n_mesh"
```

With constant kappa/tau, the weighted matrices are simply scaled
versions of the standard ones:

``` r

ns_const <- fem_matrices_nonstationary(mesh, rep(2, n), rep(3, n))
fem <- fem_matrices(mesh)

max(abs(ns_const$Ck - 4 * fem$C))  # kappa^2 = 4
#> [1] 0
max(abs(ns_const$Ct - 9 * fem$C))  # tau^2 = 9
#> [1] 2.775558e-17
```

## P2 quadratic elements

For higher accuracy with fewer mesh nodes, use 6-node quadratic
elements:

``` r

set.seed(42)
mesh <- tulpa_mesh(cbind(runif(30), runif(30)))
p2 <- fem_matrices_p2(mesh)

cat("P1 nodes:", mesh$n_vertices, "\n")
#> P1 nodes: 40
cat("P2 nodes:", p2$n_mesh, "(", p2$n_vertices, "vertices +",
    p2$n_midpoints, "midpoints)\n")
#> P2 nodes: 147 ( 40 vertices + 107 midpoints)
```

The total area is preserved:

``` r

fem_p1 <- fem_matrices(mesh)
cat("P1 total area:", sum(fem_p1$C), "\n")
#> P1 total area: 0.9660631
cat("P2 total area:", sum(p2$C), "\n")
#> P2 total area: 0.9660631
```

## Parallel FEM assembly

For large meshes (\>50K triangles), enable parallel assembly:

``` r

set.seed(42)
mesh_large <- tulpa_mesh(cbind(runif(500), runif(500)), max_edge = 0.03)
cat(mesh_large$n_triangles, "triangles\n")
#> 4014 triangles

fem_seq <- fem_matrices(mesh_large, parallel = FALSE)
fem_par <- fem_matrices(mesh_large, parallel = TRUE)

# Results are identical
max(abs(fem_seq$C - fem_par$C))
#> [1] 2.168404e-19
```
