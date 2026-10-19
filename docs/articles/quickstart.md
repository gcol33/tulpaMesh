# Quick Start

## Installation

``` r

install.packages("tulpaMesh")            # CRAN

install.packages("pak")                  # development version
pak::pak("gcol33/tulpaMesh")
```

## Your first mesh

tulpaMesh takes point coordinates and returns a triangulated mesh with
FEM matrices ready for SPDE models.

``` r

set.seed(42)
coords <- cbind(x = runif(100), y = runif(100))
mesh <- tulpa_mesh(coords)
mesh
#> tulpa_mesh:
#>   Vertices:   113 
#>   Triangles:  211 
#>   Edges:      323
```

The mesh extends slightly beyond the convex hull of your points
(controlled by `extend`). Plot it:

``` r

plot(mesh, vertex_col = "steelblue", main = "Basic mesh")
```

![](quickstart_files/figure-html/plot-basic-1.svg)

## Controlling mesh density

Use `max_edge` to add refinement points. The mesh generator places a
hexagonal lattice of points at this spacing, producing near-equilateral
triangles.

``` r

mesh_fine <- tulpa_mesh(coords, max_edge = 0.08)
mesh_fine
#> tulpa_mesh:
#>   Vertices:   301 
#>   Triangles:  511 
#>   Edges:      773
plot(mesh_fine, main = "Refined mesh (max_edge = 0.08)")
```

![](quickstart_files/figure-html/refined-mesh-1.svg)

## Getting FEM matrices

[`fem_matrices()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices.md)
returns the three sparse matrices needed for SPDE models:

- **C**: mass matrix (consistent, symmetric positive definite)

- **G**: stiffness matrix (symmetric, zero row sums)

- **A**: projection matrix mapping mesh vertices to observation
  locations

``` r

fem <- fem_matrices(mesh_fine, obs_coords = coords)
dim(fem$C)
#> [1] 301 301
dim(fem$A)
#> [1] 100 301

# Verify key properties
all(Matrix::diag(fem$C) > 0)        # positive diagonal
#> [1] FALSE
max(abs(Matrix::rowSums(fem$G)))     # row sums ~ 0
#> [1] 2.842171e-14
range(Matrix::rowSums(fem$A))        # row sums = 1
#> [1] 1 1
```

For the SPDE Q-builder, you typically need the lumped (diagonal) mass
matrix:

``` r

fem_l <- fem_matrices(mesh_fine, obs_coords = coords, lumped = TRUE)
Matrix::isDiagonal(fem_l$C0)
#> [1] TRUE
```

## Using a formula interface

If your coordinates live in a data.frame, use a formula:

``` r

df <- data.frame(lon = runif(50), lat = runif(50), y = rnorm(50))
mesh_f <- tulpa_mesh(~ lon + lat, data = df)
mesh_f
#> tulpa_mesh:
#>   Vertices:   60 
#>   Triangles:  108 
#>   Edges:      167
```

## Mesh quality

Check triangle quality with
[`mesh_quality()`](https://github.com/gcol33/tulpaMesh/reference/mesh_quality.md)
and
[`mesh_summary()`](https://github.com/gcol33/tulpaMesh/reference/mesh_summary.md):

``` r

mesh_summary(mesh_fine)
#> Mesh quality summary:
#>   Triangles:      511
#>   Min angle:      0.4 / 40.9 / 60.0 (min / median / max)
#>   Max angle:      60.0 / 78.6 / 179.2
#>   Aspect ratio:   1.00 / 1.48 / 35.34
#>   Area:           0.0000 / 0.0020 / 0.0151
#>   WARNING: 29 triangles with min angle < 10 degrees
#>   WARNING: 17 triangles with min angle < 5 degrees (slivers)
```

Color triangles by minimum angle:

``` r

plot(mesh_fine, color = "quality", main = "Colored by minimum angle")
```

![](quickstart_files/figure-html/quality-plot-1.svg)

## Ruppert refinement

For guaranteed minimum angles, use `min_angle`:

``` r

mesh_r <- tulpa_mesh(coords, min_angle = 25, max_edge = 0.15)
mesh_summary(mesh_r)
#> Mesh quality summary:
#>   Triangles:      759
#>   Min angle:      25.0 / 39.9 / 60.0 (min / median / max)
#>   Max angle:      60.0 / 76.9 / 128.9
#>   Aspect ratio:   1.00 / 1.52 / 2.31
#>   Area:           0.0000 / 0.0009 / 0.0116
```

## Next steps

- [Spatial
  Workflows](https://github.com/gcol33/tulpaMesh/articles/workflows.md)
  – boundary constraints, barrier models, sf integration

- [Spherical and Temporal
  Meshes](https://github.com/gcol33/tulpaMesh/articles/advanced.md) –
  global meshes, space-time, metric graphs
