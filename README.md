<!-- badges: start -->
[![R-CMD-check](https://github.com/gcol33/tulpaMesh/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/gcol33/tulpaMesh/actions/workflows/R-CMD-check.yml)
[![Codecov test coverage](https://codecov.io/gh/gcol33/tulpaMesh/graph/badge.svg)](https://app.codecov.io/gh/gcol33/tulpaMesh)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

# tulpaMesh

**Spatial meshes and FEM matrices for SPDE models, from point clouds to precision matrices in three lines.**

## Quick Start

```r
library(tulpaMesh)

coords <- cbind(runif(200), runif(200))
mesh <- tulpa_mesh(coords, max_edge = 0.1)
fem  <- fem_matrices(mesh, obs_coords = coords)
# fem$C (mass), fem$G (stiffness), fem$A (projection) ready for SPDE
```

## Statement of Need

Spatial models based on the SPDE approach (Lindgren, Rue & Lindstrom 2011) need a triangulated mesh and finite element matrices. tulpaMesh handles mesh generation, quality refinement, and FEM matrix assembly with no external dependencies beyond Rcpp and Matrix.

Common applications:

- Species distribution models with spatial random effects
- Geostatistical models with Matern covariance via SPDE
- Occupancy models with spatially structured detection
- Any Bayesian hierarchical model needing a spatial mesh

## Features

### Mesh Generation

- `tulpa_mesh()` - Constrained Delaunay from points, formulas, or sf boundaries
- `tulpa_mesh_sphere()` - Icosahedral geodesic meshes for global-scale models
- `tulpa_mesh_1d()` - Temporal meshes for space-time SPDE
- `tulpa_mesh_graph()` - Network meshes for river/road SPDE

### FEM Assembly

- `fem_matrices()` - Mass (C), stiffness (G), and projection (A) matrices
- `fem_matrices(lumped = TRUE)` - Diagonal mass for SPDE Q-builder
- `fem_matrices(parallel = TRUE)` - Multithreaded assembly via TBB
- `fem_matrices(barrier = ...)` - Barrier models (Bakka et al. 2019)
- `fem_matrices_p2()` - Quadratic (P2) elements for higher accuracy
- `fem_matrices_nonstationary()` - Spatially varying range and variance

### Quality & Diagnostics

- `mesh_quality()` - Per-triangle angles, aspect ratio, area
- `mesh_summary()` - Min/median/max quality statistics
- `plot()` - Mesh visualization with optional quality coloring
- Ruppert refinement via `min_angle` and `max_area` parameters

### Mesh Operations

- `subdivide_mesh()` - Uniform 1-to-4 triangle subdivision
- `subset_mesh()` - Extract submesh from triangle indices
- `mesh_components()` - Find disconnected components
- `refine_mesh()` - Adaptive refinement from error indicators
- `as_tulpa_mesh()` - Convert fmesher/INLA meshes
- `barrier_triangles()` - Identify barrier regions from polygons

## Installation

Install from CRAN:

```r
install.packages("tulpaMesh")
```

Development version:

```r
# install.packages("pak")
pak::pak("gcol33/tulpaMesh")
```

## Usage Examples

### Mesh with Ruppert Refinement

```r
library(tulpaMesh)

set.seed(42)
coords <- cbind(runif(100), runif(100))
mesh <- tulpa_mesh(coords, min_angle = 25, max_edge = 0.15)
mesh_summary(mesh)
plot(mesh, color = "quality", vertex_col = "black")
```

### Barrier Model (Coastline)

```r
library(sf)

# Barrier polygon (e.g., land mass that field can't cross)
barrier <- st_polygon(list(rbind(c(3,3), c(7,3), c(7,7), c(3,7), c(3,3))))
bt <- barrier_triangles(mesh, st_sfc(barrier))
fem <- fem_matrices(mesh, barrier = bt)
# fem$G has zero stiffness across barrier triangles
```

### Spherical Mesh for Global Models

```r
globe <- tulpa_mesh_sphere(subdivisions = 4)
globe
#> tulpa_mesh (sphere, radius = 1):
#>   Vertices:   2562
#>   Triangles:  5120
#>   Edges:      7680
fem <- fem_matrices(globe, lumped = TRUE)
# Total surface area approximates 4*pi
sum(fem$va)
```

### Migrate from fmesher

```r
library(fmesher)
fm <- fm_mesh_2d(loc = coords, max.edge = c(0.3, 0.6))
tm <- as_tulpa_mesh(fm)  # direct conversion
fem <- fem_matrices(tm)  # same C, G, A as fmesher
```

## Documentation

- [Quick Start](https://gillescolling.com/tulpaMesh/articles/quickstart.html)
- [Spatial Workflows](https://gillescolling.com/tulpaMesh/articles/workflows.html) -- boundaries, barriers, sf integration
- [Spherical and Temporal Meshes](https://gillescolling.com/tulpaMesh/articles/advanced.html) -- globe, space-time, graphs, P2
- [Full Reference](https://gillescolling.com/tulpaMesh/reference/)

## Support

> "Software is like sex: it's better when it's free." -- Linus Torvalds

I'm a PhD student who builds R packages in my free time because I believe good tools should be free and open. I started these projects for my own work and figured others might find them useful too.

If this package saved you some time, buying me a coffee is a nice way to say thanks. It helps with my coffee addiction.

[![Buy Me A Coffee](https://img.shields.io/badge/-Buy%20me%20a%20coffee-FFDD00?logo=buymeacoffee&logoColor=black)](https://buymeacoffee.com/gcol33)

## License

MIT (see the LICENSE.md file)

## Citation

```bibtex
@software{tulpaMesh,
  author = {Colling, Gilles},
  title = {tulpaMesh: Constrained Delaunay Triangulation Meshes for Spatial SPDE Models},
  year = {2026},
  url = {https://CRAN.R-project.org/package=tulpaMesh},
  doi = {10.32614/CRAN.package.tulpaMesh}
}
```
