# tulpaMesh

*a mesh from your scattered points*

<!-- badges: start -->
[![R-CMD-check](https://github.com/gcol33/tulpaMesh/actions/workflows/R-CMD-check.yml/badge.svg)](https://github.com/gcol33/tulpaMesh/actions/workflows/R-CMD-check.yml)
[![Codecov test coverage](https://codecov.io/gh/gcol33/tulpaMesh/graph/badge.svg)](https://app.codecov.io/gh/gcol33/tulpaMesh)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
<!-- badges: end -->

**Spatial meshes and finite-element matrices for SPDE models, built on a constrained Delaunay triangulation in C++.**

Hand it a point cloud. `tulpaMesh` triangulates it, refines the triangles to a
guaranteed minimum angle (Ruppert), and assembles the mass, stiffness, and
projection matrices the SPDE approach (Lindgren, Rue & Lindstrom 2011) needs.
The triangulation runs on the [CDT](https://github.com/artem-ogre/CDT)
header-only library with exact geometric predicates; the FEM assembly has no
dependency beyond Rcpp and Matrix. It is the mesh backend for the `tulpa`
Bayesian engine and runs standalone for any spatial triangulation task.

```r
library(tulpaMesh)

coords <- cbind(runif(200), runif(200))
mesh   <- tulpa_mesh(coords, max_edge = 0.1)
fem    <- fem_matrices(mesh, obs_coords = coords)
# fem$C (mass), fem$G (stiffness), fem$A (projection), ready for SPDE
```

## Drop-in for fmesher

If you already build SPDE meshes with `fmesher`, convert an existing mesh and
get the same C, G, and A matrices:

```r
library(fmesher)
fm  <- fm_mesh_2d(loc = coords, max.edge = c(0.3, 0.6))
tm  <- as_tulpa_mesh(fm)        # direct conversion
fem <- fem_matrices(tm)         # same C, G, A
```

## Mesh generation

- `tulpa_mesh()`: constrained Delaunay from points, a `~ x + y` formula, or an
  `sf` boundary (POLYGON, MULTIPOLYGON, holes). Ruppert refinement via
  `min_angle` and `max_area`; hexagonal-lattice density via `max_edge`.
- `tulpa_mesh_sphere()`: icosahedral geodesic meshes with 3D surface FEM, for
  global-scale models.
- `tulpa_mesh_1d()`: temporal meshes with tridiagonal FEM, for space-time SPDE.
- `tulpa_mesh_graph()`: metric-graph meshes for river and road networks.

## FEM assembly

`fem_matrices()` returns the mass (C), stiffness (G), and projection (A) sparse
matrices. The options cover the cases SPDE work runs into:

- `lumped = TRUE`: diagonal mass matrix for the SPDE Q-builder.
- `parallel = TRUE`: multithreaded assembly via RcppParallel/TBB.
- `barrier = ...`: barrier models (Bakka et al. 2019), where the field cannot
  cross a region such as a coastline.
- `fem_matrices_p2()`: quadratic (P2) elements on 6-node triangles for higher
  accuracy.
- `fem_matrices_nonstationary()`: spatially varying range and variance.

## Barrier model on a coastline

```r
library(sf)

barrier <- st_polygon(list(rbind(c(3,3), c(7,3), c(7,7), c(3,7), c(3,3))))
bt      <- barrier_triangles(mesh, st_sfc(barrier))
fem     <- fem_matrices(mesh, barrier = bt)
# fem$G carries zero stiffness across the barrier triangles
```

## Quality, diagnostics, and mesh operations

- `mesh_quality()` / `mesh_summary()`: per-triangle angles, aspect ratio, area,
  and their min/median/max.
- `plot()`: mesh visualization, with `color = "quality"` to shade by minimum
  angle.
- `subdivide_mesh()`, `subset_mesh()`, `mesh_components()`: 1-to-4 subdivision,
  submesh extraction, connected-component finding.
- `refine_mesh()`: adaptive refinement from per-triangle error indicators.
- `barrier_triangles()`, `mesh_crs()`, `set_crs()`: barrier identification and
  CRS handling.

## Installation

```r
install.packages("tulpaMesh")            # CRAN

install.packages("pak")                  # development version
pak::pak("gcol33/tulpaMesh")
```

## Documentation

- [Quick Start](https://gillescolling.com/tulpaMesh/articles/quickstart.html)
- [Spatial Workflows](https://gillescolling.com/tulpaMesh/articles/workflows.html)
- [Spherical and Temporal Meshes](https://gillescolling.com/tulpaMesh/articles/advanced.html)
- [Full Reference](https://gillescolling.com/tulpaMesh/reference/)

## Support

> "Software is like sex: it's better when it's free." — Linus Torvalds

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
