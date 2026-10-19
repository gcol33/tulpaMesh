# Changelog

## tulpaMesh 0.1.3

CRAN release: 2026-06-15

- Parallel FEM assembly (`fem_matrices(parallel = TRUE)`) no longer
  references Intel TBB symbols directly. Each worker now accumulates
  triplets into a chunk-local buffer that is merged under a mutex, so
  the package builds on platforms where `RcppParallel` falls back to its
  TinyThread backend (for example Alpine Linux with musl).
- Fixed `tulpa_mesh(max_edge = ...)` collapsing to zero triangles for
  some `(max_edge, cutoff)` settings and point counts. Vertex
  deduplication could merge extended-hull boundary vertices that fell
  within `max_edge * 0.3` of each other, deleting a constraint edge and
  leaving the boundary loop open; the constrained triangulation then
  erased every triangle. Boundary and hole vertices are now protected
  from deduplication so the constraint loop stays closed.
- [`tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh.md)
  now errors loudly when the triangulation yields zero triangles (for
  example a collinear point set) instead of returning an empty mesh that
  silently produces all-zero FEM matrices.

## tulpaMesh 0.1.1

CRAN release: 2026-04-03

- Removed single quotes from person names and algorithm names in
  DESCRIPTION Title and Description fields per CRAN policy.
- Added Artem Amirkhanov (CDT library) and William C. Lenthe
  (predicates.h) to <Authors@R> with ctb and cph roles.
- Updated inst/COPYRIGHTS to separately document predicates.h
  (BSD-3-Clause).

## tulpaMesh 0.1.0

Initial CRAN release.

### Mesh Generation

- [`tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh.md):
  constrained Delaunay triangulation with formula interface, sf boundary
  support (POLYGON, MULTIPOLYGON, holes), hexagonal lattice refinement,
  Ruppert refinement (`min_angle`, `max_area`), and vertex
  deduplication.
- [`tulpa_mesh_sphere()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh_sphere.md):
  icosahedral geodesic meshes with 3D surface FEM.
- [`tulpa_mesh_1d()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh_1d.md):
  temporal meshes with tridiagonal FEM matrices.
- [`tulpa_mesh_graph()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh_graph.md):
  metric graph meshes for network SPDE models.

### FEM Assembly

- [`fem_matrices()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices.md):
  P1 linear FEM producing mass (C), stiffness (G), and projection (A)
  sparse matrices. Supports lumped mass, barrier models, and parallel
  assembly via RcppParallel.
- [`fem_matrices_p2()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices_p2.md):
  P2 quadratic FEM with 6-node triangular elements.
- [`fem_matrices_nonstationary()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices_nonstationary.md):
  spatially varying kappa/tau weighting.

### Mesh Operations

- [`mesh_quality()`](https://github.com/gcol33/tulpaMesh/reference/mesh_quality.md),
  [`mesh_summary()`](https://github.com/gcol33/tulpaMesh/reference/mesh_summary.md),
  [`plot.tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/plot.tulpa_mesh.md):
  diagnostics.
- [`subdivide_mesh()`](https://github.com/gcol33/tulpaMesh/reference/subdivide_mesh.md),
  [`subset_mesh()`](https://github.com/gcol33/tulpaMesh/reference/subset_mesh.md),
  [`mesh_components()`](https://github.com/gcol33/tulpaMesh/reference/mesh_components.md):
  mesh manipulation.
- [`refine_mesh()`](https://github.com/gcol33/tulpaMesh/reference/refine_mesh.md):
  adaptive refinement from error indicators.
- [`as_tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/as_tulpa_mesh.md):
  convert fmesher/INLA meshes.
- [`barrier_triangles()`](https://github.com/gcol33/tulpaMesh/reference/barrier_triangles.md):
  identify barrier regions from polygons.
- [`mesh_crs()`](https://github.com/gcol33/tulpaMesh/reference/mesh_crs.md),
  [`set_crs()`](https://github.com/gcol33/tulpaMesh/reference/mesh_crs.md):
  CRS support.
