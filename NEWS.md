# tulpaMesh 0.1.1

* Removed single quotes from person names and algorithm names in DESCRIPTION
  Title and Description fields per CRAN policy.
* Added Artem Amirkhanov (CDT library) and William C. Lenthe (predicates.h)
  to Authors@R with ctb and cph roles.
* Updated inst/COPYRIGHTS to separately document predicates.h (BSD-3-Clause).

# tulpaMesh 0.1.0

Initial CRAN release.

## Mesh Generation

* `tulpa_mesh()`: constrained Delaunay triangulation with formula interface,
  sf boundary support (POLYGON, MULTIPOLYGON, holes), hexagonal lattice
  refinement, Ruppert refinement (`min_angle`, `max_area`), and vertex
  deduplication.
* `tulpa_mesh_sphere()`: icosahedral geodesic meshes with 3D surface FEM.
* `tulpa_mesh_1d()`: temporal meshes with tridiagonal FEM matrices.
* `tulpa_mesh_graph()`: metric graph meshes for network SPDE models.

## FEM Assembly

* `fem_matrices()`: P1 linear FEM producing mass (C), stiffness (G), and
  projection (A) sparse matrices. Supports lumped mass, barrier models,
  and parallel assembly via RcppParallel.
* `fem_matrices_p2()`: P2 quadratic FEM with 6-node triangular elements.
* `fem_matrices_nonstationary()`: spatially varying kappa/tau weighting.

## Mesh Operations

* `mesh_quality()`, `mesh_summary()`, `plot.tulpa_mesh()`: diagnostics.
* `subdivide_mesh()`, `subset_mesh()`, `mesh_components()`: mesh manipulation.
* `refine_mesh()`: adaptive refinement from error indicators.
* `as_tulpa_mesh()`: convert fmesher/INLA meshes.
* `barrier_triangles()`: identify barrier regions from polygons.
* `mesh_crs()`, `set_crs()`: CRS support.
