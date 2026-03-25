# Package index

## Mesh Generation

Create triangulated meshes from point coordinates

- [`tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh.md)
  [`print(`*`<tulpa_mesh>`*`)`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh.md)
  : Create a Triangular Mesh for SPDE Spatial Models
- [`tulpa_mesh_sphere()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh_sphere.md)
  : Create a Triangular Mesh on the Sphere
- [`tulpa_mesh_1d()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh_1d.md)
  [`print(`*`<tulpa_mesh_1d>`*`)`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh_1d.md)
  : Create a 1D Mesh for Temporal SPDE Models
- [`tulpa_mesh_graph()`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh_graph.md)
  [`print(`*`<tulpa_mesh_graph>`*`)`](https://github.com/gcol33/tulpaMesh/reference/tulpa_mesh_graph.md)
  : Create a Metric Graph Mesh

## FEM Assembly

Finite element matrices for SPDE models

- [`fem_matrices()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices.md)
  : Compute FEM Matrices from a Mesh
- [`fem_matrices_p2()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices_p2.md)
  : Compute P2 (Quadratic) FEM Matrices
- [`fem_matrices_nonstationary()`](https://github.com/gcol33/tulpaMesh/reference/fem_matrices_nonstationary.md)
  : Non-stationary FEM Matrices with Spatially Varying Parameters

## Mesh Quality

Diagnostics and visualization

- [`mesh_quality()`](https://github.com/gcol33/tulpaMesh/reference/mesh_quality.md)
  : Per-Triangle Mesh Quality Metrics
- [`mesh_summary()`](https://github.com/gcol33/tulpaMesh/reference/mesh_summary.md)
  : Print Mesh Quality Summary
- [`plot(`*`<tulpa_mesh>`*`)`](https://github.com/gcol33/tulpaMesh/reference/plot.tulpa_mesh.md)
  : Plot a Triangular Mesh

## Mesh Operations

Manipulate and refine meshes

- [`subdivide_mesh()`](https://github.com/gcol33/tulpaMesh/reference/subdivide_mesh.md)
  : Subdivide a Mesh
- [`subset_mesh()`](https://github.com/gcol33/tulpaMesh/reference/subset_mesh.md)
  : Extract a Submesh from Triangle Indices
- [`mesh_components()`](https://github.com/gcol33/tulpaMesh/reference/mesh_components.md)
  : Identify Disconnected Mesh Components
- [`refine_mesh()`](https://github.com/gcol33/tulpaMesh/reference/refine_mesh.md)
  : Adaptively Refine a Mesh Based on Error Indicators
- [`barrier_triangles()`](https://github.com/gcol33/tulpaMesh/reference/barrier_triangles.md)
  : Identify Barrier Triangles

## Conversion and Metadata

Interoperability and coordinate systems

- [`as_tulpa_mesh()`](https://github.com/gcol33/tulpaMesh/reference/as_tulpa_mesh.md)
  : Convert to a tulpa_mesh Object
- [`mesh_crs()`](https://github.com/gcol33/tulpaMesh/reference/mesh_crs.md)
  [`set_crs()`](https://github.com/gcol33/tulpaMesh/reference/mesh_crs.md)
  : Get or Set the CRS of a Mesh
