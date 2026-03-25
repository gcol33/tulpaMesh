test_that("icosahedral mesh has correct vertex/triangle counts", {
  # Level 0: 12 vertices, 20 triangles
  mesh0 <- tulpa_mesh_sphere(subdivisions = 0)
  expect_equal(mesh0$n_vertices, 12)
  expect_equal(mesh0$n_triangles, 20)

  # Level 1: 42 vertices, 80 triangles
  mesh1 <- tulpa_mesh_sphere(subdivisions = 1)
  expect_equal(mesh1$n_vertices, 42)
  expect_equal(mesh1$n_triangles, 80)

  # Level 2: 162 vertices, 320 triangles
  mesh2 <- tulpa_mesh_sphere(subdivisions = 2)
  expect_equal(mesh2$n_vertices, 162)
  expect_equal(mesh2$n_triangles, 320)
})

test_that("all vertices lie on the unit sphere", {
  mesh <- tulpa_mesh_sphere(subdivisions = 2)
  radii <- sqrt(rowSums(mesh$vertices^2))
  expect_equal(radii, rep(1, mesh$n_vertices), tolerance = 1e-12)
})

test_that("custom radius works", {
  mesh <- tulpa_mesh_sphere(subdivisions = 1, radius = 6371)
  radii <- sqrt(rowSums(mesh$vertices^2))
  expect_equal(radii, rep(6371, mesh$n_vertices), tolerance = 1e-8)
})

test_that("lonlat coordinates are valid", {
  mesh <- tulpa_mesh_sphere(subdivisions = 2)
  expect_equal(ncol(mesh$lonlat), 2)
  expect_equal(nrow(mesh$lonlat), mesh$n_vertices)
  # Latitude in [-90, 90]
  expect_true(all(mesh$lonlat[, 2] >= -90 - 1e-10))
  expect_true(all(mesh$lonlat[, 2] <= 90 + 1e-10))
  # Longitude in [-180, 180]
  expect_true(all(mesh$lonlat[, 1] >= -180 - 1e-10))
  expect_true(all(mesh$lonlat[, 1] <= 180 + 1e-10))
})

test_that("data points can be inserted", {
  obs <- cbind(lon = c(0, 45, -90, 180), lat = c(0, 45, -30, 60))
  mesh <- tulpa_mesh_sphere(subdivisions = 2, coords = obs)
  # Should have extra vertices
  expect_equal(mesh$n_vertices, 162 + 4)
})

test_that("FEM matrices work on spherical mesh", {
  mesh <- tulpa_mesh_sphere(subdivisions = 2)
  fem <- fem_matrices(mesh)

  expect_equal(nrow(fem$C), mesh$n_vertices)
  expect_equal(ncol(fem$C), mesh$n_vertices)
  expect_true(Matrix::isSymmetric(fem$C, tol = 1e-10))
  expect_true(Matrix::isSymmetric(fem$G, tol = 1e-10))
  expect_true(all(Matrix::diag(fem$C) > 0))

  # Stiffness row sums should be ~0
  g_rowsums <- as.numeric(Matrix::rowSums(fem$G))
  expect_equal(g_rowsums, rep(0, length(g_rowsums)), tolerance = 1e-10)
})

test_that("total surface area of unit sphere mesh approximates 4*pi", {
  mesh <- tulpa_mesh_sphere(subdivisions = 3)
  fem <- fem_matrices(mesh, lumped = TRUE)

  total_area <- sum(fem$va)
  expect_equal(total_area, 4 * pi, tolerance = 0.01)  # ~0.3% error at level 3
})

test_that("sphere projection matrix has row sums = 1", {
  mesh <- tulpa_mesh_sphere(subdivisions = 2)
  obs <- cbind(lon = runif(10, -180, 180), lat = runif(10, -60, 60))
  fem <- fem_matrices(mesh, obs_coords = obs)

  row_sums <- as.numeric(Matrix::rowSums(fem$A))
  expect_equal(row_sums, rep(1, 10), tolerance = 1e-6)
})

test_that("print method shows sphere info", {
  mesh <- tulpa_mesh_sphere(subdivisions = 1)
  expect_output(print(mesh), "sphere")
})

test_that("Euler characteristic is 2 for sphere", {
  mesh <- tulpa_mesh_sphere(subdivisions = 2)
  # V - E + T = 2 for a closed sphere
  euler <- mesh$n_vertices - mesh$n_edges + mesh$n_triangles
  expect_equal(euler, 2)
})
