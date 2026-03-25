test_that("tulpa_mesh creates valid triangulation from random points", {
  set.seed(42)
  coords <- cbind(runif(50), runif(50))
  mesh <- tulpa_mesh(coords)

  expect_s3_class(mesh, "tulpa_mesh")
  expect_true(mesh$n_vertices >= 50)
  expect_true(mesh$n_triangles > 0)
  expect_true(mesh$n_edges > 0)
  expect_equal(ncol(mesh$vertices), 2)
  expect_equal(ncol(mesh$triangles), 3)
  expect_equal(ncol(mesh$edges), 2)

  # All triangle indices should be valid
  expect_true(all(mesh$triangles >= 1))
  expect_true(all(mesh$triangles <= mesh$n_vertices))
})

test_that("tulpa_mesh works with formula interface", {
  set.seed(42)
  df <- data.frame(x = runif(30), y = runif(30), z = rnorm(30))
  mesh <- tulpa_mesh(~ x + y, data = df)

  expect_s3_class(mesh, "tulpa_mesh")
  expect_true(mesh$n_vertices >= 30)
})

test_that("fem_matrices returns valid sparse matrices", {
  set.seed(42)
  coords <- cbind(runif(30), runif(30))
  mesh <- tulpa_mesh(coords)

  fem <- fem_matrices(mesh)

  expect_true(inherits(fem$C, "Matrix"))
  expect_true(inherits(fem$G, "Matrix"))
  expect_equal(nrow(fem$C), mesh$n_vertices)
  expect_equal(ncol(fem$C), mesh$n_vertices)
  expect_equal(nrow(fem$G), mesh$n_vertices)

  # Mass matrix should be positive (all diagonal entries > 0)
  expect_true(all(Matrix::diag(fem$C) > 0))

  # Stiffness matrix should be symmetric
  expect_true(Matrix::isSymmetric(fem$G, tol = 1e-10))
})

test_that("fem_matrices computes projection matrix", {
  set.seed(42)
  coords <- cbind(runif(30), runif(30))
  mesh <- tulpa_mesh(coords)

  # Observation points (some inside mesh)
  obs <- cbind(runif(10, 0.1, 0.9), runif(10, 0.1, 0.9))
  fem <- fem_matrices(mesh, obs_coords = obs)

  expect_equal(nrow(fem$A), 10)
  expect_equal(ncol(fem$A), mesh$n_vertices)

  # Each row of A should sum to ~1 (barycentric weights)
  row_sums <- Matrix::rowSums(fem$A)
  expect_equal(row_sums, rep(1, 10), tolerance = 1e-6)
})

test_that("print method works", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))
  expect_output(print(mesh), "tulpa_mesh")
})
