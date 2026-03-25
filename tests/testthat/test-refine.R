test_that("refine_mesh adds vertices for high-error triangles", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))

  # Simulate error: large for first half of triangles
  indicators <- rep(0.1, mesh$n_triangles)
  indicators[1:floor(mesh$n_triangles / 2)] <- 1.0

  refined <- refine_mesh(mesh, indicators, threshold = 0.5)

  expect_s3_class(refined, "tulpa_mesh")
  expect_gt(refined$n_vertices, mesh$n_vertices)
  expect_gt(refined$n_triangles, mesh$n_triangles)
})

test_that("refine_mesh with fraction refines specified proportion", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))
  indicators <- runif(mesh$n_triangles)

  refined <- refine_mesh(mesh, indicators, fraction = 0.3)

  expect_gt(refined$n_vertices, mesh$n_vertices)
})

test_that("refine_mesh with no high-error triangles returns same size", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))
  indicators <- rep(0.1, mesh$n_triangles)

  refined <- refine_mesh(mesh, indicators, threshold = 1.0)

  # No refinement — same number of input points
  expect_equal(refined$n_vertices, mesh$n_vertices)
})

test_that("refine_mesh respects min_area", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))
  indicators <- rep(1.0, mesh$n_triangles)

  # Set min_area very large — no triangle should be refined
  refined <- refine_mesh(mesh, indicators, threshold = 0.5, min_area = 100)

  expect_equal(refined$n_vertices, mesh$n_vertices)
})

test_that("FEM matrices work on refined mesh", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))
  indicators <- runif(mesh$n_triangles)

  refined <- refine_mesh(mesh, indicators, fraction = 0.5)
  fem <- fem_matrices(refined)

  expect_equal(nrow(fem$C), refined$n_vertices)
  expect_true(Matrix::isSymmetric(fem$G, tol = 1e-10))
})

test_that("refine_mesh validates inputs", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))

  expect_error(refine_mesh(mesh, c(1, 2, 3)), "n_triangles")
  expect_error(refine_mesh("not a mesh", rep(1, 10)), "tulpa_mesh")
})
