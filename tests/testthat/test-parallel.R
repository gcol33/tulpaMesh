test_that("parallel FEM matches sequential FEM", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(100), runif(100)))

  fem_seq <- fem_matrices(mesh, parallel = FALSE)
  fem_par <- fem_matrices(mesh, parallel = TRUE)

  expect_lt(max(abs(fem_seq$C - fem_par$C)), 1e-12)
  expect_lt(max(abs(fem_seq$G - fem_par$G)), 1e-12)
})

test_that("parallel FEM matrices are symmetric", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(50), runif(50)))
  fem <- fem_matrices(mesh, parallel = TRUE)

  expect_true(Matrix::isSymmetric(fem$C, tol = 1e-10))
  expect_true(Matrix::isSymmetric(fem$G, tol = 1e-10))
})

test_that("parallel FEM stiffness row sums are zero", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(50), runif(50)))
  fem <- fem_matrices(mesh, parallel = TRUE)

  g_rowsums <- as.numeric(Matrix::rowSums(fem$G))
  expect_equal(g_rowsums, rep(0, mesh$n_vertices), tolerance = 1e-10)
})

test_that("parallel FEM with lumped mass works", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(50), runif(50)))

  fem_par <- fem_matrices(mesh, parallel = TRUE, lumped = TRUE)
  fem_seq <- fem_matrices(mesh, parallel = FALSE, lumped = TRUE)

  expect_lt(max(abs(fem_par$C0 - fem_seq$C0)), 1e-12)
  expect_equal(fem_par$va, fem_seq$va, tolerance = 1e-12)
})

test_that("parallel FEM cross-validates against fmesher", {
  skip_if_not_installed("fmesher")

  set.seed(42)
  pts <- cbind(runif(50), runif(50))
  tm <- tulpa_mesh(pts)

  fm <- fmesher::fm_rcdt_2d(loc = tm$vertices)
  fm$loc <- cbind(tm$vertices, 0)
  fm$n <- tm$n_vertices
  fm$graph$tv <- tm$triangles

  tp_fem <- fem_matrices(tm, parallel = TRUE)
  fm_fem <- fmesher::fm_fem(fm, order = 2)

  expect_lt(max(abs(tp_fem$C - fm_fem$c1)), 1e-12)
  expect_lt(max(abs(tp_fem$G - fm_fem$g1)), 1e-12)
})
