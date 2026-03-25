test_that("non-stationary FEM with constant kappa/tau equals scaled stationary", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))
  n <- mesh$n_vertices

  # Constant kappa = 2, tau = 3
  kappa <- rep(2, n)
  tau <- rep(3, n)

  ns <- fem_matrices_nonstationary(mesh, kappa, tau)
  fem <- fem_matrices(mesh)

  # Ck should be kappa^2 * C = 4 * C
  expect_lt(max(abs(ns$Ck - 4 * fem$C)), 1e-12)

  # Gk should be kappa^2 * G = 4 * G
  expect_lt(max(abs(ns$Gk - 4 * fem$G)), 1e-12)

  # Ct should be tau^2 * C = 9 * C
  expect_lt(max(abs(ns$Ct - 9 * fem$C)), 1e-12)
})

test_that("non-stationary matrices are symmetric", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))
  n <- mesh$n_vertices

  kappa <- runif(n, 0.5, 5)
  tau <- runif(n, 0.1, 2)

  ns <- fem_matrices_nonstationary(mesh, kappa, tau)

  expect_true(Matrix::isSymmetric(ns$Ck, tol = 1e-10))
  expect_true(Matrix::isSymmetric(ns$Gk, tol = 1e-10))
  expect_true(Matrix::isSymmetric(ns$Ct, tol = 1e-10))
})

test_that("non-stationary mass matrices are positive definite", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))
  n <- mesh$n_vertices

  kappa <- runif(n, 1, 3)
  tau <- runif(n, 0.5, 2)

  ns <- fem_matrices_nonstationary(mesh, kappa, tau)

  expect_true(all(Matrix::diag(ns$Ck) > 0))
  expect_true(all(Matrix::diag(ns$Ct) > 0))
})

test_that("non-stationary stiffness row sums near zero", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))
  n <- mesh$n_vertices

  kappa <- runif(n, 0.5, 5)
  tau <- runif(n, 0.1, 2)

  ns <- fem_matrices_nonstationary(mesh, kappa, tau)

  gk_rowsums <- as.numeric(Matrix::rowSums(ns$Gk))
  expect_equal(gk_rowsums, rep(0, n), tolerance = 1e-10)
})

test_that("non-stationary FEM validates inputs", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))

  expect_error(fem_matrices_nonstationary(mesh, 1:3, 1:3), "n_vertices")
  expect_error(fem_matrices_nonstationary("x", 1, 1), "tulpa_mesh")
})

test_that("non-stationary returns standard matrices too", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))
  n <- mesh$n_vertices

  ns <- fem_matrices_nonstationary(mesh, rep(1, n), rep(1, n))
  fem <- fem_matrices(mesh)

  # Standard C and G should match
  expect_lt(max(abs(ns$C - fem$C)), 1e-12)
  expect_lt(max(abs(ns$G - fem$G)), 1e-12)
  expect_true(Matrix::isDiagonal(ns$C0))
})
