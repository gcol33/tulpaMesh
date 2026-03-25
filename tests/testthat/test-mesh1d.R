test_that("tulpa_mesh_1d creates valid 1D mesh", {
  m <- tulpa_mesh_1d(1:10)
  expect_s3_class(m, "tulpa_mesh_1d")
  expect_true(m$n >= 10)
  expect_equal(length(m$knots), m$n)
  expect_true(all(diff(m$knots) > 0))  # sorted
})

test_that("1D FEM mass matrix is symmetric positive definite", {
  m <- tulpa_mesh_1d(seq(0, 1, by = 0.1), n_extend = 0)
  expect_true(Matrix::isSymmetric(m$C, tol = 1e-12))
  expect_true(all(Matrix::diag(m$C) > 0))
  eigs <- eigen(as.matrix(m$C), symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > 0))
})

test_that("1D FEM stiffness matrix has zero row sums", {
  m <- tulpa_mesh_1d(seq(0, 1, by = 0.1), n_extend = 0)
  expect_true(Matrix::isSymmetric(m$G, tol = 1e-12))
  g_rowsums <- as.numeric(Matrix::rowSums(m$G))
  expect_equal(g_rowsums, rep(0, m$n), tolerance = 1e-12)
})

test_that("1D lumped mass equals row sums of consistent mass", {
  m <- tulpa_mesh_1d(seq(0, 5, by = 0.5), n_extend = 0)
  c_rowsums <- as.numeric(Matrix::rowSums(m$C))
  c0_diag <- as.numeric(Matrix::diag(m$C0))
  expect_equal(c0_diag, c_rowsums, tolerance = 1e-12)
})

test_that("1D total mass equals domain length", {
  m <- tulpa_mesh_1d(seq(0, 10, by = 1), n_extend = 0)
  total <- sum(m$C)
  expect_equal(total, 10, tolerance = 1e-12)
})

test_that("extension adds knots beyond boundary", {
  m <- tulpa_mesh_1d(1:5, n_extend = 3)
  expect_true(min(m$knots) < 1)
  expect_true(max(m$knots) > 5)
  expect_true(m$n > 5)
})

test_that("n_extend = 0 gives exact knots", {
  m <- tulpa_mesh_1d(c(0, 1, 3, 7), n_extend = 0)
  expect_equal(m$knots, c(0, 1, 3, 7))
  expect_equal(m$n, 4)
})

test_that("irregular spacing works", {
  knots <- c(0, 0.1, 0.5, 2, 10)
  m <- tulpa_mesh_1d(knots, n_extend = 0)
  expect_equal(m$n, 5)
  expect_true(Matrix::isSymmetric(m$C, tol = 1e-12))
  expect_true(Matrix::isSymmetric(m$G, tol = 1e-12))
})

test_that("print method works", {
  m <- tulpa_mesh_1d(1:5)
  expect_output(print(m), "tulpa_mesh_1d")
})

test_that("matrices are tridiagonal", {
  m <- tulpa_mesh_1d(seq(0, 1, by = 0.1), n_extend = 0)
  C_dense <- as.matrix(m$C)
  G_dense <- as.matrix(m$G)
  n <- m$n
  # All entries more than 1 away from diagonal should be zero
  for (i in 1:n) {
    for (j in 1:n) {
      if (abs(i - j) > 1) {
        expect_equal(C_dense[i, j], 0)
        expect_equal(G_dense[i, j], 0)
      }
    }
  }
})
