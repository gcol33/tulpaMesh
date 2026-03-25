test_that("P2 FEM matrices have correct dimensions", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))
  p2 <- fem_matrices_p2(mesh)

  # Should have more nodes than P1
  expect_gt(p2$n_mesh, mesh$n_vertices)
  expect_equal(p2$n_vertices, mesh$n_vertices)
  expect_gt(p2$n_midpoints, 0)

  # Matrices should be n_total x n_total
  expect_equal(nrow(p2$C), p2$n_mesh)
  expect_equal(ncol(p2$C), p2$n_mesh)
  expect_equal(nrow(p2$G), p2$n_mesh)
})

test_that("P2 mass matrix is symmetric positive definite", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(15), runif(15)))
  p2 <- fem_matrices_p2(mesh)

  expect_true(Matrix::isSymmetric(p2$C, tol = 1e-10))

  # Check positive definiteness
  skip_if(p2$n_mesh > 500, "too large for dense eigendecomposition")
  eigs <- eigen(as.matrix(p2$C), symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > -1e-10),
              label = sprintf("min eigenvalue = %.2e", min(eigs)))
})

test_that("P2 stiffness matrix is symmetric with zero row sums", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))
  p2 <- fem_matrices_p2(mesh)

  expect_true(Matrix::isSymmetric(p2$G, tol = 1e-10))

  g_rowsums <- as.numeric(Matrix::rowSums(p2$G))
  expect_equal(g_rowsums, rep(0, p2$n_mesh), tolerance = 1e-10)
})

test_that("P2 total area equals P1 total area", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))

  p1 <- fem_matrices(mesh)
  p2 <- fem_matrices_p2(mesh)

  # sum(C) = total area for both P1 and P2
  area_p1 <- sum(p1$C)
  area_p2 <- sum(p2$C)
  expect_equal(area_p2, area_p1, tolerance = 1e-10)
})

test_that("P2 triangles6 has correct structure", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(15), runif(15)))
  p2 <- fem_matrices_p2(mesh)

  expect_equal(ncol(p2$triangles6), 6)
  expect_equal(nrow(p2$triangles6), mesh$n_triangles)

  # First 3 columns should be the original triangle vertices
  expect_equal(p2$triangles6[, 1:3], mesh$triangles)

  # Midpoint indices should be > n_vertices
  expect_true(all(p2$triangles6[, 4:6] > mesh$n_vertices))
})

test_that("midpoints are at edge centers", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(10), runif(10)))
  p2 <- fem_matrices_p2(mesh)

  # Check first triangle's midpoints
  t1 <- p2$triangles6[1, ]
  v0 <- p2$vertices[t1[1], ]
  v1 <- p2$vertices[t1[2], ]
  v2 <- p2$vertices[t1[3], ]
  m01 <- p2$vertices[t1[4], ]
  m12 <- p2$vertices[t1[5], ]
  m20 <- p2$vertices[t1[6], ]

  expect_equal(m01, (v0 + v1) / 2, tolerance = 1e-12)
  expect_equal(m12, (v1 + v2) / 2, tolerance = 1e-12)
  expect_equal(m20, (v2 + v0) / 2, tolerance = 1e-12)
})
