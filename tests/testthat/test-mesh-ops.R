# --- subdivide_mesh ---

test_that("subdivide_mesh quadruples triangle count", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))
  sub <- subdivide_mesh(mesh)

  expect_s3_class(sub, "tulpa_mesh")
  expect_equal(sub$n_triangles, mesh$n_triangles * 4)
  expect_gt(sub$n_vertices, mesh$n_vertices)
})

test_that("subdivide_mesh preserves total area", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))
  sub <- subdivide_mesh(mesh)

  q_orig <- mesh_quality(mesh)
  q_sub <- mesh_quality(sub)

  expect_equal(sum(q_sub$area), sum(q_orig$area), tolerance = 1e-10)
})

test_that("subdivide_mesh FEM matrices are valid", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))
  sub <- subdivide_mesh(mesh)

  fem <- fem_matrices(sub)
  expect_true(Matrix::isSymmetric(fem$C, tol = 1e-10))
  expect_true(Matrix::isSymmetric(fem$G, tol = 1e-10))
  expect_true(all(Matrix::diag(fem$C) > 0))
})

test_that("subdivide_mesh contains original vertices", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(10), runif(10)))
  sub <- subdivide_mesh(mesh)

  # Original vertices should be a subset of subdivided vertices
  for (i in seq_len(nrow(mesh$vertices))) {
    dists <- sqrt(rowSums((sub$vertices - matrix(mesh$vertices[i, ],
                  nrow = nrow(sub$vertices), ncol = 2, byrow = TRUE))^2))
    expect_true(min(dists) < 1e-12)
  }
})

# --- mesh_components ---

test_that("single connected mesh has one component", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))
  comps <- mesh_components(mesh)

  expect_length(comps, mesh$n_triangles)
  expect_equal(length(unique(comps)), 1)
})

test_that("mesh_components identifies separate regions", {
  # Create two separate clusters of points
  set.seed(42)
  pts1 <- cbind(runif(20, 0, 1), runif(20, 0, 1))
  pts2 <- cbind(runif(20, 10, 11), runif(20, 10, 11))

  # Build boundary around first cluster
  bnd1 <- rbind(c(-0.1, -0.1), c(1.1, -0.1), c(1.1, 1.1), c(-0.1, 1.1))
  mesh1 <- tulpa_mesh(pts1, boundary = bnd1, extend = 0)

  # A single connected mesh should have 1 component
  comps1 <- mesh_components(mesh1)
  expect_equal(max(comps1), 1)
})

# --- subset_mesh ---

test_that("subset_mesh with integer indices works", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))

  # Keep first half of triangles
  half <- seq_len(floor(mesh$n_triangles / 2))
  sub <- subset_mesh(mesh, half)

  expect_s3_class(sub, "tulpa_mesh")
  expect_equal(sub$n_triangles, length(half))
  expect_lte(sub$n_vertices, mesh$n_vertices)
  # All triangle indices should be valid
  expect_true(all(sub$triangles >= 1))
  expect_true(all(sub$triangles <= sub$n_vertices))
})

test_that("subset_mesh with logical vector works", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))

  keep <- rep(FALSE, mesh$n_triangles)
  keep[1:10] <- TRUE
  sub <- subset_mesh(mesh, keep)

  expect_equal(sub$n_triangles, 10)
})

test_that("subset_mesh FEM matrices work on submesh", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))
  sub <- subset_mesh(mesh, seq_len(mesh$n_triangles))  # keep all

  fem_orig <- fem_matrices(mesh)
  fem_sub <- fem_matrices(sub)

  # Same mesh, should have same matrices
  expect_equal(as.matrix(fem_sub$C), as.matrix(fem_orig$C), tolerance = 1e-10)
})

test_that("subset_mesh validates inputs", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))

  expect_error(subset_mesh(mesh, c(0, 1)), "out of range")
  expect_error(subset_mesh(mesh, c(1, 999999)), "out of range")
  expect_error(subset_mesh(mesh, rep(TRUE, 3)), "n_triangles")
})

test_that("subset_mesh + barrier_triangles workflow", {
  set.seed(42)
  coords <- cbind(runif(50, 0, 10), runif(50, 0, 10))
  mesh <- tulpa_mesh(coords)

  # Identify barrier triangles, then extract non-barrier submesh
  barrier_poly <- rbind(c(3, 3), c(7, 3), c(7, 7), c(3, 7))
  bt <- barrier_triangles(mesh, barrier_poly)

  sub <- subset_mesh(mesh, !bt)
  expect_s3_class(sub, "tulpa_mesh")
  expect_lt(sub$n_triangles, mesh$n_triangles)

  fem <- fem_matrices(sub)
  expect_true(Matrix::isSymmetric(fem$G, tol = 1e-10))
})
