# Cross-validate tulpaMesh against fmesher.
# tulpaMesh builds the mesh; fmesher's fm_fem / fm_basis provide reference
# FEM matrices on that same triangulation.

skip_if_not_installed("fmesher")

# ---------------------------------------------------------------------------
# Helper: wrap a tulpa_mesh as an fm_mesh_2d so fmesher can compute on it
# ---------------------------------------------------------------------------
tulpa_to_fm <- function(tm) {
  # fm_fem needs an fm_mesh_2d with $loc (N x 3), $n, $graph$tv
  # Build a minimal dummy, then overwrite its internals
  fm <- fmesher::fm_rcdt_2d(loc = tm$vertices)
  fm$loc      <- cbind(tm$vertices, 0)
  fm$n        <- tm$n_vertices
  fm$graph$tv <- tm$triangles
  fm
}

# ---------------------------------------------------------------------------
# Shared fixture: tulpaMesh builds the mesh, fmesher cross-checks FEM
# ---------------------------------------------------------------------------
set.seed(42)
pts <- cbind(runif(50), runif(50))
tm   <- tulpa_mesh(pts)
fm   <- tulpa_to_fm(tm)

tp_fem <- fem_matrices(tm)
fm_fem <- fmesher::fm_fem(fm, order = 2)

# ---------------------------------------------------------------------------
# Mass matrix: tulpaMesh C == fmesher c1 (consistent mass)
# ---------------------------------------------------------------------------
test_that("consistent mass matrix C matches fmesher c1", {
  diff <- max(abs(tp_fem$C - fm_fem$c1))
  expect_lt(diff, 1e-12,
            label = sprintf("max |C - c1| = %.2e", diff))
})

# ---------------------------------------------------------------------------
# Stiffness matrix: tulpaMesh G == fmesher g1
# ---------------------------------------------------------------------------
test_that("stiffness matrix G matches fmesher g1", {
  diff <- max(abs(tp_fem$G - fm_fem$g1))
  expect_lt(diff, 1e-12,
            label = sprintf("max |G - g1| = %.2e", diff))
})

# ---------------------------------------------------------------------------
# Lumped mass: diagonal of fmesher c0 == rowSums(tulpaMesh C)
# ---------------------------------------------------------------------------
test_that("row sums of C equal fmesher lumped mass c0 diagonal", {
  c0_diag   <- Matrix::diag(fm_fem$c0)
  c_rowsums <- as.numeric(Matrix::rowSums(tp_fem$C))
  diff <- max(abs(c_rowsums - c0_diag))
  expect_lt(diff, 1e-12,
            label = sprintf("max |rowSums(C) - diag(c0)| = %.2e", diff))
})

# ---------------------------------------------------------------------------
# Stiffness row sums ~ 0 (divergence theorem)
# ---------------------------------------------------------------------------
test_that("stiffness matrix G row sums are near zero", {
  g_rowsums <- as.numeric(Matrix::rowSums(tp_fem$G))
  expect_equal(g_rowsums, rep(0, length(g_rowsums)), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Projection matrix: tulpaMesh A == fmesher fm_basis for interior points
# ---------------------------------------------------------------------------
test_that("projection matrix A matches fmesher fm_basis for interior points", {
  set.seed(123)
  obs <- cbind(runif(20, 0.1, 0.9), runif(20, 0.1, 0.9))

  tp_A <- fem_matrices(tm, obs_coords = obs)$A
  fm_A <- fmesher::fm_basis(fm, loc = obs)

  diff <- max(abs(tp_A - fm_A))
  expect_lt(diff, 1e-12,
            label = sprintf("max |A_tulpa - A_fmesher| = %.2e", diff))
})

# ---------------------------------------------------------------------------
# Projection at mesh vertices: should be identity
# ---------------------------------------------------------------------------
test_that("projection at mesh vertices recovers identity", {
  tp_A <- fem_matrices(tm, obs_coords = tm$vertices)$A
  fm_A <- fmesher::fm_basis(fm, loc = tm$vertices)

  diff <- max(abs(tp_A - fm_A))
  expect_lt(diff, 1e-12,
            label = sprintf("max |A_vertex_tulpa - A_vertex_fmesher| = %.2e", diff))
})

# ---------------------------------------------------------------------------
# Projection row sums == 1
# ---------------------------------------------------------------------------
test_that("projection matrix row sums are 1", {
  set.seed(456)
  obs <- cbind(runif(30, 0.05, 0.95), runif(30, 0.05, 0.95))
  tp_A <- fem_matrices(tm, obs_coords = obs)$A
  row_sums <- as.numeric(Matrix::rowSums(tp_A))
  expect_equal(row_sums, rep(1, 30), tolerance = 1e-10)
})

# ---------------------------------------------------------------------------
# Lumped mass: C0 diagonal matches fmesher c0
# ---------------------------------------------------------------------------
test_that("lumped mass C0 matches fmesher c0", {
  tp_lumped <- fem_matrices(tm, lumped = TRUE)

  expect_true(Matrix::isDiagonal(tp_lumped$C0))
  diff <- max(abs(tp_lumped$C0 - fm_fem$c0))
  expect_lt(diff, 1e-12,
            label = sprintf("max |C0 - c0| = %.2e", diff))
})

# ---------------------------------------------------------------------------
# Vertex areas match fmesher va
# ---------------------------------------------------------------------------
test_that("vertex areas match fmesher va", {
  tp_lumped <- fem_matrices(tm, lumped = TRUE)

  diff <- max(abs(tp_lumped$va - fm_fem$va))
  expect_lt(diff, 1e-12,
            label = sprintf("max |va_tulpa - va_fmesher| = %.2e", diff))
})

# ---------------------------------------------------------------------------
# Triangle areas match fmesher ta
# ---------------------------------------------------------------------------
test_that("triangle areas match fmesher ta", {
  tp_lumped <- fem_matrices(tm, lumped = TRUE)

  diff <- max(abs(tp_lumped$ta - fm_fem$ta))
  expect_lt(diff, 1e-12,
            label = sprintf("max |ta_tulpa - ta_fmesher| = %.2e", diff))
})

# ---------------------------------------------------------------------------
# Symmetry
# ---------------------------------------------------------------------------
test_that("C and G are symmetric", {
  expect_true(Matrix::isSymmetric(tp_fem$C, tol = 1e-12))
  expect_true(Matrix::isSymmetric(tp_fem$G, tol = 1e-12))
})

# ---------------------------------------------------------------------------
# Mass matrix positive definite
# ---------------------------------------------------------------------------
test_that("mass matrix C is positive definite", {
  skip_if(tm$n_vertices > 500, "mesh too large for dense eigendecomposition")
  eigs <- eigen(as.matrix(tp_fem$C), symmetric = TRUE, only.values = TRUE)$values
  expect_true(all(eigs > 0),
              label = sprintf("min eigenvalue = %.2e", min(eigs)))
})
