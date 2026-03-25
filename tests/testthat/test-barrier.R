test_that("barrier_triangles identifies triangles inside barrier", {
  set.seed(42)
  coords <- cbind(runif(50, 0, 10), runif(50, 0, 10))
  mesh <- tulpa_mesh(coords)

  # Barrier polygon in center of domain
  barrier_poly <- rbind(c(3, 3), c(7, 3), c(7, 7), c(3, 7))
  bt <- barrier_triangles(mesh, barrier_poly)

  expect_type(bt, "logical")
  expect_length(bt, mesh$n_triangles)
  # Some triangles should be inside, some outside
  expect_true(any(bt))
  expect_true(any(!bt))
})

test_that("fem_matrices with barrier zeros out stiffness for barrier triangles", {
  set.seed(42)
  coords <- cbind(runif(50, 0, 10), runif(50, 0, 10))
  mesh <- tulpa_mesh(coords)

  barrier_poly <- rbind(c(3, 3), c(7, 3), c(7, 7), c(3, 7))
  bt <- barrier_triangles(mesh, barrier_poly)

  fem_full <- fem_matrices(mesh)
  fem_barr <- fem_matrices(mesh, barrier = bt)

  # Mass matrix should be unchanged
  expect_equal(as.matrix(fem_barr$C), as.matrix(fem_full$C))

  # Stiffness should have fewer nonzeros (barrier triangles excluded)
  nnz_full <- length(fem_full$G@x)
  nnz_barr <- length(fem_barr$G@x)
  if (any(bt)) {
    expect_lt(nnz_barr, nnz_full)
  }

  # Stiffness should still be symmetric
  expect_true(Matrix::isSymmetric(fem_barr$G, tol = 1e-10))
})

test_that("barrier_triangles with sf polygon", {
  skip_if_not_installed("sf")

  set.seed(42)
  coords <- cbind(runif(50, 0, 10), runif(50, 0, 10))
  mesh <- tulpa_mesh(coords)

  barrier_sf <- sf::st_polygon(list(
    rbind(c(3, 3), c(7, 3), c(7, 7), c(3, 7), c(3, 3))
  ))
  barrier_sfc <- sf::st_sfc(barrier_sf)

  bt <- barrier_triangles(mesh, barrier_sfc)
  expect_type(bt, "logical")
  expect_length(bt, mesh$n_triangles)
  expect_true(any(bt))
})

test_that("barrier with no barrier triangles returns full stiffness", {
  set.seed(42)
  coords <- cbind(runif(30, 0, 1), runif(30, 0, 1))
  mesh <- tulpa_mesh(coords)

  # Barrier outside mesh domain
  bt <- rep(FALSE, mesh$n_triangles)

  fem_full <- fem_matrices(mesh)
  fem_barr <- fem_matrices(mesh, barrier = bt)

  expect_equal(as.matrix(fem_barr$G), as.matrix(fem_full$G))
})

test_that("barrier with all triangles returns zero stiffness", {
  set.seed(42)
  coords <- cbind(runif(20, 0, 1), runif(20, 0, 1))
  mesh <- tulpa_mesh(coords)

  bt <- rep(TRUE, mesh$n_triangles)
  fem_barr <- fem_matrices(mesh, barrier = bt)

  # All stiffness should be zero
  expect_equal(max(abs(fem_barr$G)), 0)
})
