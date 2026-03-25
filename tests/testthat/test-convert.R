skip_if_not_installed("fmesher")

test_that("as_tulpa_mesh converts fm_mesh_2d correctly", {
  set.seed(42)
  pts <- cbind(runif(50), runif(50))
  fm <- fmesher::fm_mesh_2d(loc = pts, max.edge = c(0.3, 0.6))

  tm <- as_tulpa_mesh(fm)

  expect_s3_class(tm, "tulpa_mesh")
  expect_equal(tm$n_vertices, fm$n)
  expect_equal(tm$n_triangles, nrow(fm$graph$tv))
  expect_equal(ncol(tm$vertices), 2)
  expect_equal(ncol(tm$triangles), 3)
  expect_equal(ncol(tm$edges), 2)
})

test_that("FEM matrices from converted mesh match fmesher", {
  set.seed(42)
  pts <- cbind(runif(50), runif(50))
  fm <- fmesher::fm_mesh_2d(loc = pts, max.edge = c(0.3, 0.6))
  tm <- as_tulpa_mesh(fm)

  tp_fem <- fem_matrices(tm)
  fm_fem <- fmesher::fm_fem(fm, order = 2)

  expect_lt(max(abs(tp_fem$C - fm_fem$c1)), 1e-12)
  expect_lt(max(abs(tp_fem$G - fm_fem$g1)), 1e-12)
})

test_that("projection matrix from converted mesh matches fmesher", {
  set.seed(42)
  pts <- cbind(runif(50), runif(50))
  fm <- fmesher::fm_mesh_2d(loc = pts, max.edge = c(0.3, 0.6))
  tm <- as_tulpa_mesh(fm)

  obs <- cbind(runif(15, 0.1, 0.9), runif(15, 0.1, 0.9))
  tp_A <- fem_matrices(tm, obs_coords = obs)$A
  fm_A <- fmesher::fm_basis(fm, loc = obs)

  expect_lt(max(abs(tp_A - fm_A)), 1e-12)
})
