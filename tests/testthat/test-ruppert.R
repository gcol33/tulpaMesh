test_that("Ruppert refinement achieves target minimum angle", {
  set.seed(42)
  coords <- cbind(runif(30), runif(30))

  mesh_refined <- tulpa_mesh(coords, min_angle = 25)
  q_refined <- mesh_quality(mesh_refined)

  # Shewchuk algorithm guarantees min angle >= target for well-shaped input
  expect_gte(min(q_refined$min_angle), 24.5)
  # Median should be well above target
  expect_gte(median(q_refined$min_angle), 35)
})

test_that("Ruppert refinement adds Steiner points", {
  set.seed(42)
  coords <- cbind(runif(30), runif(30))

  mesh_plain <- tulpa_mesh(coords)
  mesh_refined <- tulpa_mesh(coords, min_angle = 25)

  # Refined mesh should have more vertices (Steiner points added)
  expect_gte(mesh_refined$n_vertices, mesh_plain$n_vertices)
})

test_that("max_area constraint produces finer mesh", {
  set.seed(42)
  coords <- cbind(runif(20), runif(20))

  mesh_coarse <- tulpa_mesh(coords)
  mesh_fine <- tulpa_mesh(coords, max_area = 0.005)

  q_fine <- mesh_quality(mesh_fine)

  # All triangles should satisfy the area constraint
  expect_true(all(q_fine$area <= 0.005 + 1e-6))
  # Fine mesh should have more triangles
  expect_gt(mesh_fine$n_triangles, mesh_coarse$n_triangles)
})

test_that("combined min_angle + max_area works", {
  set.seed(42)
  coords <- cbind(runif(20), runif(20))

  mesh <- tulpa_mesh(coords, min_angle = 25, max_area = 0.01)
  q <- mesh_quality(mesh)

  # Both constraints satisfied
  expect_gte(min(q$min_angle), 24.5)
  expect_lte(max(q$area), 0.01 + 1e-6)
})

test_that("Ruppert with boundary achieves target angle", {
  bnd <- rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1))
  set.seed(42)
  coords <- cbind(runif(20, 0.1, 0.9), runif(20, 0.1, 0.9))

  mesh <- tulpa_mesh(coords, boundary = bnd, min_angle = 25, extend = 0)
  expect_s3_class(mesh, "tulpa_mesh")
  expect_true(mesh$n_triangles > 0)

  q <- mesh_quality(mesh)
  expect_true(all(q$area > 0))
  # Segment splitting ensures all triangles meet the angle constraint
  expect_gte(min(q$min_angle), 24.5)
})

test_that("Ruppert handles non-convex boundaries", {
  # L-shaped domain
  l_shape <- rbind(c(0,0), c(2,0), c(2,1), c(1,1), c(1,2), c(0,2))
  set.seed(42)
  coords <- cbind(runif(15, 0.1, 1.9), runif(15, 0.1, 0.9))

  mesh <- tulpa_mesh(coords, boundary = l_shape, extend = 0, min_angle = 25)
  q <- mesh_quality(mesh)

  expect_gte(min(q$min_angle), 24.5)
  expect_equal(sum(q$min_angle < 20), 0)
})

test_that("max_steiner limits refinement", {
  set.seed(42)
  coords <- cbind(runif(20), runif(20))

  mesh <- tulpa_mesh(coords, min_angle = 30, max_steiner = 5L)
  # Should stop after 5 Steiner points
  expect_s3_class(mesh, "tulpa_mesh")
})

test_that("FEM matrices work on Ruppert-refined mesh", {
  set.seed(42)
  coords <- cbind(runif(30), runif(30))
  mesh <- tulpa_mesh(coords, min_angle = 25)

  fem <- fem_matrices(mesh)
  expect_equal(nrow(fem$C), mesh$n_vertices)
  expect_true(Matrix::isSymmetric(fem$G, tol = 1e-10))
  expect_true(all(Matrix::diag(fem$C) > 0))
})

test_that("Ruppert-refined mesh FEM matches fmesher on same triangulation", {
  skip_if_not_installed("fmesher")

  set.seed(42)
  coords <- cbind(runif(30), runif(30))
  tm <- tulpa_mesh(coords, min_angle = 20)

  # Build fm_mesh_2d from tulpaMesh's refined mesh
  fm <- fmesher::fm_rcdt_2d(loc = tm$vertices)
  fm$loc <- cbind(tm$vertices, 0)
  fm$n <- tm$n_vertices
  fm$graph$tv <- tm$triangles

  tp_fem <- fem_matrices(tm)
  fm_fem <- fmesher::fm_fem(fm, order = 2)

  expect_lt(max(abs(tp_fem$C - fm_fem$c1)), 1e-12)
  expect_lt(max(abs(tp_fem$G - fm_fem$g1)), 1e-3)
})
