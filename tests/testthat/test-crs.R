test_that("mesh_crs returns NULL for meshes without CRS", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))
  expect_null(mesh_crs(mesh))
})

test_that("set_crs attaches and mesh_crs retrieves", {
  skip_if_not_installed("sf")

  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))
  mesh <- set_crs(mesh, 4326)

  crs <- mesh_crs(mesh)
  expect_false(is.null(crs))
  expect_equal(crs$epsg, 4326L)
})

test_that("set_crs with NULL removes CRS", {
  skip_if_not_installed("sf")

  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))
  mesh <- set_crs(mesh, 4326)
  mesh <- set_crs(mesh, NULL)
  expect_null(mesh_crs(mesh))
})

test_that("CRS propagates from sf boundary", {
  skip_if_not_installed("sf")

  poly <- sf::st_polygon(list(
    rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))
  ))
  sfc <- sf::st_sfc(poly, crs = 32633)

  set.seed(42)
  pts <- cbind(runif(20, 0.1, 0.9), runif(20, 0.1, 0.9))
  mesh <- tulpa_mesh(pts, boundary = sfc, extend = 0)

  crs <- mesh_crs(mesh)
  expect_false(is.null(crs))
  expect_equal(crs$epsg, 32633L)
})

test_that("CRS propagates through as_tulpa_mesh", {
  skip_if_not_installed("sf")
  skip_if_not_installed("fmesher")

  set.seed(42)
  pts <- cbind(runif(30), runif(30))
  fm <- fmesher::fm_mesh_2d(loc = pts, max.edge = c(0.3, 0.6))
  fm$crs <- sf::st_crs(4326)

  tm <- as_tulpa_mesh(fm)
  crs <- mesh_crs(tm)
  expect_false(is.null(crs))
  expect_equal(crs$epsg, 4326L)
})
