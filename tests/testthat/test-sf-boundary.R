skip_if_not_installed("sf")

test_that("sf POLYGON boundary works", {
  # Simple square polygon
  poly <- sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))))
  sfc <- sf::st_sfc(poly)

  set.seed(42)
  pts <- cbind(runif(30, 0.1, 0.9), runif(30, 0.1, 0.9))
  mesh <- tulpa_mesh(pts, boundary = sfc, extend = 0)

  expect_s3_class(mesh, "tulpa_mesh")
  expect_true(mesh$n_vertices >= 30)

  # All vertices should be within the polygon (or very close to boundary)
  expect_true(all(mesh$vertices[, 1] >= -0.01))
  expect_true(all(mesh$vertices[, 1] <= 1.01))
  expect_true(all(mesh$vertices[, 2] >= -0.01))
  expect_true(all(mesh$vertices[, 2] <= 1.01))
})

test_that("sf POLYGON with hole creates constrained mesh", {
  # Outer square with inner square hole
  outer <- rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0))
  hole <- rbind(c(3, 3), c(7, 3), c(7, 7), c(3, 7), c(3, 3))
  poly <- sf::st_polygon(list(outer, hole))
  sfc <- sf::st_sfc(poly)

  set.seed(42)
  pts <- cbind(runif(40, 0.5, 9.5), runif(40, 0.5, 9.5))
  # Remove points inside the hole
  inside_hole <- pts[, 1] > 3 & pts[, 1] < 7 & pts[, 2] > 3 & pts[, 2] < 7
  pts <- pts[!inside_hole, , drop = FALSE]

  mesh <- tulpa_mesh(pts, boundary = sfc, extend = 0)
  expect_s3_class(mesh, "tulpa_mesh")
  expect_true(mesh$n_triangles > 0)
})

test_that("sf MULTIPOLYGON works", {
  # Two separate squares
  poly1 <- sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))))
  poly2 <- sf::st_polygon(list(rbind(c(3, 0), c(4, 0), c(4, 1), c(3, 1), c(3, 0))))
  mpoly <- sf::st_multipolygon(list(poly1, poly2))
  sfc <- sf::st_sfc(mpoly)

  set.seed(42)
  pts <- cbind(c(runif(15, 0.1, 0.9), runif(15, 3.1, 3.9)),
               runif(30, 0.1, 0.9))
  mesh <- tulpa_mesh(pts, boundary = sfc, extend = 0)

  expect_s3_class(mesh, "tulpa_mesh")
  expect_true(mesh$n_triangles > 0)
})

test_that("sf_to_boundary extracts correct structure", {
  # Polygon with hole
  outer <- rbind(c(0, 0), c(10, 0), c(10, 10), c(0, 10), c(0, 0))
  hole <- rbind(c(3, 3), c(7, 3), c(7, 7), c(3, 7), c(3, 3))
  poly <- sf::st_polygon(list(outer, hole))
  sfc <- sf::st_sfc(poly)

  result <- tulpaMesh:::sf_to_boundary(sfc)

  expect_type(result, "list")
  expect_named(result, c("boundary", "holes"))

  # Outer ring should have 4 vertices (closing point removed)
  expect_equal(nrow(result$boundary), 4)
  # One hole
  expect_length(result$holes, 1)
  expect_equal(nrow(result$holes[[1]]), 4)
})

test_that("sf sfg object works directly", {
  poly <- sf::st_polygon(list(rbind(c(0, 0), c(1, 0), c(1, 1), c(0, 1), c(0, 0))))

  set.seed(42)
  pts <- cbind(runif(20, 0.1, 0.9), runif(20, 0.1, 0.9))
  mesh <- tulpa_mesh(pts, boundary = poly, extend = 0)

  expect_s3_class(mesh, "tulpa_mesh")
})
