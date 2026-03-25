test_that("tulpa_mesh_graph creates valid graph mesh", {
  # Simple cross-shaped network
  edges <- list(
    cbind(c(0, 1), c(0.5, 0.5)),  # horizontal
    cbind(c(0.5, 0.5), c(0, 1))   # vertical
  )
  g <- tulpa_mesh_graph(edges)

  expect_s3_class(g, "tulpa_mesh_graph")
  expect_gt(g$n_vertices, 0)
  expect_gt(g$n_segments, 0)
})

test_that("junction detection works", {
  # T-junction: 3 edges meeting at (0.5, 0.5)
  edges <- list(
    cbind(c(0, 0.5), c(0.5, 0.5)),
    cbind(c(0.5, 1), c(0.5, 0.5)),
    cbind(c(0.5, 0.5), c(0.5, 1))
  )
  g <- tulpa_mesh_graph(edges)

  # The junction point should have degree 3
  expect_true(any(g$degree == 3))
})

test_that("graph FEM mass matrix is symmetric", {
  edges <- list(
    cbind(c(0, 5), c(0, 0)),
    cbind(c(0, 0), c(0, 5)),
    cbind(c(0, 5), c(0, 5))
  )
  g <- tulpa_mesh_graph(edges, max_edge = 1)

  expect_true(Matrix::isSymmetric(g$C, tol = 1e-12))
  expect_true(all(Matrix::diag(g$C) > 0))
})

test_that("graph FEM stiffness has zero row sums", {
  edges <- list(
    cbind(seq(0, 10, by = 1), rep(0, 11)),
    cbind(rep(5, 6), seq(0, 5, by = 1))
  )
  g <- tulpa_mesh_graph(edges)

  expect_true(Matrix::isSymmetric(g$G, tol = 1e-12))
  g_rowsums <- as.numeric(Matrix::rowSums(g$G))
  expect_equal(g_rowsums, rep(0, g$n_vertices), tolerance = 1e-10)
})

test_that("max_edge subdivides long segments", {
  edges <- list(cbind(c(0, 10), c(0, 0)))
  g_coarse <- tulpa_mesh_graph(edges)
  g_fine <- tulpa_mesh_graph(edges, max_edge = 1)

  expect_gt(g_fine$n_vertices, g_coarse$n_vertices)
  expect_gt(g_fine$n_segments, g_coarse$n_segments)
})

test_that("total mass equals total edge length", {
  edges <- list(
    cbind(c(0, 3), c(0, 0)),   # length 3
    cbind(c(0, 0), c(0, 4))    # length 4
  )
  g <- tulpa_mesh_graph(edges)

  total_mass <- sum(g$C)
  expect_equal(total_mass, 7, tolerance = 1e-10)
})

test_that("print method works for graph", {
  edges <- list(cbind(c(0, 1), c(0, 0)))
  g <- tulpa_mesh_graph(edges)
  expect_output(print(g), "tulpa_mesh_graph")
})

test_that("sf LINESTRING input works", {
  skip_if_not_installed("sf")

  line1 <- sf::st_linestring(cbind(c(0, 1, 2), c(0, 1, 0)))
  line2 <- sf::st_linestring(cbind(c(1, 1), c(1, 2)))
  sfc <- sf::st_sfc(line1, line2)

  g <- tulpa_mesh_graph(sfc)
  expect_s3_class(g, "tulpa_mesh_graph")
  expect_gt(g$n_vertices, 0)
})

test_that("snap_tolerance merges close endpoints", {
  edges <- list(
    cbind(c(0, 1), c(0, 0)),
    cbind(c(1 + 1e-10, 2), c(0, 0))  # almost touching
  )
  g <- tulpa_mesh_graph(edges, snap_tolerance = 1e-8)

  # Should have 3 vertices (0,0), (1,0), (2,0) — not 4
  expect_equal(g$n_vertices, 3)
})
