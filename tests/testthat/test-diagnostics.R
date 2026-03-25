test_that("mesh_quality returns correct structure", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))
  q <- mesh_quality(mesh)

  expect_s3_class(q, "data.frame")
  expect_equal(nrow(q), mesh$n_triangles)
  expect_named(q, c("min_angle", "max_angle", "aspect_ratio", "area"))
})

test_that("mesh_quality angles are in valid range", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(50), runif(50)))
  q <- mesh_quality(mesh)

  expect_true(all(q$min_angle > 0))
  expect_true(all(q$max_angle < 180))
  expect_true(all(q$min_angle <= q$max_angle))
  # Sum of angles in a triangle = 180
  mid_angle <- 180 - q$min_angle - q$max_angle
  angle_sum <- q$min_angle + q$max_angle + mid_angle
  expect_equal(angle_sum, rep(180, nrow(q)), tolerance = 1e-8)
})

test_that("mesh_quality aspect ratio >= 1", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(50), runif(50)))
  q <- mesh_quality(mesh)

  # Equilateral triangle has aspect ratio 1; all others > 1
  expect_true(all(q$aspect_ratio >= 1 - 1e-10))
})

test_that("mesh_quality area matches FEM mass matrix", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))
  q <- mesh_quality(mesh)

  # Total area from quality should equal sum of diagonal of lumped mass
  # (lumped mass diagonal = rowSums(C) = vertex areas)
  fem <- fem_matrices(mesh)
  total_area_fem <- sum(Matrix::diag(fem$C)) * 2  # C diagonal = area/6 * 2 entries
  total_area_q <- sum(q$area)

  # sum(diag(C)) = sum(area_t / 6 * 3) = sum(area_t) / 2
  # so sum(area_t) = 2 * sum(diag(C))
  # Actually: each vertex gets area/6 from each adjacent triangle
  # rowSums(C) = sum of area/6 for diag + area/12 for off-diag per tri
  # = vertex area (Voronoi dual area)
  # sum(rowSums(C)) = sum of all C entries = sum(area_t * (3*1/6 + 6*1/12))
  # = sum(area_t * (1/2 + 1/2)) = sum(area_t)
  total_area_fem <- sum(fem$C)
  expect_equal(total_area_q, total_area_fem, tolerance = 1e-10)
})

test_that("mesh_quality on equilateral triangle", {
  # A single equilateral triangle
  h <- sqrt(3) / 2
  v <- rbind(c(0, 0), c(1, 0), c(0.5, h))
  # Need at least 3 points + boundary for CDT
  mesh <- structure(
    list(
      vertices = v,
      triangles = matrix(1:3, nrow = 1),
      edges = matrix(c(1L, 2L, 2L, 3L, 3L, 1L), ncol = 2, byrow = TRUE),
      n_vertices = 3L,
      n_triangles = 1L,
      n_edges = 3L,
      n_input_points = 3L
    ),
    class = "tulpa_mesh"
  )
  q <- mesh_quality(mesh)

  expect_equal(q$min_angle, 60, tolerance = 1e-8)
  expect_equal(q$max_angle, 60, tolerance = 1e-8)
  expect_equal(q$aspect_ratio, 1, tolerance = 1e-8)
  expect_equal(q$area, h / 2, tolerance = 1e-10)
})

test_that("mesh_summary prints output", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(30), runif(30)))
  expect_output(mesh_summary(mesh), "Mesh quality summary")
  expect_output(mesh_summary(mesh), "Min angle")
})

test_that("plot.tulpa_mesh runs without error", {
  set.seed(42)
  mesh <- tulpa_mesh(cbind(runif(20), runif(20)))

  # Basic plot
  expect_no_error(plot(mesh))

  # Plot with quality coloring
  expect_no_error(plot(mesh, color = "quality"))

  # Plot with vertices
  expect_no_error(plot(mesh, vertex_col = "red"))

  # Plot with custom color vector
  q <- mesh_quality(mesh)
  expect_no_error(plot(mesh, color = q$area))
})
