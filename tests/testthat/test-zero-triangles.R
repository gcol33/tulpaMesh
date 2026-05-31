# Regression tests for the zero-triangle collapse (issue #3).
#
# tulpa_mesh(max_edge = ...) used to drop to 0 triangles non-monotonically
# across (N, max_edge, cutoff): the dedup pass merged extended-hull boundary
# vertices that fell within max_edge * 0.3 of each other, which deleted a
# constraint edge and left the boundary loop open. eraseOuterTrianglesAndHoles
# then classified every triangle as outer and erased the whole mesh. The
# resulting empty mesh produced an all-zero C and an empty G silently.

test_that("previously-collapsing (N, max_edge, cutoff) cases now triangulate", {
  cases <- list(
    list(n = 300, max_edge = c(0.1, 0.3), cutoff = 0),
    list(n = 300, max_edge = c(0.1, 0.3), cutoff = 0.02),
    list(n = 80,  max_edge = c(0.3, 0.6), cutoff = 0),
    list(n = 80,  max_edge = c(0.1, 0.3), cutoff = 0)
  )
  for (cs in cases) {
    set.seed(1)
    coords <- cbind(runif(cs$n), runif(cs$n))
    mesh <- tulpa_mesh(coords = coords, max_edge = cs$max_edge,
                       cutoff = cs$cutoff)
    expect_gt(mesh$n_triangles, 0L)
  }
})

test_that("max_edge meshes triangulate across N and constraint settings", {
  for (n in c(50, 100, 300)) {
    for (me in list(0.05, 0.1, 0.2, 0.3, c(0.1, 0.3), c(0.3, 0.6))) {
      for (cut in c(0, 0.02)) {
        set.seed(7)
        coords <- cbind(runif(n), runif(n))
        mesh <- tulpa_mesh(coords = coords, max_edge = me, cutoff = cut)
        expect_gt(mesh$n_triangles, 0L)
      }
    }
  }
})

test_that("fem_matrices on a max_edge mesh is non-degenerate", {
  set.seed(1)
  coords <- cbind(runif(300), runif(300))
  mesh <- tulpa_mesh(coords = coords, max_edge = c(0.1, 0.3), cutoff = 0)
  fem <- fem_matrices(mesh, obs_coords = coords)

  # Mass matrix C must have a strictly positive diagonal sum (total area > 0).
  expect_gt(sum(Matrix::diag(fem$C)), 0)
  # Stiffness matrix G must carry entries.
  G <- methods::as(fem$G, "CsparseMatrix")
  expect_gt(length(G@x), 0L)
})

test_that("tulpa_mesh errors loudly on a genuinely degenerate (collinear) input", {
  collinear <- cbind(seq(0, 1, length.out = 20), seq(0, 1, length.out = 20))
  expect_error(tulpa_mesh(collinear), "0 triangles")
  expect_error(tulpa_mesh(collinear, extend = 0), "0 triangles")
})

test_that("fem_matrices errors loudly on a 0-triangle mesh", {
  empty_mesh <- structure(
    list(vertices = cbind(c(0, 1, 2), c(0, 0, 0)),
         triangles = matrix(integer(0), nrow = 0, ncol = 3),
         n_triangles = 0L),
    class = "tulpa_mesh"
  )
  expect_error(fem_matrices(empty_mesh), "0 triangles")
})
