#' Identify Barrier Triangles
#'
#' Determines which mesh triangles fall inside barrier regions (e.g.,
#' coastlines, rivers, lakes). A triangle is marked as a barrier if its
#' centroid falls inside any barrier polygon.
#'
#' @param mesh A `tulpa_mesh` object.
#' @param barriers An sf/sfc object defining barrier regions, or a list
#'   of N x 2 coordinate matrices (each defining a closed polygon).
#'
#' @return A logical vector of length `n_triangles`. `TRUE` indicates
#'   the triangle is inside a barrier region.
#'
#' @export
barrier_triangles <- function(mesh, barriers) {
  if (!inherits(mesh, "tulpa_mesh")) {
    stop("mesh must be a tulpa_mesh object")
  }

  # Compute triangle centroids
  v <- mesh$vertices
  tri <- mesh$triangles
  centroids <- cbind(
    (v[tri[, 1], 1] + v[tri[, 2], 1] + v[tri[, 3], 1]) / 3,
    (v[tri[, 1], 2] + v[tri[, 2], 2] + v[tri[, 3], 2]) / 3
  )

  if (inherits(barriers, c("sf", "sfc", "sfg"))) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("sf package required for sf barrier objects")
    }
    barrier_geom <- sf::st_geometry(barriers)
    pts_sf <- sf::st_as_sf(
      data.frame(x = centroids[, 1], y = centroids[, 2]),
      coords = c("x", "y"),
      crs = sf::st_crs(barrier_geom)
    )
    inside <- as.logical(sf::st_intersects(pts_sf, sf::st_union(barrier_geom),
                                            sparse = FALSE))
  } else {
    # List of coordinate matrices — use ray casting
    if (is.matrix(barriers)) {
      barriers <- list(barriers)
    }
    inside <- rep(FALSE, nrow(centroids))
    for (poly in barriers) {
      poly <- as.matrix(poly)
      for (i in seq_len(nrow(centroids))) {
        if (!inside[i]) {
          inside[i] <- point_in_poly(centroids[i, 1], centroids[i, 2], poly)
        }
      }
    }
  }

  inside
}

# Ray casting point-in-polygon test
point_in_poly <- function(px, py, poly) {
  n <- nrow(poly)
  crossings <- 0L
  for (i in seq_len(n)) {
    j <- if (i == n) 1L else i + 1L
    y0 <- poly[i, 2]
    y1 <- poly[j, 2]
    x0 <- poly[i, 1]
    x1 <- poly[j, 1]
    if ((y0 <= py && y1 > py) || (y1 <= py && y0 > py)) {
      t <- (py - y0) / (y1 - y0)
      if (px < x0 + t * (x1 - x0)) {
        crossings <- crossings + 1L
      }
    }
  }
  (crossings %% 2L) == 1L
}
