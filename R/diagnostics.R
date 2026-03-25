#' Per-Triangle Mesh Quality Metrics
#'
#' Computes quality metrics for each triangle in a mesh: minimum angle,
#' maximum angle, aspect ratio, and area.
#'
#' @param mesh A `tulpa_mesh` object.
#'
#' @return A data.frame with one row per triangle and columns:
#'   \itemize{
#'     \item `min_angle`: minimum interior angle (degrees)
#'     \item `max_angle`: maximum interior angle (degrees)
#'     \item `aspect_ratio`: ratio of circumradius to twice the inradius
#'       (1 for equilateral, larger for worse quality)
#'     \item `area`: triangle area
#'   }
#'
#' @export
mesh_quality <- function(mesh) {
  if (!inherits(mesh, "tulpa_mesh")) {
    stop("mesh must be a tulpa_mesh object")
  }

  v <- mesh$vertices
  tri <- mesh$triangles
  n_tri <- nrow(tri)

  min_angle <- numeric(n_tri)
  max_angle <- numeric(n_tri)
  aspect_ratio <- numeric(n_tri)
  area <- numeric(n_tri)

  for (t in seq_len(n_tri)) {
    p0 <- v[tri[t, 1], ]
    p1 <- v[tri[t, 2], ]
    p2 <- v[tri[t, 3], ]

    # Edge vectors
    e0 <- p1 - p0
    e1 <- p2 - p1
    e2 <- p0 - p2

    # Edge lengths
    l0 <- sqrt(sum(e0^2))
    l1 <- sqrt(sum(e1^2))
    l2 <- sqrt(sum(e2^2))

    # Angles via dot product (angle at vertex i is between edges meeting there)
    a0 <- acos(pmin(1, pmax(-1, sum(e0 * (-e2)) / (l0 * l2))))
    a1 <- acos(pmin(1, pmax(-1, sum((-e0) * e1) / (l0 * l1))))
    a2 <- pi - a0 - a1

    angles <- c(a0, a1, a2) * 180 / pi
    min_angle[t] <- min(angles)
    max_angle[t] <- max(angles)

    # Area via cross product
    a <- abs(e0[1] * (-e2[2]) - e0[2] * (-e2[1])) / 2
    area[t] <- a

    # Aspect ratio: circumradius / (2 * inradius)
    # circumradius R = (l0 * l1 * l2) / (4 * area)
    # inradius r = 2 * area / (l0 + l1 + l2)
    # ratio = R / (2r) = (l0 * l1 * l2 * (l0 + l1 + l2)) / (16 * area^2)
    if (a > 0) {
      aspect_ratio[t] <- (l0 * l1 * l2 * (l0 + l1 + l2)) / (16 * a^2)
    } else {
      aspect_ratio[t] <- Inf
    }
  }

  data.frame(
    min_angle = min_angle,
    max_angle = max_angle,
    aspect_ratio = aspect_ratio,
    area = area
  )
}

#' Print Mesh Quality Summary
#'
#' Prints min/median/max of mesh quality metrics.
#'
#' @param mesh A `tulpa_mesh` object.
#'
#' @return Invisible data.frame of per-triangle quality metrics.
#'
#' @export
mesh_summary <- function(mesh) {
  q <- mesh_quality(mesh)

  cat("Mesh quality summary (", nrow(q), " triangles):\n", sep = "")
  cat("  Min angle:    ",
      sprintf("min=%.1f  median=%.1f  max=%.1f deg",
              min(q$min_angle), median(q$min_angle), max(q$min_angle)), "\n")
  cat("  Max angle:    ",
      sprintf("min=%.1f  median=%.1f  max=%.1f deg",
              min(q$max_angle), median(q$max_angle), max(q$max_angle)), "\n")
  cat("  Aspect ratio: ",
      sprintf("min=%.2f  median=%.2f  max=%.2f",
              min(q$aspect_ratio), median(q$aspect_ratio), max(q$aspect_ratio)),
      "  (1 = equilateral)\n")
  cat("  Area:         ",
      sprintf("min=%.2e  median=%.2e  max=%.2e",
              min(q$area), median(q$area), max(q$area)), "\n")

  n_bad <- sum(q$min_angle < 20)
  if (n_bad > 0) {
    cat("  Warning: ", n_bad, " triangles with min angle < 20 deg (",
        sprintf("%.1f%%", 100 * n_bad / nrow(q)), ")\n", sep = "")
  }

  invisible(q)
}

#' Plot a Triangular Mesh
#'
#' Draws the mesh using base R graphics: edges as line segments,
#' vertices as points. Optionally colors triangles by a quality metric.
#'
#' @param x A `tulpa_mesh` object.
#' @param color Optional per-triangle numeric vector to color triangles
#'   (e.g., output of `mesh_quality()$min_angle`). If `"quality"`, uses
#'   minimum angle. If `NULL`, draws edges only.
#' @param border Edge color. Default `"grey50"`.
#' @param vertex_col Vertex point color. Default `NULL` (no vertices drawn).
#' @param vertex_cex Vertex point size. Default 0.5.
#' @param palette Color palette function for triangle fill. Default
#'   `grDevices::hcl.colors`.
#' @param n_colors Number of colors in palette. Default 100.
#' @param main Plot title.
#' @param ... Additional arguments passed to [plot.default()].
#'
#' @return The `tulpa_mesh` object `x`, returned invisibly. Called for
#'   the side effect of producing a plot.
#'
#' @export
plot.tulpa_mesh <- function(x, color = NULL, border = "grey50",
                            vertex_col = NULL, vertex_cex = 0.5,
                            palette = grDevices::hcl.colors,
                            n_colors = 100L, main = "tulpa_mesh", ...) {
  v <- x$vertices
  tri <- x$triangles

  # Set up plot area
  plot(v, type = "n", asp = 1, xlab = "x", ylab = "y", main = main, ...)

  # Resolve color
  if (identical(color, "quality")) {
    color <- mesh_quality(x)$min_angle
  }

  if (!is.null(color)) {
    pal <- palette(n_colors)
    rng <- range(color, finite = TRUE)
    if (rng[1] == rng[2]) {
      col_idx <- rep(1L, length(color))
    } else {
      col_idx <- as.integer((color - rng[1]) / (rng[2] - rng[1]) * (n_colors - 1)) + 1L
    }

    for (t in seq_len(nrow(tri))) {
      idx <- tri[t, ]
      graphics::polygon(v[idx, 1], v[idx, 2],
                         col = pal[col_idx[t]], border = border)
    }
  } else {
    # Draw edges only
    for (t in seq_len(nrow(tri))) {
      idx <- c(tri[t, ], tri[t, 1])
      graphics::lines(v[idx, 1], v[idx, 2], col = border)
    }
  }

  # Draw vertices
  if (!is.null(vertex_col)) {
    graphics::points(v[, 1], v[, 2], pch = 16, col = vertex_col, cex = vertex_cex)
  }

  invisible(x)
}
