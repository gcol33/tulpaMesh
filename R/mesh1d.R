#' Create a 1D Mesh for Temporal SPDE Models
#'
#' Generates a 1D mesh (interval partition) for temporal components of
#' spatio-temporal SPDE models. Returns 1D FEM matrices (tridiagonal
#' mass and stiffness) that can be combined with 2D spatial FEM matrices
#' via Kronecker products in tulpa's space-time Q-builder.
#'
#' @param knots Numeric vector of mesh knot locations (e.g., time points).
#'   Will be sorted and deduplicated.
#' @param boundary Two-element numeric vector `c(lower, upper)` defining
#'   the mesh domain. Defaults to `range(knots)`.
#' @param n_extend Number of extra knots to add beyond each boundary,
#'   spaced at the median knot interval. Default 3.
#'
#' @return A `tulpa_mesh_1d` object with components:
#'   \itemize{
#'     \item `knots`: sorted numeric vector of mesh knot locations
#'     \item `n`: number of knots
#'     \item `C`: consistent mass matrix (tridiagonal, n x n sparse)
#'     \item `G`: stiffness matrix (tridiagonal, n x n sparse)
#'     \item `C0`: lumped (diagonal) mass matrix
#'   }
#'
#' @export
tulpa_mesh_1d <- function(knots, boundary = NULL, n_extend = 3L) {
  knots <- sort(unique(as.numeric(knots)))
  if (length(knots) < 2) {
    stop("need at least 2 unique knot locations")
  }

  # Extend beyond boundaries
  if (n_extend > 0) {
    h_med <- median(diff(knots))
    if (is.null(boundary)) {
      boundary <- range(knots)
    }
    ext_lo <- boundary[1] - seq(n_extend, 1) * h_med
    ext_hi <- boundary[2] + seq(1, n_extend) * h_med
    knots <- sort(unique(c(ext_lo, knots, ext_hi)))
  }

  n <- length(knots)
  h <- diff(knots)  # n-1 interval lengths

  # 1D P1 FEM matrices (tridiagonal)
  # Mass matrix C: consistent mass
  #   C[i,i] += h[k]/3 for each adjacent interval k
  #   C[i,i+1] = C[i+1,i] = h[k]/6
  # Stiffness matrix G:
  #   G[i,i] += 1/h[k] for each adjacent interval k
  #   G[i,i+1] = G[i+1,i] = -1/h[k]

  C_i <- integer(0)
  C_j <- integer(0)
  C_x <- numeric(0)
  G_i <- integer(0)
  G_j <- integer(0)
  G_x <- numeric(0)

  for (k in seq_along(h)) {
    i <- k
    j <- k + 1L

    # Mass
    C_i <- c(C_i, i, j, i, j)
    C_j <- c(C_j, i, j, j, i)
    C_x <- c(C_x, h[k] / 3, h[k] / 3, h[k] / 6, h[k] / 6)

    # Stiffness
    G_i <- c(G_i, i, j, i, j)
    G_j <- c(G_j, i, j, j, i)
    G_x <- c(G_x, 1 / h[k], 1 / h[k], -1 / h[k], -1 / h[k])
  }

  C <- Matrix::sparseMatrix(i = C_i, j = C_j, x = C_x,
                             dims = c(n, n), repr = "C")
  G <- Matrix::sparseMatrix(i = G_i, j = G_j, x = G_x,
                             dims = c(n, n), repr = "C")
  C0 <- Matrix::Diagonal(x = as.numeric(Matrix::rowSums(C)))

  structure(
    list(
      knots = knots,
      n = n,
      C = C,
      G = G,
      C0 = C0
    ),
    class = "tulpa_mesh_1d"
  )
}

#' @rdname tulpa_mesh_1d
#' @param x A `tulpa_mesh_1d` object.
#' @param ... Additional arguments (ignored).
#' @return The `tulpa_mesh_1d` object `x`, returned invisibly.
#' @export
print.tulpa_mesh_1d <- function(x, ...) {
  cat("tulpa_mesh_1d:\n")
  cat("  Knots:    ", x$n, "\n")
  cat("  Range:    [", min(x$knots), ", ", max(x$knots), "]\n", sep = "")
  cat("  Intervals:", x$n - 1, "\n")
  invisible(x)
}
