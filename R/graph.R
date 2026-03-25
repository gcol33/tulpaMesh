#' Create a Metric Graph Mesh
#'
#' Builds a 1D FEM mesh along the edges of a spatial network (roads,
#' rivers, coastlines). Each edge is discretized into segments with
#' P1 linear elements. Junction nodes where edges meet are shared.
#' Based on the Whittle-Matern SPDE formulation on metric graphs
#' (Bolin, Simas & Wallin, 2024).
#'
#' @param edges A list of N x 2 coordinate matrices, each defining a
#'   polyline edge of the network. Or an sf object with LINESTRING
#'   geometries.
#' @param max_edge Maximum segment length along edges. Edges longer
#'   than this are subdivided. Default `NULL` (use vertices as-is).
#' @param snap_tolerance Distance below which endpoints are snapped
#'   to the same junction node. Default `1e-8`.
#'
#' @return A `tulpa_mesh_graph` object with components:
#'   \itemize{
#'     \item `vertices`: N x 2 matrix of node coordinates
#'     \item `segments`: M x 2 integer matrix of segment connectivity (1-based)
#'     \item `n_vertices`, `n_segments`: counts
#'     \item `C`: consistent mass matrix (tridiagonal-block sparse)
#'     \item `G`: stiffness matrix (tridiagonal-block sparse)
#'     \item `C0`: lumped (diagonal) mass matrix
#'     \item `degree`: integer vector of node degrees (junction = degree > 2)
#'   }
#'
#' @export
tulpa_mesh_graph <- function(edges, max_edge = NULL, snap_tolerance = 1e-8) {
  # Extract coordinate lists from sf
  if (inherits(edges, c("sf", "sfc"))) {
    if (!requireNamespace("sf", quietly = TRUE)) {
      stop("sf package required for sf input")
    }
    geom <- sf::st_geometry(edges)
    edges <- lapply(geom, function(g) {
      sf::st_coordinates(g)[, 1:2, drop = FALSE]
    })
  }

  if (!is.list(edges)) {
    stop("edges must be a list of coordinate matrices or an sf LINESTRING object")
  }

  # Subdivide long edges
  if (!is.null(max_edge)) {
    edges <- lapply(edges, function(e) {
      subdivide_polyline(as.matrix(e), max_edge)
    })
  }

  # Collect all vertices, snap close endpoints
  all_pts <- do.call(rbind, edges)
  n_total <- nrow(all_pts)

  # Assign unique vertex IDs with snapping
  vertex_id <- integer(n_total)
  unique_verts <- list()
  next_id <- 1L

  for (i in seq_len(n_total)) {
    found <- FALSE
    for (j in seq_along(unique_verts)) {
      d <- sqrt(sum((all_pts[i, ] - unique_verts[[j]])^2))
      if (d < snap_tolerance) {
        vertex_id[i] <- j
        found <- TRUE
        break
      }
    }
    if (!found) {
      unique_verts[[next_id]] <- all_pts[i, ]
      vertex_id[i] <- next_id
      next_id <- next_id + 1L
    }
  }

  verts <- do.call(rbind, unique_verts)
  n_verts <- nrow(verts)

  # Build segments: consecutive pairs within each edge
  segments <- list()
  row_offset <- 0L
  for (e in edges) {
    n_pts <- nrow(e)
    for (k in seq_len(n_pts - 1)) {
      v1 <- vertex_id[row_offset + k]
      v2 <- vertex_id[row_offset + k + 1]
      if (v1 != v2) {  # skip degenerate segments
        segments <- c(segments, list(c(v1, v2)))
      }
    }
    row_offset <- row_offset + n_pts
  }
  seg_mat <- do.call(rbind, segments)
  n_seg <- nrow(seg_mat)

  # Compute node degree
  degree <- integer(n_verts)
  for (s in seq_len(n_seg)) {
    degree[seg_mat[s, 1]] <- degree[seg_mat[s, 1]] + 1L
    degree[seg_mat[s, 2]] <- degree[seg_mat[s, 2]] + 1L
  }

  # Build 1D FEM matrices along graph edges
  C_i <- integer(0); C_j <- integer(0); C_x <- numeric(0)
  G_i <- integer(0); G_j <- integer(0); G_x <- numeric(0)

  for (s in seq_len(n_seg)) {
    v1 <- seg_mat[s, 1]
    v2 <- seg_mat[s, 2]
    h <- sqrt(sum((verts[v1, ] - verts[v2, ])^2))
    if (h < 1e-15) next

    # 1D P1 mass
    C_i <- c(C_i, v1, v2, v1, v2)
    C_j <- c(C_j, v1, v2, v2, v1)
    C_x <- c(C_x, h / 3, h / 3, h / 6, h / 6)

    # 1D P1 stiffness
    G_i <- c(G_i, v1, v2, v1, v2)
    G_j <- c(G_j, v1, v2, v2, v1)
    G_x <- c(G_x, 1 / h, 1 / h, -1 / h, -1 / h)
  }

  C <- Matrix::sparseMatrix(i = C_i, j = C_j, x = C_x,
                             dims = c(n_verts, n_verts), repr = "C")
  G <- Matrix::sparseMatrix(i = G_i, j = G_j, x = G_x,
                             dims = c(n_verts, n_verts), repr = "C")
  C0 <- Matrix::Diagonal(x = as.numeric(Matrix::rowSums(C)))

  structure(
    list(
      vertices = verts,
      segments = seg_mat,
      n_vertices = n_verts,
      n_segments = n_seg,
      C = C, G = G, C0 = C0,
      degree = degree
    ),
    class = "tulpa_mesh_graph"
  )
}

#' @rdname tulpa_mesh_graph
#' @param x A `tulpa_mesh_graph` object.
#' @param ... Additional arguments (ignored).
#' @return The `tulpa_mesh_graph` object `x`, returned invisibly.
#' @export
print.tulpa_mesh_graph <- function(x, ...) {
  cat("tulpa_mesh_graph:\n")
  cat("  Vertices: ", x$n_vertices, "\n")
  cat("  Segments: ", x$n_segments, "\n")
  n_junc <- sum(x$degree > 2)
  n_end <- sum(x$degree == 1)
  cat("  Junctions:", n_junc, "  Endpoints:", n_end, "\n")
  invisible(x)
}

# Subdivide a polyline so no segment exceeds max_length
subdivide_polyline <- function(coords, max_length) {
  result <- list(coords[1, , drop = FALSE])
  for (i in seq_len(nrow(coords) - 1)) {
    p1 <- coords[i, ]
    p2 <- coords[i + 1, ]
    d <- sqrt(sum((p2 - p1)^2))
    if (d > max_length) {
      n_sub <- ceiling(d / max_length)
      for (k in seq_len(n_sub - 1)) {
        frac <- k / n_sub
        result <- c(result, list(matrix(p1 + frac * (p2 - p1), nrow = 1)))
      }
    }
    result <- c(result, list(matrix(p2, nrow = 1)))
  }
  do.call(rbind, result)
}
