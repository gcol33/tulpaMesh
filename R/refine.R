#' Adaptively Refine a Mesh Based on Error Indicators
#'
#' Refines triangles where error indicators exceed a threshold by inserting
#' their centroids as new vertices and re-triangulating. Designed for
#' iterative solve-refine-re-solve workflows where error indicators come
#' from an SPDE solver (e.g., tulpa's posterior variance).
#'
#' @param mesh A `tulpa_mesh` object (2D only).
#' @param indicators Numeric vector of per-triangle error indicators
#'   (length `n_triangles`). Higher values trigger refinement.
#' @param threshold Triangles with indicator above this value are refined.
#'   Default: `median(indicators)` (refine worst half).
#' @param fraction Alternative to `threshold`: refine this fraction of
#'   triangles with the highest indicators. Default `NULL` (use threshold).
#' @param max_iter Maximum number of refine-retriangulate iterations.
#'   Default 1 (single refinement pass).
#' @param min_area Minimum triangle area below which triangles are never
#'   refined, regardless of indicator. Default 0 (no limit).
#'
#' @return A new `tulpa_mesh` object with refined triangulation.
#'
#' @export
refine_mesh <- function(mesh, indicators, threshold = NULL,
                        fraction = NULL, max_iter = 1L,
                        min_area = 0) {
  if (!inherits(mesh, "tulpa_mesh")) {
    stop("mesh must be a tulpa_mesh object")
  }
  if (ncol(mesh$vertices) != 2) {
    stop("refine_mesh only supports 2D meshes")
  }
  if (length(indicators) != mesh$n_triangles) {
    stop(sprintf(
      "indicators must have length n_triangles (%d), got %d",
      mesh$n_triangles, length(indicators)
    ))
  }

  current_verts <- mesh$vertices
  current_tris <- mesh$triangles

  for (iter in seq_len(max_iter)) {
    n_tri <- nrow(current_tris)

    # Recompute indicators if we've re-triangulated (only first iter uses input)
    if (iter > 1 && length(indicators) != n_tri) {
      # After re-triangulation, indicators don't map to new triangles.
      # Only single-pass refinement is well-defined without a solver callback.
      break
    }

    # Determine which triangles to refine
    if (!is.null(fraction)) {
      n_refine <- max(1L, as.integer(ceiling(fraction * n_tri)))
      cutoff <- sort(indicators, decreasing = TRUE)[n_refine]
      refine <- indicators >= cutoff
    } else {
      if (is.null(threshold)) {
        threshold <- median(indicators)
      }
      refine <- indicators > threshold
    }

    # Don't refine triangles below min_area
    if (min_area > 0) {
      v <- current_verts
      tri <- current_tris
      for (t in which(refine)) {
        p0 <- v[tri[t, 1], ]
        p1 <- v[tri[t, 2], ]
        p2 <- v[tri[t, 3], ]
        area <- abs((p1[1] - p0[1]) * (p2[2] - p0[2]) -
                     (p2[1] - p0[1]) * (p1[2] - p0[2])) / 2
        if (area < min_area) refine[t] <- FALSE
      }
    }

    if (!any(refine)) break

    # Insert centroids of marked triangles
    centroids <- cbind(
      (current_verts[current_tris[refine, 1], 1] +
       current_verts[current_tris[refine, 2], 1] +
       current_verts[current_tris[refine, 3], 1]) / 3,
      (current_verts[current_tris[refine, 1], 2] +
       current_verts[current_tris[refine, 2], 2] +
       current_verts[current_tris[refine, 3], 2]) / 3
    )

    all_points <- rbind(current_verts, centroids)
    storage.mode(all_points) <- "double"

    # Build convex hull boundary to preserve mesh extent
    hull_idx <- grDevices::chull(all_points)
    n_hull <- length(hull_idx)
    hull_edges <- matrix(0L, nrow = n_hull, ncol = 2)
    for (h in seq_len(n_hull)) {
      hull_edges[h, 1] <- hull_idx[h]
      hull_edges[h, 2] <- hull_idx[(h %% n_hull) + 1]
    }

    # Re-triangulate with all points and hull boundary
    result <- cpp_triangulate(
      points = all_points,
      edges_nullable = hull_edges,
      min_dist_tolerance = 1e-10
    )

    current_verts <- result$vertices
    current_tris <- result$triangles
  }

  # Build final mesh object
  n_verts <- nrow(current_verts)
  n_tri <- nrow(current_tris)

  # Rebuild edge set
  edge_list <- list()
  for (t in seq_len(n_tri)) {
    for (j in 1:3) {
      v1 <- current_tris[t, j]
      v2 <- current_tris[t, (j %% 3) + 1]
      key <- paste(min(v1, v2), max(v1, v2))
      edge_list[[key]] <- c(min(v1, v2), max(v1, v2))
    }
  }
  edges <- do.call(rbind, edge_list)
  storage.mode(edges) <- "integer"

  structure(
    list(
      vertices = current_verts,
      triangles = current_tris,
      edges = edges,
      n_vertices = n_verts,
      n_triangles = n_tri,
      n_edges = nrow(edges),
      n_input_points = mesh$n_input_points
    ),
    class = "tulpa_mesh"
  )
}
