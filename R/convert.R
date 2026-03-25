#' Convert to a tulpa_mesh Object
#'
#' Generic function to convert mesh objects from other packages into
#' `tulpa_mesh` objects. Currently supports `fm_mesh_2d` objects from
#' the fmesher package.
#'
#' @param x Object to convert.
#' @param ... Additional arguments (currently unused).
#'
#' @return A `tulpa_mesh` object.
#'
#' @export
as_tulpa_mesh <- function(x, ...) {
  UseMethod("as_tulpa_mesh")
}

#' @rdname as_tulpa_mesh
#' @export
as_tulpa_mesh.fm_mesh_2d <- function(x, ...) {
  verts <- x$loc[, 1:2, drop = FALSE]
  tris <- x$graph$tv

  # Build edge set from triangles
  n_tri <- nrow(tris)
  edge_set <- list()
  for (t in seq_len(n_tri)) {
    for (j in 1:3) {
      v1 <- tris[t, j]
      v2 <- tris[t, (j %% 3) + 1]
      key <- paste(min(v1, v2), max(v1, v2))
      edge_set[[key]] <- c(min(v1, v2), max(v1, v2))
    }
  }
  edges <- do.call(rbind, edge_set)
  storage.mode(edges) <- "integer"

  result <- structure(
    list(
      vertices = verts,
      triangles = tris,
      edges = edges,
      n_vertices = nrow(verts),
      n_triangles = n_tri,
      n_edges = nrow(edges),
      n_input_points = nrow(verts)
    ),
    class = "tulpa_mesh"
  )

  # Propagate CRS from fmesher mesh
  if (!is.null(x$crs) && requireNamespace("sf", quietly = TRUE)) {
    result$crs <- sf::st_crs(x$crs)
  }

  result
}

#' @rdname as_tulpa_mesh
#' @export
as_tulpa_mesh.inla.mesh <- as_tulpa_mesh.fm_mesh_2d
