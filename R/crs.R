#' Get or Set the CRS of a Mesh
#'
#' Access or assign a coordinate reference system to a `tulpa_mesh`
#' or `tulpa_mesh_graph` object. The CRS is stored as metadata and
#' propagated through mesh operations.
#'
#' @param x A `tulpa_mesh`, `tulpa_mesh_graph`, or `tulpa_mesh_1d` object.
#' @param value A CRS specification: an integer EPSG code, a PROJ string,
#'   a WKT string, an `sf::st_crs()` object, or `NULL` to remove.
#'
#' @return `mesh_crs()` returns the CRS (an `sf::crs` object or `NULL`).
#'   `set_crs()` returns the mesh with CRS attached.
#'
#' @export
mesh_crs <- function(x) {
  x$crs
}

#' @rdname mesh_crs
#' @export
set_crs <- function(x, value) {
  if (is.null(value)) {
    x$crs <- NULL
    return(x)
  }

  if (requireNamespace("sf", quietly = TRUE)) {
    x$crs <- sf::st_crs(value)
  } else {
    # Store as-is if sf not available
    x$crs <- value
  }
  x
}
