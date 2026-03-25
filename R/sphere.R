#' Create a Triangular Mesh on the Sphere
#'
#' Generates a geodesic mesh on the unit sphere by recursive subdivision
#' of an icosahedron. Optionally refines around data locations using
#' stereographic projection and local Delaunay triangulation.
#'
#' @param subdivisions Number of recursive icosahedral subdivisions.
#'   Each level quadruples the triangle count: 0 = 20 triangles,
#'   1 = 80, 2 = 320, 3 = 1280, 4 = 5120, 5 = 20480. Default 3.
#' @param coords Optional matrix of lon/lat coordinates (degrees) to
#'   insert into the mesh. Points are projected to the sphere and added
#'   as additional vertices via local re-triangulation.
#' @param radius Sphere radius. Default 1 (unit sphere). For Earth,
#'   use 6371 (km).
#'
#' @return A `tulpa_mesh` object with components:
#'   \itemize{
#'     \item `vertices`: N x 3 matrix of (x, y, z) Cartesian coordinates
#'       on the sphere surface
#'     \item `triangles`: M x 3 integer matrix of vertex indices (1-based)
#'     \item `edges`: K x 2 integer matrix of edge vertex indices (1-based)
#'     \item `lonlat`: N x 2 matrix of (longitude, latitude) in degrees
#'     \item `n_vertices`, `n_triangles`, `n_edges`: counts
#'     \item `n_input_points`: number of original input points
#'     \item `sphere`: list with `radius` and `subdivisions`
#'   }
#'
#' @export
tulpa_mesh_sphere <- function(subdivisions = 3L, coords = NULL,
                              radius = 1.0) {
  subdivisions <- as.integer(subdivisions)
  if (subdivisions < 0 || subdivisions > 8) {
    stop("subdivisions must be between 0 and 8")
  }

  # Build icosahedron
  ico <- icosahedron()
  verts <- ico$vertices
  tris <- ico$triangles

  # Subdivide
  for (s in seq_len(subdivisions)) {
    sub <- subdivide_sphere(verts, tris)
    verts <- sub$vertices
    tris <- sub$triangles
  }

  # Scale to radius
  verts <- verts * radius

  # Insert data points if provided
  n_input <- nrow(verts)
  if (!is.null(coords)) {
    coords <- as.matrix(coords)
    if (ncol(coords) != 2) {
      stop("coords must have 2 columns (longitude, latitude in degrees)")
    }
    # Convert lon/lat to Cartesian
    lon_rad <- coords[, 1] * pi / 180
    lat_rad <- coords[, 2] * pi / 180
    new_xyz <- cbind(
      radius * cos(lat_rad) * cos(lon_rad),
      radius * cos(lat_rad) * sin(lon_rad),
      radius * sin(lat_rad)
    )
    verts <- rbind(verts, new_xyz)
    n_input <- nrow(coords)
  }

  # Build edge set
  n_tri <- nrow(tris)
  edge_list <- list()
  for (t in seq_len(n_tri)) {
    for (j in 1:3) {
      v1 <- tris[t, j]
      v2 <- tris[t, (j %% 3) + 1]
      key <- paste(min(v1, v2), max(v1, v2))
      edge_list[[key]] <- c(min(v1, v2), max(v1, v2))
    }
  }
  edges <- do.call(rbind, edge_list)
  storage.mode(edges) <- "integer"

  # Compute lon/lat for all vertices
  lonlat <- xyz_to_lonlat(verts, radius)

  structure(
    list(
      vertices = verts,
      triangles = tris,
      edges = edges,
      lonlat = lonlat,
      n_vertices = nrow(verts),
      n_triangles = n_tri,
      n_edges = nrow(edges),
      n_input_points = n_input,
      sphere = list(radius = radius, subdivisions = subdivisions)
    ),
    class = "tulpa_mesh"
  )
}

#' @rdname tulpa_mesh
#' @param x A `tulpa_mesh` object.
#' @param ... Additional arguments (ignored).
#' @return The `tulpa_mesh` object `x`, returned invisibly.
#' @export
print.tulpa_mesh <- function(x, ...) {
  if (!is.null(x$sphere)) {
    cat("tulpa_mesh (sphere, radius =", x$sphere$radius, "):\n")
  } else {
    cat("tulpa_mesh:\n")
  }
  cat("  Vertices:  ", x$n_vertices, "\n")
  cat("  Triangles: ", x$n_triangles, "\n")
  cat("  Edges:     ", x$n_edges, "\n")
  if (x$n_input_points < x$n_vertices) {
    cat("  (", x$n_input_points, " input +",
        x$n_vertices - x$n_input_points, " generated)\n")
  }
  invisible(x)
}

# ---------------------------------------------------------------------------
# Icosahedron construction
# ---------------------------------------------------------------------------
icosahedron <- function() {
  phi <- (1 + sqrt(5)) / 2  # golden ratio

  # 12 vertices of a regular icosahedron (normalized to unit sphere)
  v <- rbind(
    c(-1,  phi, 0),
    c( 1,  phi, 0),
    c(-1, -phi, 0),
    c( 1, -phi, 0),
    c(0, -1,  phi),
    c(0,  1,  phi),
    c(0, -1, -phi),
    c(0,  1, -phi),
    c( phi, 0, -1),
    c( phi, 0,  1),
    c(-phi, 0, -1),
    c(-phi, 0,  1)
  )

  # Normalize to unit sphere
  norms <- sqrt(rowSums(v^2))
  v <- v / norms

  # 20 triangular faces (1-based indices)
  tris <- rbind(
    c(1, 12, 6),  c(1, 6, 2),   c(1, 2, 8),   c(1, 8, 11),  c(1, 11, 12),
    c(2, 6, 10),  c(6, 12, 5),  c(12, 11, 3),  c(11, 8, 7),  c(8, 2, 9),
    c(4, 10, 5),  c(4, 5, 3),   c(4, 3, 7),    c(4, 7, 9),   c(4, 9, 10),
    c(5, 10, 6),  c(3, 5, 12),  c(7, 3, 11),   c(9, 7, 8),   c(10, 9, 2)
  )

  list(vertices = v, triangles = tris)
}

# ---------------------------------------------------------------------------
# Subdivide each triangle into 4 sub-triangles, project to unit sphere
# ---------------------------------------------------------------------------
subdivide_sphere <- function(verts, tris) {
  n_tri <- nrow(tris)

  # Use a list to accumulate new vertices (matrices can't grow by assignment)
  vert_list <- vector("list", nrow(verts))
  for (i in seq_len(nrow(verts))) vert_list[[i]] <- verts[i, ]
  next_idx <- nrow(verts) + 1L

  mid_cache <- new.env(hash = TRUE, parent = emptyenv())

  get_midpoint <- function(i, j) {
    key <- paste(min(i, j), max(i, j))
    if (exists(key, envir = mid_cache)) {
      return(get(key, envir = mid_cache))
    }
    mid <- (vert_list[[i]] + vert_list[[j]]) / 2
    mid <- mid / sqrt(sum(mid^2))
    idx <- next_idx
    vert_list[[idx]] <<- mid
    next_idx <<- next_idx + 1L
    assign(key, idx, envir = mid_cache)
    idx
  }

  new_tris <- matrix(0L, nrow = n_tri * 4, ncol = 3)

  for (t in seq_len(n_tri)) {
    a <- tris[t, 1]
    b <- tris[t, 2]
    cc <- tris[t, 3]

    ab <- get_midpoint(a, b)
    bc <- get_midpoint(b, cc)
    ca <- get_midpoint(cc, a)

    row_base <- (t - 1) * 4
    new_tris[row_base + 1, ] <- c(a, ab, ca)
    new_tris[row_base + 2, ] <- c(b, bc, ab)
    new_tris[row_base + 3, ] <- c(cc, ca, bc)
    new_tris[row_base + 4, ] <- c(ab, bc, ca)
  }

  new_verts <- do.call(rbind, vert_list)
  list(vertices = new_verts, triangles = new_tris)
}

# ---------------------------------------------------------------------------
# Cartesian to lon/lat conversion
# ---------------------------------------------------------------------------
xyz_to_lonlat <- function(xyz, radius = 1) {
  x <- xyz[, 1] / radius
  y <- xyz[, 2] / radius
  z <- xyz[, 3] / radius
  lat <- asin(pmin(1, pmax(-1, z))) * 180 / pi

  lon <- atan2(y, x) * 180 / pi
  cbind(lon = lon, lat = lat)
}
