#' Create a Triangular Mesh for SPDE Spatial Models
#'
#' Generates a constrained Delaunay triangulation from point coordinates,
#' optionally with boundary constraints. The resulting mesh can be used
#' directly with tulpa's SPDE spatial fields.
#'
#' @param coords A matrix or data.frame with columns x and y, or a formula
#'   like `~ x + y` evaluated in `data`.
#' @param data Optional data.frame for formula evaluation.
#' @param boundary Optional boundary specification: a matrix of boundary
#'   vertex coordinates (N x 2), an sf polygon, or NULL for convex hull.
#' @param max_edge Maximum edge length. A single value or a vector of two
#'   values `c(inner, outer)` where `inner` controls the study area and
#'   `outer` controls the extension region.
#' @param cutoff Minimum distance between mesh vertices. Points closer than
#'   this are merged. Default 0 (no merging).
#' @param extend Numeric extension factor beyond the boundary. Default 0.1
#'   (10% of domain diameter). Set to 0 for no extension.
#' @param min_angle Minimum angle (degrees) for Ruppert refinement. If
#'   specified, Steiner points are inserted at circumcenters of triangles
#'   with angles below this threshold. Theoretical maximum is ~20.7 degrees;
#'   values up to 30 usually work. Default `NULL` (no refinement).
#' @param max_area Maximum triangle area for refinement. Triangles larger
#'   than this are refined regardless of angle quality. Default `NULL`
#'   (no area constraint).
#' @param max_steiner Maximum number of Steiner points to insert during
#'   Ruppert refinement. Default 10000.
#'
#' @return A `tulpa_mesh` object with components:
#'   \itemize{
#'     \item `vertices`: N x 2 matrix of vertex coordinates
#'     \item `triangles`: M x 3 integer matrix of vertex indices (1-based)
#'     \item `edges`: K x 2 integer matrix of edge vertex indices (1-based)
#'     \item `n_vertices`: number of mesh vertices
#'     \item `n_triangles`: number of triangles
#'     \item `n_edges`: number of edges
#'     \item `n_input_points`: number of original input points
#'   }
#'
#' @export
#'
#' @examples
#' # Simple mesh from random points
#' set.seed(42)
#' coords <- cbind(runif(50), runif(50))
#' mesh <- tulpa_mesh(coords)
#' print(mesh)
#'
#' # Mesh with Ruppert refinement (min angle 25 degrees)
#' mesh_refined <- tulpa_mesh(coords, min_angle = 25)
tulpa_mesh <- function(coords, data = NULL, boundary = NULL,
                       max_edge = NULL, cutoff = 0, extend = 0.1,
                       min_angle = NULL, max_area = NULL,
                       max_steiner = 10000L) {
  # Resolve coordinates
  if (inherits(coords, "formula")) {
    if (is.null(data)) stop("data must be provided when coords is a formula")
    vars <- all.vars(coords)
    if (length(vars) != 2) stop("formula must have exactly 2 variables (x and y)")
    coord_mat <- cbind(data[[vars[1]]], data[[vars[2]]])
  } else if (is.data.frame(coords)) {
    coord_mat <- as.matrix(coords[, 1:2])
  } else {
    coord_mat <- as.matrix(coords)
  }

  if (ncol(coord_mat) != 2) stop("coords must have exactly 2 columns")
  storage.mode(coord_mat) <- "double"

  # Extract CRS from sf boundary if available
  input_crs <- NULL
  if (!is.null(boundary) && inherits(boundary, c("sf", "sfc")) &&
      requireNamespace("sf", quietly = TRUE)) {
    input_crs <- sf::st_crs(boundary)
    if (is.na(input_crs)) input_crs <- NULL
  }

  # Build boundary edges if provided
  boundary_edges <- NULL
  all_points <- coord_mat

  if (!is.null(boundary)) {
    holes <- NULL
    if (inherits(boundary, c("sf", "sfc", "sfg"))) {
      sf_result <- sf_to_boundary(boundary)
      boundary <- sf_result$boundary
      holes <- sf_result$holes
    }
    boundary <- as.matrix(boundary)
    if (ncol(boundary) != 2) stop("boundary must have 2 columns (x, y)")

    # Add boundary vertices to point set
    n_orig <- nrow(coord_mat)
    n_bnd <- nrow(boundary)
    all_points <- rbind(coord_mat, boundary)

    # Create boundary edge loop (connecting consecutive boundary vertices)
    bnd_start <- n_orig + 1L
    boundary_edges <- make_edge_loop(bnd_start, n_bnd)

    # Add hole constraint edges
    if (length(holes) > 0) {
      for (hole in holes) {
        hole <- as.matrix(hole)
        n_hole <- nrow(hole)
        hole_start <- nrow(all_points) + 1L
        all_points <- rbind(all_points, hole)
        boundary_edges <- rbind(boundary_edges,
                                make_edge_loop(hole_start, n_hole))
      }
    }
  } else if (extend > 0) {
    # Create convex hull boundary with extension
    hull_idx <- grDevices::chull(coord_mat)
    hull_pts <- coord_mat[hull_idx, , drop = FALSE]

    # Extend outward
    center <- colMeans(hull_pts)
    hull_extended <- t(apply(hull_pts, 1, function(p) {
      center + (1 + extend) * (p - center)
    }))

    n_orig <- nrow(coord_mat)
    n_bnd <- nrow(hull_extended)
    all_points <- rbind(coord_mat, hull_extended)

    bnd_start <- n_orig + 1L
    boundary_edges <- make_edge_loop(bnd_start, n_bnd)
  }

  # Add refinement points if max_edge is specified
  if (!is.null(max_edge)) {
    me <- max_edge[1]
    xr <- range(all_points[, 1])
    yr <- range(all_points[, 2])
    # Hexagonal lattice: produces more equilateral triangles than rectangular
    grid_pts <- hex_lattice(xr[1], xr[2], yr[1], yr[2], spacing = me)
    all_points <- rbind(all_points, grid_pts)
  }

  # Remove near-duplicate points (always deduplicate when max_edge is used,
  # since hex lattice can produce exact duplicates with boundary vertices)
  dedup_tol <- if (!is.null(max_edge) && cutoff == 0) 1e-10 else cutoff
  if (dedup_tol > 0) {
    keep <- rep(TRUE, nrow(all_points))
    for (i in 2:nrow(all_points)) {
      if (!keep[i]) next
      dists <- sqrt(rowSums((all_points[1:(i-1), , drop = FALSE] -
                               matrix(all_points[i, ], nrow = i-1, ncol = 2, byrow = TRUE))^2))
      if (any(dists[keep[1:(i-1)]] < dedup_tol)) keep[i] <- FALSE
    }
    # Remap boundary edges if points were removed
    if (!is.null(boundary_edges) && any(!keep)) {
      remap <- cumsum(keep)
      remap[!keep] <- 0L
      boundary_edges[, 1] <- remap[boundary_edges[, 1]]
      boundary_edges[, 2] <- remap[boundary_edges[, 2]]
      # Remove edges referencing removed points
      valid <- boundary_edges[, 1] > 0 & boundary_edges[, 2] > 0
      boundary_edges <- boundary_edges[valid, , drop = FALSE]
    }
    all_points <- all_points[keep, , drop = FALSE]
  }

  # Triangulate
  use_ruppert <- !is.null(min_angle) || !is.null(max_area)

  if (use_ruppert) {
    result <- cpp_ruppert_refine(
      points = all_points,
      edges_nullable = boundary_edges,
      min_angle_deg = if (!is.null(min_angle)) min_angle else 0,
      max_area = if (!is.null(max_area)) max_area else 0,
      max_steiner = as.integer(max_steiner),
      min_dist_tolerance = max(cutoff, 1e-10)
    )
  } else {
    result <- cpp_triangulate(
      points = all_points,
      edges_nullable = boundary_edges,
      min_dist_tolerance = max(cutoff, 1e-10)
    )
  }

  result <- structure(result, class = "tulpa_mesh")
  if (!is.null(input_crs)) result$crs <- input_crs
  result
}

#' Compute FEM Matrices from a Mesh
#'
#' Computes the finite element mass (C) and stiffness (G) matrices from
#' a triangular mesh, plus the projection matrix (A) that maps mesh
#' vertices to observation locations.
#'
#' @param mesh A `tulpa_mesh` object.
#' @param obs_coords Observation coordinates (N x 2 matrix). If NULL, the
#'   projection matrix A is the identity (observations at mesh nodes).
#' @param barrier Optional logical vector of length `n_triangles`. Triangles
#'   marked `TRUE` are treated as physical barriers (coastlines, rivers):
#'   their stiffness contributions are zeroed so the spatial field cannot
#'   smooth across them. Based on Bakka et al. (2019).
#' @param parallel Logical. If `TRUE`, uses parallel FEM assembly via
#'   RcppParallel (thread-local triplet accumulation). Beneficial for
#'   meshes with >50K triangles. Default `FALSE`.
#' @param lumped Logical. If `TRUE`, returns a diagonal lumped mass matrix
#'   C0 (vertex areas) in addition to the consistent mass matrix C. The
#'   lumped mass inverse is trivial and needed for the SPDE Q-builder.
#'   Default `FALSE`.
#'
#' @return A list with sparse matrices (dgCMatrix class from Matrix package):
#'   \itemize{
#'     \item `C`: consistent mass matrix (n_vertices x n_vertices)
#'     \item `G`: stiffness matrix (n_vertices x n_vertices)
#'     \item `A`: projection matrix (n_obs x n_vertices)
#'     \item `n_mesh`: number of mesh vertices
#'     \item `C0`: (only if `lumped = TRUE`) diagonal lumped mass matrix
#'     \item `va`: (only if `lumped = TRUE`) vertex areas (numeric vector)
#'     \item `ta`: (only if `lumped = TRUE`) triangle areas (numeric vector)
#'   }
#'
#' @export
fem_matrices <- function(mesh, obs_coords = NULL, barrier = NULL,
                        parallel = FALSE, lumped = FALSE) {
  if (!inherits(mesh, "tulpa_mesh")) {
    stop("mesh must be a tulpa_mesh object")
  }

  # Compute C and G — dispatch based on mesh type and parallel flag
  is_3d <- ncol(mesh$vertices) == 3
  if (is_3d) {
    fem <- cpp_fem_matrices_3d(mesh$vertices, mesh$triangles)
  } else if (parallel) {
    fem <- cpp_fem_matrices_parallel(mesh$vertices, mesh$triangles)
  } else {
    fem <- cpp_fem_matrices(mesh$vertices, mesh$triangles)
  }

  # Assemble into sparse matrices (0-based triplets to 1-based for Matrix)
  n <- fem$n
  C <- Matrix::sparseMatrix(
    i = fem$C_i + 1L, j = fem$C_j + 1L, x = fem$C_x,
    dims = c(n, n), repr = "C"
  )
  G <- Matrix::sparseMatrix(
    i = fem$G_i + 1L, j = fem$G_j + 1L, x = fem$G_x,
    dims = c(n, n), repr = "C"
  )

  # Projection matrix A
  if (!is.null(obs_coords)) {
    obs_coords <- as.matrix(obs_coords)
    storage.mode(obs_coords) <- "double"
    if (is_3d) {
      # For spherical meshes, convert lon/lat obs_coords to xyz if needed
      if (ncol(obs_coords) == 2) {
        radius <- if (!is.null(mesh$sphere)) mesh$sphere$radius else 1
        lon_rad <- obs_coords[, 1] * pi / 180
        lat_rad <- obs_coords[, 2] * pi / 180
        obs_coords <- cbind(
          radius * cos(lat_rad) * cos(lon_rad),
          radius * cos(lat_rad) * sin(lon_rad),
          radius * sin(lat_rad)
        )
      }
      proj <- cpp_projection_matrix_3d(obs_coords, mesh$vertices,
                                        mesh$triangles)
    } else {
      proj <- cpp_projection_matrix(obs_coords, mesh$vertices, mesh$triangles)
    }
    A <- Matrix::sparseMatrix(
      i = proj$i + 1L, j = proj$j + 1L, x = proj$x,
      dims = c(proj$nrow, proj$ncol), repr = "C"
    )
  } else {
    A <- Matrix::Diagonal(n)
  }

  # Barrier model: zero out stiffness contributions from barrier triangles
  # This prevents the spatial field from smoothing across barriers.
  # See Bakka et al. (2019) for the barrier model formulation.
  if (!is.null(barrier)) {
    barrier <- as.logical(barrier)
    if (length(barrier) != nrow(mesh$triangles)) {
      stop(sprintf(
        "barrier must be a logical vector of length n_triangles (%d), got %d",
        nrow(mesh$triangles), length(barrier)
      ))
    }
    # Rebuild G excluding barrier triangles
    # Re-run C++ FEM but with a triangle mask
    fem_barrier <- cpp_fem_matrices(mesh$vertices, mesh$triangles)
    # Identify which triplet entries come from barrier triangles
    # Each triangle contributes exactly 9 entries to both C and G
    n_tri <- nrow(mesh$triangles)
    tri_mask <- rep(!barrier, each = 9)

    # Rebuild G with barrier triangles zeroed out
    G <- Matrix::sparseMatrix(
      i = fem_barrier$G_i[tri_mask] + 1L,
      j = fem_barrier$G_j[tri_mask] + 1L,
      x = fem_barrier$G_x[tri_mask],
      dims = c(n, n), repr = "C"
    )
  }

  result <- list(C = C, G = G, A = A, n_mesh = n)

  if (lumped) {
    va <- as.numeric(Matrix::rowSums(C))
    result$C0 <- Matrix::Diagonal(x = va)
    result$va <- va

    # Triangle areas
    v <- mesh$vertices
    tri <- mesh$triangles
    ta <- numeric(nrow(tri))
    for (t in seq_len(nrow(tri))) {
      p0 <- v[tri[t, 1], ]
      p1 <- v[tri[t, 2], ]
      p2 <- v[tri[t, 3], ]
      ta[t] <- abs((p1[1] - p0[1]) * (p2[2] - p0[2]) -
                    (p2[1] - p0[1]) * (p1[2] - p0[2])) / 2
    }
    result$ta <- ta
  }

  result
}

# Internal: convert sf geometry to boundary coords + hole constraints.
# Returns a list with:
#   $boundary: N x 2 matrix of outer boundary coordinates
#   $holes: list of M x 2 matrices of hole boundary coordinates (may be empty)
sf_to_boundary <- function(geom) {
  if (!requireNamespace("sf", quietly = TRUE)) {
    stop("sf package required for polygon boundaries")
  }

  # Normalize to a single sfc_POLYGON or sfc_MULTIPOLYGON
  geom <- sf::st_geometry(geom)
  if (length(geom) > 1) {
    geom <- sf::st_union(geom)
  }
  geom <- geom[[1]]  # extract the raw sfg

  # Collect all polygon rings
  if (inherits(geom, "MULTIPOLYGON")) {
    # MULTIPOLYGON: list of polygons, each a list of rings
    polygons <- lapply(geom, identity)
  } else if (inherits(geom, "POLYGON")) {
    polygons <- list(geom)
  } else {
    stop(sprintf("unsupported geometry type: %s", class(geom)[1]))
  }

  # Extract outer rings and holes from all component polygons
  outer_rings <- list()
  hole_rings <- list()

  for (poly in polygons) {
    rings <- lapply(poly, identity)  # each ring is a matrix
    # First ring is outer boundary
    outer_rings <- c(outer_rings, list(close_ring(rings[[1]])))
    # Remaining rings are holes
    if (length(rings) > 1) {
      for (k in 2:length(rings)) {
        hole_rings <- c(hole_rings, list(close_ring(rings[[k]])))
      }
    }
  }

  # For the boundary, combine all outer rings
  # (for MULTIPOLYGON with separate components, the user should union first
  #  or we take the largest component as outer boundary)
  if (length(outer_rings) == 1) {
    boundary <- outer_rings[[1]]
  } else {
    # Use the ring with the largest area as the outer boundary,
    # treat other outer rings as additional constraints
    areas <- vapply(outer_rings, ring_area, numeric(1))
    main_idx <- which.max(areas)
    boundary <- outer_rings[[main_idx]]
    # Other outer rings become "island" constraints (treated like holes
    # but with reversed winding, so they remain filled)
    for (idx in seq_along(outer_rings)) {
      if (idx != main_idx) {
        hole_rings <- c(hole_rings, list(outer_rings[[idx]]))
      }
    }
  }

  list(boundary = boundary, holes = hole_rings)
}

# Remove closing point from a ring matrix (last row == first row)
close_ring <- function(ring) {
  ring <- ring[, 1:2, drop = FALSE]
  if (nrow(ring) > 1 && all(ring[1, ] == ring[nrow(ring), ])) {
    ring <- ring[-nrow(ring), , drop = FALSE]
  }
  ring
}

# Shoelace area of a ring (absolute value)
ring_area <- function(ring) {
  n <- nrow(ring)
  x <- ring[, 1]
  y <- ring[, 2]
  abs(sum(x * c(y[-1], y[1]) - c(x[-1], x[1]) * y)) / 2
}

# Build an edge loop connecting consecutive 1-based vertex indices
# start_idx: 1-based index of first vertex in the loop
# n: number of vertices in the loop
make_edge_loop <- function(start_idx, n) {
  edges <- matrix(0L, nrow = n, ncol = 2)
  for (i in seq_len(n)) {
    edges[i, 1] <- start_idx + i - 1L
    edges[i, 2] <- start_idx + (i %% n)
  }
  edges
}

# Generate hexagonal lattice points within a bounding box.
# Row spacing = spacing * sqrt(3)/2; odd rows offset by spacing/2.
# Produces more equilateral triangles than a rectangular grid.
hex_lattice <- function(x_min, x_max, y_min, y_max, spacing) {
  row_height <- spacing * sqrt(3) / 2
  ys <- seq(y_min, y_max, by = row_height)
  pts <- vector("list", length(ys))
  for (i in seq_along(ys)) {
    offset <- if (i %% 2 == 0) spacing / 2 else 0
    xs <- seq(x_min + offset, x_max, by = spacing)
    pts[[i]] <- cbind(xs, rep(ys[i], length(xs)))
  }
  do.call(rbind, pts)
}
