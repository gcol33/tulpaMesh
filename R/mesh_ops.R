#' Subdivide a Mesh
#'
#' Splits each triangle into 4 sub-triangles by inserting edge midpoints,
#' then re-triangulates. Useful for multi-resolution workflows where a
#' coarser mesh needs uniform refinement.
#'
#' @param mesh A `tulpa_mesh` object (2D only).
#'
#' @return A new `tulpa_mesh` object with approximately 4x as many triangles.
#'
#' @export
subdivide_mesh <- function(mesh) {
  if (!inherits(mesh, "tulpa_mesh")) {
    stop("mesh must be a tulpa_mesh object")
  }
  if (ncol(mesh$vertices) != 2) {
    stop("subdivide_mesh only supports 2D meshes (use tulpa_mesh_sphere subdivisions for spherical)")
  }

  verts <- mesh$vertices
  tris <- mesh$triangles
  n_tri <- nrow(tris)

  # Midpoint cache
  mid_cache <- new.env(hash = TRUE, parent = emptyenv())
  vert_list <- vector("list", nrow(verts))
  for (i in seq_len(nrow(verts))) vert_list[[i]] <- verts[i, ]
  next_idx <- nrow(verts) + 1L

  get_midpoint <- function(i, j) {
    key <- paste(min(i, j), max(i, j))
    if (exists(key, envir = mid_cache)) {
      return(get(key, envir = mid_cache))
    }
    mid <- (vert_list[[i]] + vert_list[[j]]) / 2
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

  # Build edge set
  edge_list <- list()
  for (t in seq_len(nrow(new_tris))) {
    for (j in 1:3) {
      v1 <- new_tris[t, j]
      v2 <- new_tris[t, (j %% 3) + 1]
      key <- paste(min(v1, v2), max(v1, v2))
      edge_list[[key]] <- c(min(v1, v2), max(v1, v2))
    }
  }
  edges <- do.call(rbind, edge_list)
  storage.mode(edges) <- "integer"

  structure(
    list(
      vertices = new_verts,
      triangles = new_tris,
      edges = edges,
      n_vertices = nrow(new_verts),
      n_triangles = nrow(new_tris),
      n_edges = nrow(edges),
      n_input_points = mesh$n_input_points
    ),
    class = "tulpa_mesh"
  )
}

#' Identify Disconnected Mesh Components
#'
#' Finds connected components of a mesh via triangle adjacency (two
#' triangles are connected if they share an edge).
#'
#' @param mesh A `tulpa_mesh` object.
#'
#' @return An integer vector of length `n_triangles` giving the
#'   component ID (1, 2, ...) for each triangle.
#'
#' @export
mesh_components <- function(mesh) {
  if (!inherits(mesh, "tulpa_mesh")) {
    stop("mesh must be a tulpa_mesh object")
  }

  tris <- mesh$triangles
  n_tri <- nrow(tris)

  # Build adjacency: map each edge to the triangles that share it
  edge_to_tri <- new.env(hash = TRUE, parent = emptyenv())
  for (t in seq_len(n_tri)) {
    for (j in 1:3) {
      v1 <- tris[t, j]
      v2 <- tris[t, (j %% 3) + 1]
      key <- paste(min(v1, v2), max(v1, v2))
      existing <- if (exists(key, envir = edge_to_tri)) {
        get(key, envir = edge_to_tri)
      } else {
        integer(0)
      }
      assign(key, c(existing, t), envir = edge_to_tri)
    }
  }

  # Build triangle adjacency list
  adj <- vector("list", n_tri)
  for (t in seq_len(n_tri)) adj[[t]] <- integer(0)

  for (key in ls(edge_to_tri)) {
    tri_ids <- get(key, envir = edge_to_tri)
    if (length(tri_ids) == 2) {
      adj[[tri_ids[1]]] <- c(adj[[tri_ids[1]]], tri_ids[2])
      adj[[tri_ids[2]]] <- c(adj[[tri_ids[2]]], tri_ids[1])
    }
  }

  # BFS to find components
  component <- integer(n_tri)
  comp_id <- 0L

  for (start in seq_len(n_tri)) {
    if (component[start] > 0) next
    comp_id <- comp_id + 1L
    queue <- start
    component[start] <- comp_id
    while (length(queue) > 0) {
      current <- queue[1]
      queue <- queue[-1]
      for (neighbor in adj[[current]]) {
        if (component[neighbor] == 0) {
          component[neighbor] <- comp_id
          queue <- c(queue, neighbor)
        }
      }
    }
  }

  component
}

#' Extract a Submesh from Triangle Indices
#'
#' Creates a new `tulpa_mesh` containing only the specified triangles,
#' with vertex indices remapped.
#'
#' @param mesh A `tulpa_mesh` object.
#' @param triangles Integer vector of triangle indices to keep, or a
#'   logical vector of length `n_triangles`.
#'
#' @return A new `tulpa_mesh` object containing only the selected
#'   triangles and their vertices.
#'
#' @export
subset_mesh <- function(mesh, triangles) {
  if (!inherits(mesh, "tulpa_mesh")) {
    stop("mesh must be a tulpa_mesh object")
  }

  if (is.logical(triangles)) {
    if (length(triangles) != mesh$n_triangles) {
      stop(sprintf("logical vector must have length n_triangles (%d)",
                    mesh$n_triangles))
    }
    triangles <- which(triangles)
  }

  triangles <- as.integer(triangles)
  if (any(triangles < 1) || any(triangles > mesh$n_triangles)) {
    stop("triangle indices out of range")
  }

  # Extract selected triangles
  old_tris <- mesh$triangles[triangles, , drop = FALSE]

  # Find unique vertices used
  used_verts <- sort(unique(as.integer(old_tris)))
  n_new_verts <- length(used_verts)

  # Remap indices
  remap <- integer(mesh$n_vertices)
  remap[used_verts] <- seq_len(n_new_verts)

  new_verts <- mesh$vertices[used_verts, , drop = FALSE]
  new_tris <- matrix(remap[old_tris], ncol = 3)

  # Build edge set
  edge_list <- list()
  for (t in seq_len(nrow(new_tris))) {
    for (j in 1:3) {
      v1 <- new_tris[t, j]
      v2 <- new_tris[t, (j %% 3) + 1]
      key <- paste(min(v1, v2), max(v1, v2))
      edge_list[[key]] <- c(min(v1, v2), max(v1, v2))
    }
  }
  edges <- do.call(rbind, edge_list)
  storage.mode(edges) <- "integer"

  structure(
    list(
      vertices = new_verts,
      triangles = new_tris,
      edges = edges,
      n_vertices = n_new_verts,
      n_triangles = nrow(new_tris),
      n_edges = nrow(edges),
      n_input_points = n_new_verts
    ),
    class = "tulpa_mesh"
  )
}
