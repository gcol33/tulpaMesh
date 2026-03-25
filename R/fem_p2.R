#' Compute P2 (Quadratic) FEM Matrices
#'
#' Computes finite element matrices using 6-node quadratic triangular
#' elements. Each triangle edge gets a midpoint node, giving 6 basis
#' functions per triangle instead of 3. Quadratic elements provide
#' better approximation accuracy with fewer mesh nodes.
#'
#' @param mesh A `tulpa_mesh` object (2D only).
#'
#' @return A list with:
#'   \itemize{
#'     \item `C`: consistent mass matrix (n_total x n_total sparse)
#'     \item `G`: stiffness matrix (n_total x n_total sparse)
#'     \item `n_mesh`: total number of nodes (vertices + midpoints)
#'     \item `n_vertices`: number of original mesh vertices
#'     \item `n_midpoints`: number of added midpoint nodes
#'     \item `vertices`: N x 2 matrix of all node coordinates
#'     \item `triangles6`: M x 6 connectivity matrix
#'       (columns: v0, v1, v2, mid01, mid12, mid20)
#'   }
#'
#' @export
fem_matrices_p2 <- function(mesh) {
  if (!inherits(mesh, "tulpa_mesh")) {
    stop("mesh must be a tulpa_mesh object")
  }
  if (ncol(mesh$vertices) != 2) {
    stop("P2 FEM only supports 2D meshes")
  }

  verts <- mesh$vertices
  tris <- mesh$triangles
  n_verts <- nrow(verts)
  n_tri <- nrow(tris)

  # Create midpoint nodes — one per unique edge
  mid_cache <- new.env(hash = TRUE, parent = emptyenv())
  vert_list <- vector("list", n_verts)
  for (i in seq_len(n_verts)) vert_list[[i]] <- verts[i, ]
  next_idx <- n_verts + 1L

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

  # Build 6-column connectivity
  triangles6 <- matrix(0L, nrow = n_tri, ncol = 6)
  for (t in seq_len(n_tri)) {
    v0 <- tris[t, 1]
    v1 <- tris[t, 2]
    v2 <- tris[t, 3]
    m01 <- get_midpoint(v0, v1)
    m12 <- get_midpoint(v1, v2)
    m20 <- get_midpoint(v2, v0)
    triangles6[t, ] <- c(v0, v1, v2, m01, m12, m20)
  }

  all_verts <- do.call(rbind, vert_list)
  n_total <- nrow(all_verts)
  storage.mode(all_verts) <- "double"

  # Compute P2 FEM matrices
  fem <- cpp_fem_matrices_p2(all_verts, triangles6)

  n <- fem$n
  C <- Matrix::sparseMatrix(
    i = fem$C_i + 1L, j = fem$C_j + 1L, x = fem$C_x,
    dims = c(n, n), repr = "C"
  )
  G <- Matrix::sparseMatrix(
    i = fem$G_i + 1L, j = fem$G_j + 1L, x = fem$G_x,
    dims = c(n, n), repr = "C"
  )

  list(
    C = C, G = G,
    n_mesh = n_total,
    n_vertices = n_verts,
    n_midpoints = n_total - n_verts,
    vertices = all_verts,
    triangles6 = triangles6
  )
}
