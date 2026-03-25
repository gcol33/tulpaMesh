#' Non-stationary FEM Matrices with Spatially Varying Parameters
#'
#' Computes weighted FEM matrices for non-stationary SPDE models where
#' the range and variance parameters vary spatially. The weights are
#' per-vertex kappa(s) and tau(s) values, interpolated within each
#' triangle using the element-average approximation.
#'
#' @param mesh A `tulpa_mesh` object (2D only).
#' @param kappa Numeric vector of length `n_vertices` giving the
#'   spatial scale parameter kappa(s) = sqrt(8*nu) / range(s).
#' @param tau Numeric vector of length `n_vertices` giving the
#'   precision scaling tau(s).
#'
#' @return A list with sparse matrices:
#'   \itemize{
#'     \item `Ck`: kappa²-weighted mass matrix
#'     \item `Gk`: kappa²-weighted stiffness matrix
#'     \item `Ct`: tau²-weighted mass matrix
#'     \item `C`: unweighted consistent mass matrix
#'     \item `G`: unweighted stiffness matrix
#'     \item `C0`: lumped (diagonal) mass matrix
#'     \item `n_mesh`: number of mesh vertices
#'   }
#'
#' @export
fem_matrices_nonstationary <- function(mesh, kappa, tau) {
  if (!inherits(mesh, "tulpa_mesh")) {
    stop("mesh must be a tulpa_mesh object")
  }
  if (ncol(mesh$vertices) != 2) {
    stop("non-stationary FEM only supports 2D meshes")
  }

  n <- mesh$n_vertices
  if (length(kappa) != n) {
    stop(sprintf("kappa must have length n_vertices (%d), got %d", n, length(kappa)))
  }
  if (length(tau) != n) {
    stop(sprintf("tau must have length n_vertices (%d), got %d", n, length(tau)))
  }

  # Standard FEM matrices
  fem <- cpp_fem_matrices(mesh$vertices, mesh$triangles)
  C <- Matrix::sparseMatrix(
    i = fem$C_i + 1L, j = fem$C_j + 1L, x = fem$C_x,
    dims = c(n, n), repr = "C"
  )
  G <- Matrix::sparseMatrix(
    i = fem$G_i + 1L, j = fem$G_j + 1L, x = fem$G_x,
    dims = c(n, n), repr = "C"
  )
  C0 <- Matrix::Diagonal(x = as.numeric(Matrix::rowSums(C)))

  # Non-stationary weighted matrices
  ns <- cpp_fem_matrices_nonstationary(
    mesh$vertices, mesh$triangles,
    as.numeric(kappa), as.numeric(tau)
  )

  Ck <- Matrix::sparseMatrix(
    i = ns$Ck_i + 1L, j = ns$Ck_j + 1L, x = ns$Ck_x,
    dims = c(n, n), repr = "C"
  )
  Gk <- Matrix::sparseMatrix(
    i = ns$Gk_i + 1L, j = ns$Gk_j + 1L, x = ns$Gk_x,
    dims = c(n, n), repr = "C"
  )
  Ct <- Matrix::sparseMatrix(
    i = ns$Ct_i + 1L, j = ns$Ct_j + 1L, x = ns$Ct_x,
    dims = c(n, n), repr = "C"
  )

  list(
    Ck = Ck, Gk = Gk, Ct = Ct,
    C = C, G = G, C0 = C0,
    n_mesh = n
  )
}
