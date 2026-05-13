## Trace the exact arguments cpp_triangulate gets for the failing case.
library(tulpaMesh)

## Monkeypatch tulpaMesh:::cpp_triangulate to capture inputs
orig <- tulpaMesh:::cpp_triangulate
captured <- list()
wrap <- function(points, edges_nullable = NULL, min_dist_tolerance = 0) {
  captured$points <<- points
  captured$edges <<- edges_nullable
  captured$tol <<- min_dist_tolerance
  res <- orig(points = points, edges_nullable = edges_nullable,
              min_dist_tolerance = min_dist_tolerance)
  captured$result <<- res
  res
}
assignInNamespace("cpp_triangulate", wrap, ns = "tulpaMesh")

set.seed(42)
coords <- cbind(runif(400), runif(400))
m <- tulpa_mesh(coords = coords, max_edge = c(0.3, 0.6))
cat("Final mesh: ", m$n_vertices, "v, ", m$n_triangles, "t\n", sep="")

cat("\ncpp_triangulate was called with:\n")
cat("  points: ", nrow(captured$points), "x", ncol(captured$points), "\n", sep="")
cat("  points range x: [", min(captured$points[,1]), ",", max(captured$points[,1]), "]\n", sep="")
cat("  points range y: [", min(captured$points[,2]), ",", max(captured$points[,2]), "]\n", sep="")
cat("  edges:  ", if (is.null(captured$edges)) "NULL" else paste(nrow(captured$edges), "x 2"), "\n")
if (!is.null(captured$edges)) {
  cat("    edge indices range: [", min(captured$edges), ",", max(captured$edges), "]\n", sep="")
  cat("    first 5 edges:\n")
  print(head(captured$edges, 5))
}
cat("  tol:    ", captured$tol, "\n")

cat("\nResult from cpp_triangulate:\n")
cat("  n_vertices: ", captured$result$n_vertices, "\n")
cat("  n_triangles: ", captured$result$n_triangles, "\n")
cat("  n_edges: ", captured$result$n_edges, "\n")
