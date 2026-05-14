## Repro for issue #2 — cutoff > 0 always produces 0 triangles
suppressMessages(devtools::load_all())

set.seed(11)
coords <- cbind(runif(200), runif(200))

cases <- list(
  list("default",                            list()),
  list("max_edge=c(0.15,0.4) cutoff=0",      list(max_edge = c(0.15, 0.4))),
  list("max_edge=c(0.15,0.4) cutoff=0.05",   list(max_edge = c(0.15, 0.4), cutoff = 0.05)),
  list("max_edge=c(0.15,0.4) cutoff=0.1",    list(max_edge = c(0.15, 0.4), cutoff = 0.1)),
  list("cutoff=0.05 (no max_edge)",          list(cutoff = 0.05)),
  list("max_edge=0.4 cutoff=0",              list(max_edge = 0.4)),
  list("max_edge=0.4 cutoff=0.05",           list(max_edge = 0.4, cutoff = 0.05))
)

for (case in cases) {
  m <- do.call(tulpa_mesh, c(list(coords = coords), case[[2]]))
  cat(sprintf("[%s] n_vertices=%d n_triangles=%d\n",
              case[[1]], nrow(m$vertices),
              if (is.null(m$triangles)) 0L else nrow(m$triangles)))
}
