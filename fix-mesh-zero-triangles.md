# Bug: `tulpa_mesh(..., max_edge = X)` produces zero triangles for `max_edge >= 0.1`

## Summary

When `tulpa_mesh()` is called on a random 2D point cloud in `[0, 1]^2`
(n = 100) with `max_edge >= 0.1`, the returned mesh has the expected number
of vertices (Steiner points added) but `nrow(mesh$triangles) == 0`. The CDT
itself seems to produce points but never connects them into triangles when
`max_edge` is loose. With `max_edge = 0.05` (tighter) or `max_edge = NULL`
(no refinement constraint), triangulation succeeds normally.

## Reproduction

```r
library(tulpaMesh)
set.seed(42); n <- 100L
coords <- cbind(runif(n), runif(n))

# Loose max_edge => zero triangles
mesh_bad <- tulpa_mesh(coords = coords, max_edge = 0.15, extend = 0.2)
nrow(mesh_bad$triangles)  # => 0
nrow(mesh_bad$vertices)   # => 135  (vertices present, no triangulation)

# Default (no max_edge) works
mesh_def <- tulpa_mesh(coords = coords)
nrow(mesh_def$triangles)  # => 211

# Tight max_edge works
mesh_tight <- tulpa_mesh(coords = coords, max_edge = 0.05)
nrow(mesh_tight$triangles)  # => 991
```

Full sweep (see `dev_notes/mesh_sweep.R` in tulpaGlmm):

| max_edge | extend | vertices | triangles |
|----------|--------|----------|-----------|
| NULL     | any    | 113      | 211       |
| 0.05     | 0.05   | 568      | 991       |
| 0.05     | 0.10   | 630      | 1103      |
| 0.05     | 0.20   | 726      | 1279      |
| 0.10     | any    | ~200     | **0**     |
| 0.15     | any    | ~130     | **0**     |
| 0.20     | any    | ~90      | **0**     |
| 0.30     | any    | ~55      | **0**     |

## Expected behavior

After Steiner-point insertion to satisfy the `max_edge` constraint, the
triangulation should remain a valid CDT (Steiner points are inserted to
*refine* an existing triangulation, not to replace it). Even when the
constraint is so loose that no Steiner points need to be inserted, the
output should be the bare CDT over the input vertices.

## Impact

- **Downstream:** tulpaGlmm's SPDE entry (`tglmm(..., inference =
  "nested_laplace", spatial = list(type = "spde", ...))`) silently inherits
  a zero-stiffness `G` matrix when triangles are missing. The shim's inner
  Newton then diverges (CHOLMOD "matrix not positive definite" across all
  grid points, `n_iter == max_iter`, latent w_mesh ranges in the hundreds,
  `cor(spde_f_means, true_field) ~ 0`).
- **Workaround:** users must pass `max_edge = 0.05` (or some other value
  small enough that the refinement loop succeeds) for SPDE fits to work.
  tulpaGlmm's smoke test (`dev_notes/day-spde-smoke.R`) uses `max_edge =
  0.05` after this discovery.

## Likely location

`src/cdt/` Ruppert refinement path, or the post-refinement triangle
re-extraction. Worth checking whether the "no Steiner points needed"
branch falls through to returning an empty triangle list.

## Suggested test (regression)

```r
test_that("tulpa_mesh always returns a non-empty triangulation when input has >=3 points", {
  set.seed(42)
  coords <- cbind(runif(100), runif(100))
  for (me in c(NA, 0.05, 0.1, 0.15, 0.2, 0.3)) {
    args <- list(coords = coords)
    if (!is.na(me)) args$max_edge <- me
    mesh <- do.call(tulpa_mesh, args)
    expect_gt(nrow(mesh$triangles), 0L,
              info = sprintf("max_edge=%s produced 0 triangles", me))
  }
})
```

## Discovered

2026-05-13, during tulpaGlmm Day-41 SPDE wiring. The SPDE nested-Laplace
entry would fail silently on user-typical hyperparameter choices (`max_edge
~ 0.1`) without this fix.
