# Adaptively Refine a Mesh Based on Error Indicators

Refines triangles where error indicators exceed a threshold by inserting
their centroids as new vertices and re-triangulating. Designed for
iterative solve-refine-re-solve workflows where error indicators come
from an SPDE solver (e.g., tulpa's posterior variance).

## Usage

``` r
refine_mesh(
  mesh,
  indicators,
  threshold = NULL,
  fraction = NULL,
  max_iter = 1L,
  min_area = 0
)
```

## Arguments

- mesh:

  A `tulpa_mesh` object (2D only).

- indicators:

  Numeric vector of per-triangle error indicators (length
  `n_triangles`). Higher values trigger refinement.

- threshold:

  Triangles with indicator above this value are refined. Default:
  `median(indicators)` (refine worst half).

- fraction:

  Alternative to `threshold`: refine this fraction of triangles with the
  highest indicators. Default `NULL` (use threshold).

- max_iter:

  Maximum number of refine-retriangulate iterations. Default 1 (single
  refinement pass).

- min_area:

  Minimum triangle area below which triangles are never refined,
  regardless of indicator. Default 0 (no limit).

## Value

A new `tulpa_mesh` object with refined triangulation.
