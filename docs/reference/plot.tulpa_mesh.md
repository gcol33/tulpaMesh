# Plot a Triangular Mesh

Draws the mesh using base R graphics: edges as line segments, vertices
as points. Optionally colors triangles by a quality metric.

## Usage

``` r
# S3 method for class 'tulpa_mesh'
plot(
  x,
  color = NULL,
  border = "grey50",
  vertex_col = NULL,
  vertex_cex = 0.5,
  palette = grDevices::hcl.colors,
  n_colors = 100L,
  main = "tulpa_mesh",
  ...
)
```

## Arguments

- x:

  A `tulpa_mesh` object.

- color:

  Optional per-triangle numeric vector to color triangles (e.g., output
  of `mesh_quality()$min_angle`). If `"quality"`, uses minimum angle. If
  `NULL`, draws edges only.

- border:

  Edge color. Default `"grey50"`.

- vertex_col:

  Vertex point color. Default `NULL` (no vertices drawn).

- vertex_cex:

  Vertex point size. Default 0.5.

- palette:

  Color palette function for triangle fill. Default
  [`grDevices::hcl.colors`](https://rdrr.io/r/grDevices/palettes.html).

- n_colors:

  Number of colors in palette. Default 100.

- main:

  Plot title.

- ...:

  Additional arguments passed to
  [`plot.default()`](https://rdrr.io/r/graphics/plot.default.html).

## Value

The `tulpa_mesh` object `x`, returned invisibly. Called for the side
effect of producing a plot.
