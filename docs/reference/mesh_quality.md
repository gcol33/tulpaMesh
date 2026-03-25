# Per-Triangle Mesh Quality Metrics

Computes quality metrics for each triangle in a mesh: minimum angle,
maximum angle, aspect ratio, and area.

## Usage

``` r
mesh_quality(mesh)
```

## Arguments

- mesh:

  A `tulpa_mesh` object.

## Value

A data.frame with one row per triangle and columns:

- `min_angle`: minimum interior angle (degrees)

- `max_angle`: maximum interior angle (degrees)

- `aspect_ratio`: ratio of circumradius to twice the inradius (1 for
  equilateral, larger for worse quality)

- `area`: triangle area
