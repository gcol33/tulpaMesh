# Get or Set the CRS of a Mesh

Access or assign a coordinate reference system to a `tulpa_mesh` or
`tulpa_mesh_graph` object. The CRS is stored as metadata and propagated
through mesh operations.

## Usage

``` r
mesh_crs(x)

set_crs(x, value)
```

## Arguments

- x:

  A `tulpa_mesh`, `tulpa_mesh_graph`, or `tulpa_mesh_1d` object.

- value:

  A CRS specification: an integer EPSG code, a PROJ string, a WKT
  string, an
  [`sf::st_crs()`](https://r-spatial.github.io/sf/reference/st_crs.html)
  object, or `NULL` to remove.

## Value

`mesh_crs()` returns the CRS (an
[`sf::crs`](https://r-spatial.github.io/sf/reference/coerce-methods.html)
object or `NULL`). `set_crs()` returns the mesh with CRS attached.
