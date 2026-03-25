# Convert to a tulpa_mesh Object

Generic function to convert mesh objects from other packages into
`tulpa_mesh` objects. Currently supports `fm_mesh_2d` objects from the
fmesher package.

## Usage

``` r
as_tulpa_mesh(x, ...)

# S3 method for class 'fm_mesh_2d'
as_tulpa_mesh(x, ...)

# S3 method for class 'inla.mesh'
as_tulpa_mesh(x, ...)
```

## Arguments

- x:

  Object to convert.

- ...:

  Additional arguments (currently unused).

## Value

A `tulpa_mesh` object.
