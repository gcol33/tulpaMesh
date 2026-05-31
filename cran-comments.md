## R CMD check results

0 errors | 0 warnings | 0 notes

* Patch release fixing two latent bugs in mesh refinement.

## Changes since last submission (0.1.1 -> 0.1.2)

* Fixed `tulpa_mesh(max_edge = ...)` collapsing to zero triangles for some
  `(max_edge, cutoff)` settings and point counts. Vertex deduplication could
  merge extended-hull boundary vertices that fell within `max_edge * 0.3` of
  each other, deleting a constraint edge and leaving the boundary loop open;
  the constrained triangulation then erased every triangle. Boundary and hole
  vertices are now protected from deduplication so the constraint loop stays
  closed.
* `tulpa_mesh()` now errors loudly when the triangulation yields zero
  triangles (for example a collinear point set) instead of returning an empty
  mesh that silently produces all-zero FEM matrices.

## Test environments

* local: Windows 11 Pro, R 4.5.2
* win-builder: R-devel (pending)
* mac-builder: R-release (pending)

## Downstream dependencies

No reverse dependencies.

## Notes

* The package vendors the CDT (Constrained Delaunay Triangulation) C++
  header-only library by Artem Amirkhanov under MPL-2.0 license, documented
  in inst/COPYRIGHTS. The CDT library includes predicates.h by William C.
  Lenthe under BSD-3-Clause. Both authors are listed in Authors@R with
  ctb + cph roles.
* Uses RcppParallel for optional parallel FEM assembly; configure/configure.win
  scripts handle TBB linking portably.
