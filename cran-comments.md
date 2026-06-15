## R CMD check results

0 errors | 0 warnings | 0 notes

## Resubmission

The 0.1.2 submission was archived by the incoming pretest with one NOTE:
"Non-standard file/directory found at top level: '_build.R'". That was a local
build helper that had been swept into the tarball; it is removed and the build
ignore rules now exclude top-level helper scripts.

This version also fixes the additional issue reported for the released 0.1.1 on
the musl (Alpine Linux) flavour, where parallel FEM assembly failed to compile
because it referenced Intel TBB symbols directly. The worker no longer uses any
TBB-specific API, so the package builds on the TinyThread backend that
`RcppParallel` uses where TBB is unavailable.

## Changes since last CRAN release (0.1.1 -> 0.1.3)

* Parallel FEM assembly no longer references Intel TBB symbols directly; it
  builds on the `RcppParallel` TinyThread backend (musl) as well as TBB.
* Fixed `tulpa_mesh(max_edge = ...)` collapsing to zero triangles for some
  `(max_edge, cutoff)` settings: boundary and hole vertices are now protected
  from deduplication so the constraint loop stays closed.
* `tulpa_mesh()` now errors when the triangulation yields zero triangles
  instead of returning an empty mesh with all-zero FEM matrices.

## Test environments

* local: Windows 11 Pro, R 4.5.2
* win-builder: R-devel (pending)

## Downstream dependencies

No reverse dependencies.

## Notes

* The package vendors the CDT (Constrained Delaunay Triangulation) C++
  header-only library by Artem Amirkhanov under MPL-2.0 license, documented
  in inst/COPYRIGHTS. The CDT library includes predicates.h by William C.
  Lenthe under BSD-3-Clause. Both authors are listed in Authors@R with
  ctb + cph roles.
* Uses RcppParallel for optional parallel FEM assembly; configure/configure.win
  scripts handle TBB linking portably via RcppParallel::CxxFlags()/LdFlags().
