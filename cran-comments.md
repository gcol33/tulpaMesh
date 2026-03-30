## R CMD check results

0 errors | 0 warnings | 0 notes

* Resubmission addressing reviewer feedback (Benjamin Altmann, 2026-03-30).

## Changes since last submission

* Removed single quotes from person names and algorithm names in Title and
  Description fields. Only package/software/API names remain quoted.
* Added Artem Amirkhanov (CDT library, MPL-2.0) and William C. Lenthe
  (predicates.h, BSD-3-Clause) to Authors@R with ctb and cph roles.
* Updated inst/COPYRIGHTS to separately document William C. Lenthe's
  geometric predicates code (src/cdt/predicates.h).

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
