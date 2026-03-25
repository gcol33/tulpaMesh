## R CMD check results

0 errors | 0 warnings | 0 notes

* This is a new release.

## Test environments

* local: Windows 11 Pro, R 4.5.2
* win-builder: R-devel (pending)
* mac-builder: R-release (pending)

## Downstream dependencies

No reverse dependencies.

## Notes

* The package vendors the CDT (Constrained Delaunay Triangulation) C++
  header-only library by Artem Amirkhanov under MPL-2.0 license, documented
  in inst/COPYRIGHTS.
* Uses RcppParallel for optional parallel FEM assembly; configure/configure.win
  scripts handle TBB linking portably.
