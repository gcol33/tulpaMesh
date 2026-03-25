# Contribution Guidelines

First of all, thank you very much for taking the time to contribute
to the **tulpaMesh** project!

This document provides guidelines for contributing to tulpaMesh---its codebase and documentation.
These guidelines are meant to guide you, not to restrict you.
If in doubt, use your best judgment and feel free to propose improvements through an issue or pull request.

#### Table Of Contents

- [Code of Conduct](#code-of-conduct)
- [Installation](#installation)
  - [Obtaining the source](#obtaining-the-source)
  - [Setting up your R environment](#setting-up-your-r-environment)
  - [Installing from source](#installing-from-source)
- [Testing](#testing)
- [Documentation](#documentation)
- [Project organization](#project-organization)
- [Contributing workflow](#contributing-workflow)
- [Style guidelines](#style-guidelines)
- [Pull request checklist](#pull-request-checklist)
- [Reporting bugs](#reporting-bugs)

## Code of Conduct

This project and everyone participating in it is governed by our **Code of Conduct** (`CODE_OF_CONDUCT.md`).
By participating, you are expected to uphold this code and maintain a respectful, inclusive environment.

## Installation

This installation guide is focused on development.
For regular installation, please see the [README](./README.md).

### Obtaining the source

Clone the tulpaMesh repository:

```bash
git clone https://github.com/gcol33/tulpaMesh.git
cd tulpaMesh
```

### Setting up your R environment

tulpaMesh is an R package that uses C++17 code via Rcpp and RcppParallel.

1. **Install required tools**
   - R (>= 4.1)
   - Rtools (Windows) or Xcode Command Line Tools (macOS)
   - Git
   - An editor or IDE (RStudio, VS Code, etc.)

2. **Install development dependencies**
   ```r
   install.packages(c("devtools", "roxygen2", "testthat", "Rcpp",
                       "RcppParallel", "Matrix", "rmarkdown", "knitr"))
   ```

3. **Load the development build**
   ```r
   devtools::load_all()
   ```

### Installing from source

Build and install the package locally:

```r
devtools::install()
```

If you modify C++ code, regenerate Rcpp exports before rebuilding:

```r
Rcpp::compileAttributes()
devtools::document()
devtools::install()
```

## Testing

tulpaMesh uses **testthat** (edition 3) for testing.
All tests are located in `tests/testthat/`.

Run the full test suite:

```r
devtools::test()
```

Run a complete package check:

```r
devtools::check()
```

Run a single test file:

```r
devtools::test(filter = "sphere")
```

Guidelines:
- Keep tests fast and reproducible.
- Use `set.seed()` for random data.
- Include edge cases and expected failures.
- Prefer small examples to large datasets.

## Documentation

Build vignettes:

```r
devtools::build_vignettes()
```

Build the pkgdown site locally:

```r
source("~/.R/build_pkgdown.R")
build_pkgdown_site()
```

## Project organization

```
tulpaMesh/
├── .github/workflows/      <- CI (R-CMD-check, test-coverage)
├── R/                       <- R source files
│   ├── mesh.R               <- tulpa_mesh(), fem_matrices()
│   ├── sphere.R             <- tulpa_mesh_sphere(), print.tulpa_mesh()
│   ├── mesh1d.R             <- tulpa_mesh_1d()
│   ├── graph.R              <- tulpa_mesh_graph()
│   ├── diagnostics.R        <- mesh_quality(), mesh_summary(), plot()
│   ├── barrier.R            <- barrier_triangles()
│   ├── convert.R            <- as_tulpa_mesh()
│   ├── refine.R             <- refine_mesh()
│   ├── mesh_ops.R           <- subdivide, components, subset
│   ├── fem_p2.R             <- fem_matrices_p2()
│   ├── nonstationary.R      <- fem_matrices_nonstationary()
│   └── crs.R                <- mesh_crs(), set_crs()
├── src/                     <- C++ source
│   ├── mesh.cpp             <- All Rcpp exports
│   └── cdt/                 <- Vendored CDT library (MPL-2.0)
├── tests/testthat/          <- Unit tests (444 tests)
├── inst/                    <- CITATION, COPYRIGHTS, WORDLIST
├── man/                     <- Generated documentation (roxygen2)
├── DESCRIPTION
├── NAMESPACE
├── NEWS.md
├── README.md
└── CLAUDE.md                <- Dev instructions for Claude Code
```

## Contributing workflow

1. **Create a feature branch**
   ```bash
   git checkout -b feature/my-feature
   ```
2. **Make focused commits** with clear messages.
3. **Run tests and checks** before committing:
   ```r
   devtools::test()
   devtools::check()
   ```
4. **Update documentation** with `roxygen2` and `NEWS.md`.
5. **Update vignettes/examples** if user-facing behavior changes.
6. **Open a pull request** with a short description of your change.
7. **Respond to review feedback** constructively.

## Style guidelines

### R code

- Use descriptive names and consistent indentation.
- Prefer vectorized operations over loops.
- Validate inputs early with clear error messages.
- Document all exported functions with roxygen2.

### C++ code

- Keep headers minimal and separate interface from implementation.
- Use RAII where possible.
- Comment on algorithmic details or numerical behavior.
- Avoid unnecessary memory allocations.
- Use `Rcpp::stop()`, never `Rf_error()`.

### Tests

- Add or update tests when functionality changes.
- Keep tests minimal and reproducible.
- Avoid external dependencies unless essential.

## Pull request checklist

- [ ] Tests pass (`devtools::test()` and `devtools::check()`)
- [ ] Documentation updated (roxygen + `NEWS.md`)
- [ ] Vignettes/examples updated if needed
- [ ] No unrelated formatting changes
- [ ] PR description clearly explains the change

## Reporting bugs

When reporting an issue, please include:
- A minimal reproducible example (`reprex`)
- Output of `sessionInfo()`
- Expected vs. actual results
- R and operating system version
- Toolchain info if relevant (e.g., Rtools on Windows)

---

By contributing to tulpaMesh, you agree that your code is released under the same license as the package.
