# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`algebraic.mle` is an R package (v0.9.0) that provides an algebra over Maximum Likelihood Estimators (MLEs). It enables manipulation, combination, and extraction of statistical properties from fitted models using R's S3 object-oriented system.

## Development Commands

```r
# Load package for development/testing
devtools::load_all()

# Run package checks (tests, documentation validation)
devtools::check()

# Generate documentation from roxygen2 comments
devtools::document()

# Build documentation website
pkgdown::build_site()

# Render a specific vignette
rmarkdown::render("vignettes/statistics.Rmd")

# Install from GitHub
devtools::install_github("queelius/algebraic.mle")
```

## Architecture

### S3 Class Hierarchy

The package uses R's S3 OOP model with the following class hierarchy:

```
mle (base class)
├── mle_numerical  - Wrapper for stats::optim() results
├── mle_boot       - Bootstrap-based MLEs using boot::boot()
└── mle_weighted   - Combined/weighted MLEs from multiple sources
```

### Source Code Organization (R/)

- **algebraic.mle.R**: Package initialization and S3 generic function definitions
- **mle.R**: Core `mle` class constructor and methods (params, confint, vcov, se, aic, summary, etc.)
- **mle_numerical.R**: Wraps `optim()` output, extracts Hessian for variance-covariance
- **mle_boot.R**: Bootstrap-based inference without asymptotic assumptions
- **mle_weighted.R**: Combines MLEs using inverse Fisher Information weighting
- **mle_rmap.R**: Reparameterization mapping functions
- **utils.R**: Utility functions

### Key MLE Object Fields

- `theta.hat`: Parameter estimates
- `sigma`: Variance-covariance matrix
- `info`: Fisher Information matrix
- `loglike`: Log-likelihood function
- `score`: Score function (gradient)
- `obs`: Observation data
- `nobs`: Number of observations

### Dependencies

Requires companion package `algebraic.dist` (install via `devtools::install_github("queelius/algebraic.dist")`).

## Documentation

- roxygen2 comments generate man pages (NAMESPACE auto-generated)
- pkgdown builds website to `docs/` (deployed via GitHub Actions)
- Three vignettes in `vignettes/`: dgp.Rmd, fitting-common-dist.Rmd, statistics.Rmd
