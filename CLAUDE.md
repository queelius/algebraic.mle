# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

`algebraic.mle` is an R package that provides an algebra over Maximum Likelihood Estimators (MLEs). It defines a unified interface for working with MLEs through a type system with a base `mle` class and specialized subclasses (`mle_numerical`, `mle_boot`, `rmap_mle`).

**Dependencies**: Requires `algebraic.dist` package (github::queelius/algebraic.dist)

## Development Commands

### Build and Installation
```r
# Install dependencies first (required)
devtools::install_github("queelius/algebraic.dist")

# Install this package from GitHub
devtools::install_github("queelius/algebraic.mle")

# For development: load package without installing
devtools::load_all()
```

### Documentation
```r
# Generate documentation from roxygen comments
devtools::document()

# Build pkgdown site
pkgdown::build_site()
```

### Testing and Validation
**Note**: This package currently does not have a formal `tests/` directory with unit tests.

Validation is done through:
```r
# Run R CMD check (validates package structure, examples, vignettes)
devtools::check()

# Manually test by running vignette code
rmarkdown::render("vignettes/fitting-common-dist.Rmd")
rmarkdown::render("vignettes/statistics.Rmd")
rmarkdown::render("vignettes/dgp.Rmd")
```

When adding new features, create test cases in vignettes or consider adding a `tests/testthat/` directory structure

### Vignettes
Build vignettes using `knitr`:
```r
# Render specific vignette
rmarkdown::render("vignettes/fitting-common-dist.Rmd")
```

## Core Architecture

### MLE Type System
The package uses an S3 class hierarchy with `mle` as the base type:

- **`mle`** (R/mle.R): Base constructor and common methods for all MLE objects
  - Contains: `theta.hat` (parameters), `loglike`, `score`, `sigma` (vcov), `info` (FIM), `obs`, `nobs`
  - Assumes asymptotic normality: `theta ~ MVN(theta.hat, sigma)`

- **`mle_numerical`** (R/mle_numerical.R): Wrapper for numerical optimization results
  - Wraps `optim()` output into MLE interface
  - Computes variance-covariance from Hessian using `sigma = -ginv(hessian)`
  - Primary constructor for fitting custom models

- **`mle_boot`** (R/mle_boot.R): Bootstrap-based MLE for small samples
  - Wraps `boot` objects when asymptotic theory doesn't apply
  - Estimates bias and variance from bootstrap replicates
  - Methods use empirical distributions rather than asymptotic normality

- **`rmap_mle`** (R/mle_rmap.R): Transformed MLEs via invariance property
  - Implements `g(theta.hat)` for function `g` applied to MLE
  - Two methods: delta method (default) or Monte Carlo
  - Variance propagation: `J %*% vcov(x) %*% t(J)` where J is Jacobian

### Generic Methods
All MLE types support these methods (defined in R/mle.R unless noted):
- `params()`, `nparams()`: Extract parameters
- `vcov()`, `se()`: Variance-covariance and standard errors
- `confint()`: Confidence intervals (asymptotic or bootstrap)
- `loglik_val()`, `aic()`: Model fit statistics
- `bias()`, `mse()`: Bias and mean squared error
- `sampler()`: Generate sampling function from MLE distribution
- `obs()`, `nobs()`: Access original data
- `observed_fim()`: Fisher information matrix
- `score_val()`: Score function at MLE
- `pred()`: Predictive intervals using Monte Carlo
- `expectation()`: Expected values via simulation
- `marginal()`: Marginal distributions for parameter subsets
- `rmap()`: Transformations (invariance property)

### Typical Workflow
1. Define model via `resp()`, parameter functions, and `loglik()`
2. Fit using `optim()` with `control = list(fnscale = -1)` and `hessian = TRUE`
3. Wrap result: `sol <- mle_numerical(optim(...))`
4. Use generic methods: `summary()`, `confint()`, `aic()`, etc.
5. For small samples or non-asymptotic cases, use `mle_boot(boot(...))`

## Important Details

### Numerical Optimization
- Always use `fnscale = -1` in `optim()` control for maximization
- Always request `hessian = TRUE` to compute variance-covariance
- The Hessian is the observed Fisher Information Matrix (negative of Hessian at MLE)
- Variance-covariance: `sigma = -solve(hessian)` (or `-ginv(hessian)` for singular matrices)

### Hypothesis Testing
- Likelihood ratio tests for nested models: `LRT = -2 * (loglik_null - loglik_full)`
- Compare models using `aic()`: lower is better
- Use `confint()` for parameter-level inference

### Bootstrap MLEs
- Use when sample size is small or regularity conditions questionable
- `mle_boot` estimates bias empirically: `mean(t) - t0`
- Bootstrap CIs available: "norm", "basic", "perc", "bca"

### Transformations (rmap)
- By invariance property: if `theta.hat` is MLE, so is `g(theta.hat)`
- Delta method: faster, uses Jacobian for variance propagation
- Monte Carlo: more robust for complex transformations, requires more samples

## Package Structure
- `R/`: Core implementation files
  - `algebraic.mle.R`: Package documentation
  - `mle.R`: Base MLE class and common methods
  - `mle_numerical.R`: Numerical optimization wrapper
  - `mle_boot.R`: Bootstrap MLE implementation
  - `mle_rmap.R`: Transformation/invariance methods
  - `mle_weighted.R`: Weighted likelihood support
  - `utils.R`: Helper functions (confint_from_sigma, etc.)
- `man/`: Auto-generated roxygen documentation (do not edit manually)
- `vignettes/`: Tutorials and examples
  - `fitting-common-dist.Rmd`: Fitting Weibull/Normal to DGP
  - `statistics.Rmd`: Statistical inference with MLEs
  - `dgp.Rmd`: Data generating processes
- `fixing/`: Experimental/WIP code (not part of package build)
  - `mle_weighted_slow.R`: Alternative weighted MLE implementation
  - `hyp_tests.R`: Hypothesis testing experiments
- `docs/`: pkgdown-generated website (auto-generated)

## Common Pitfalls and Best Practices

### When Creating MLE Objects
- **Always** set `fnscale = -1` in `optim()` control (maximization, not minimization)
- **Always** set `hessian = TRUE` in `optim()` to enable variance estimation
- Check convergence: `sol$convergence == 0` before wrapping with `mle_numerical()`
- Use `MASS::ginv()` instead of `solve()` for potentially singular Hessians

### Model Definition Pattern
The standard pattern for fitting custom models:
```r
# 1. Define response extractor
resp <- function(df) df$y

# 2. Define parameter function (relates parameters to distribution)
rate <- function(df, beta) exp(beta[1] + beta[2] * df$x)

# 3. Define log-likelihood as closure
loglik <- function(df, resp, rate) {
  function(beta) sum(dexp(x = resp(df), rate = rate(df, beta), log = TRUE))
}

# 4. Fit with optim, wrap with mle_numerical
par0 <- c(0, 0)  # initial guess
names(par0) <- c("b0", "b1")
sol <- mle_numerical(optim(
  par = par0,
  fn = loglik(df, resp, rate),
  control = list(fnscale = -1),
  hessian = TRUE
))
```

### Weighted Likelihoods
- Use `mle_weighted.R` for case weights or sampling weights
- Weights must be positive and are typically normalized

## Notes
- This package extends `algebraic.dist` and follows its conventions
- All documentation uses roxygen2 format
- The package emphasizes statistical rigor and theoretical foundations
- MLE objects are designed to be composable and algebraically well-defined
- The `fixing/` directory contains experimental code excluded from builds via `.Rbuildignore`
