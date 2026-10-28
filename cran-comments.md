## Changes in 2.0.3

`combine()` now aligns named parameter estimates and information matrices by
name, preventing estimates of different quantities from being silently pooled
when parameter orders differ. Incompatible or ambiguous parameter names are
rejected; fully unnamed inputs retain positional matching. Empty collections
now produce an informative error.

`rmap()` now forwards transformation arguments to the point estimate as well
as the variance calculation. Regression tests cover both delta and Monte Carlo
methods, and parameter alignment with correlated covariance matrices.

## Test environment

* Ubuntu 24.04.3 LTS, R 4.3.3, x86_64-pc-linux-gnu
* Local source build and R CMD check --as-cran, 2026-10-28

## R CMD check results

0 errors | 0 warnings | 2 notes

* Checking for future file timestamps: unable to verify current time.
  This is an environment/network time-verification limitation.
* HTML manual validation was skipped because `tidy` is not installed and
  math rendering checks were skipped because the check tool `V8` is unavailable
  in this environment.

Package installation, examples (including --run-donttest), tests, vignette
rebuilding, and the PDF manual all passed.
