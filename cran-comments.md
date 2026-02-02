## R CMD check results

0 errors | 0 warnings | 1 note

* This is a new submission.

## Downstream dependencies

There are currently no downstream dependencies for this package.

## Notes

This package provides an algebraic framework for working with Maximum
Likelihood Estimators (MLEs). It depends on `algebraic.dist` which is
already on CRAN.

## Resubmission

This is a resubmission addressing reviewer feedback:

1. Added references to DESCRIPTION with ISBN format (Casella & Berger 2002,
   Lehmann & Casella 1998)
2. Added \value tags to all exported function documentation
3. Fixed options() save/restore in vignettes

## Test environments

* local Ubuntu 24.04, R 4.3.3
* win-builder (devel and release)
