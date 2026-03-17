## R CMD check results

0 errors | 0 warnings | 0 notes

## Breaking changes (v0.9.0 → v2.0.2)

This is a major update. The S3 class was renamed from "mle" to "mle_fit" to
resolve a name collision with the S4 class stats4::mle. Custom generics
aic(), bic(), loglik_val() were removed in favor of standard R generics
AIC(), BIC(), logLik().

## Coordinated submission

This is part of a coordinated 6-package submission. All packages are
maintained by me. Updated versions being submitted simultaneously:

- algebraic.dist 1.0.0
- algebraic.mle 2.0.2 (this package)
- likelihood.model 1.0.0
- compositional.mle 2.0.0
- flexhaz 0.5.1
- maskedcauses 0.9.3

## Test environments

* local Ubuntu 24.04, R 4.3.3
