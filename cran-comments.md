## R CMD check results

0 errors | 0 warnings | 0 notes

## Breaking changes in v2.0.0

This is a major version update addressing the CRAN rejection of v1.2.0.

The S3 class was renamed from `"mle"` to `"mle_fit"` to resolve a name
collision with `stats4::mle` (S4 class). When both packages were loaded,
S4 dispatch intercepted calls to `AIC()`, `BIC()`, `logLik()`, `coef()`,
`summary()`, `vcov()`, `confint()`, and `nobs()`, causing errors.

Additional breaking changes:
- Custom generics `aic()`, `bic()`, `loglik_val()` removed in favor of
  standard R generics `AIC()`, `BIC()`, `logLik()`
- `mle_weighted()` removed (replaced by `combine()`)

Constructor function names (`mle()`, `mle_numerical()`, `mle_boot()`) are
unchanged.

## Downstream dependencies

- `likelihood.model` (on CRAN): updated to `algebraic.mle (>= 2.0.0)`
- `compositional.mle` (not on CRAN): updated

## Test environments

* local Ubuntu 24.04, R 4.3.3
