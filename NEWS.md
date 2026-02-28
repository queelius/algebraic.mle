# algebraic.mle 1.2.0

## New features

* `joint()` composes independent MLEs with disjoint parameter sets into a
  joint MLE with block-diagonal covariance structure
* `combine()` optimally weights independent MLEs for the same parameter via
  inverse-variance (Fisher information) weighting
* `as_dist()` converts MLE objects to their asymptotic normal distributions,
  bridging to the `algebraic.dist` distribution algebra
* Distribution methods on MLE objects: `density()`, `cdf()`, `inv_cdf()`,
  `sup()`, `dim()`, `mean()`, `conditional()`

## Documentation

* Rewrote package Description around "MLE as technology" narrative
* Rewrote README to lead with the algebra (joint, combine, rmap, as_dist)
* Added "The Algebra of MLEs" vignette demonstrating the full composition
  pipeline

## Bug fixes

* Fixed name-propagation issue in `rmap()` where `c()` merged parameter names
  with transformation output names

# algebraic.mle 1.1.0

* Add `coef()` S3 method for base R compatibility (delegates to `params()`)
* Add `logLik()` S3 method returning proper `logLik` object with `df` and `nobs`
  attributes, enabling automatic `AIC()` and `BIC()` support from base R
* Fix `rmap()` to accept numeric `n` parameter (previously required integer)

# algebraic.mle 1.0.0

* Initial CRAN release
* Core MLE class (`mle`) with methods for:
  - Parameter extraction (`params`, `nparams`)
  - Variance-covariance (`vcov`, `se`)
  - Confidence intervals (`confint`)
  - Model comparison (`aic`, `loglik_val`)
  - Bias and MSE estimation (`bias`, `mse`)
  - Fisher information (`observed_fim`)
  - Sampling from MLE distribution (`sampler`)
  - Predictive intervals (`pred`)
  - Expected values (`expectation`)
  - Marginal distributions (`marginal`)
* Numerical optimization wrapper (`mle_numerical`) for `optim()` results
* Bootstrap MLE (`mle_boot`) for small samples
* MLE transformations via invariance property (`rmap`)
* Weighted combination of MLEs (`mle_weighted`)
* Three vignettes demonstrating usage:
  - Fitting common distributions to a DGP
  - Statistics and characteristics of the MLE
  - Data generating processes
