# Changelog

## algebraic.mle 1.0.0

- Initial CRAN release
- Core MLE class (`mle`) with methods for:
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
- Numerical optimization wrapper (`mle_numerical`) for
  [`optim()`](https://rdrr.io/r/stats/optim.html) results
- Bootstrap MLE (`mle_boot`) for small samples
- MLE transformations via invariance property (`rmap`)
- Weighted combination of MLEs (`mle_weighted`)
- Three vignettes demonstrating usage:
  - Fitting common distributions to a DGP
  - Statistics and characteristics of the MLE
  - Data generating processes
