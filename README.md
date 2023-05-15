
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package: `algebraic.mle`

<!-- badges: start -->

<!-- badges: end -->

An algebra over maximum likelihood estimators (MLE).

MLEs have many desirable, well-defined statistical properties. We define
an algebra over MLEs.

## Installation

You can install `algebraic.mle` from
[GitHub](https://github.com/queelius/algebraic.mle) with:

``` r
install.packages("devtools")
devtools::install_github("queelius/algebraic.mle")
```

## API

The object representing a fitted model is a type of `mle` object, the
maximum likelihood estimator of the model with respect to observed data.

The API mostly consists of generic methods with implementations for
various `mle` type objects. For a full list of functions, see the
[function
reference](https://queelius.github.io/algebraic.mle/reference/index.html)
for `algebraic.mle`.

Here is a simple minimum working example (MWE) of fitting an exponential
model to exponential data.

``` r
library(algebraic.mle)
# Fit an exponential model to exponential DGP (data generating process)
rate.hat <- mle_exp(rexp(100, 1))
summary(rate.hat)
#> Maximum likelihood estimator of type mle_exp is normally distributed.
#> The estimates of the parameters are given by:
#>      rate 
#> 0.9133631 
#> The standard error is  0.09133631 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>           2.5%    97.5%
#> rate 0.7631283 1.063598
#> The MSE of the estimator is  8.511705e-05 .
#> The log-likelihood is  -109.0622 .
#> The AIC is  220.1244 .
```

You can see tutorials for more examples of using the package in the
[vignettes](https://queelius.github.io/algebraic.mle/articles/index.html).
