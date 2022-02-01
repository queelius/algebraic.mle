Algebraic maximum likelihood estimators
================

-   [R package: `algebraic.mle`](#r-package-algebraicmle)
    -   [Installation](#installation)
    -   [API](#api)

<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package: `algebraic.mle`

<!-- badges: start -->
<!-- badges: end -->

An algebra over maximum likelihood estimators (MLE).

MLEs have many desirable, well-defined statistical properties. We define
an algebra over MLEs.

## Installation

You can install the development version of `algebraic.mle` from
[GitHub](https://github.com/) with:

``` r
install.packages("devtools")
devtools::install_github("queelius/algebraic.mle")
```

## API

The object representing a fitted model is a type of `mle` object, the
maximum likelihood estimator of the model with respect to observed data.

We define several generic functions with default implementations for
objects that inherit from `mle`:

1.  `vcov(mle)` returns the variance-covariance matrix of the model’s
    parameter estimates.

2.  `point(mle)` returns the point that maximizes the likelihood of the
    model.

3.  `confint(mle)` returns confidence intervals of the point estimates
    of the `mle` object.

4.  `sample(mle)` maps to a function that may be used to sample from the
    sampling distribution of `mle` object.

5.  `mse(mle)` computes the asymptotic mean squared error of the `mle`
    object.

6.  `mle_weighted(mles)` computes the *weighted* maximum likelihood
    estimate from independent `mle` objects in the list argument `mles`.

7.  `fisher_info(mle)` returns the Fisher information matrix of the
    model’s parameters.

8.  `rmap(mle,g)` provides an approximation of the maximum likelihood
    estimator of `g(mle)`.

9.  Finally, since normal distributions are closed under linear
    transformations, then letting *A* be a *p*-by-*q* matrix and *θ̂* be
    a *q* dimensional random vector (MLE for estmiating *θ*), then
    *A**θ̂* is a *p* dimensional multivariate normal with a mean *A**θ*
    and a variance-covariance *A*′vcov (*θ̂*)*A*. Furthermore, by the
    invariance property of the MLE, *A**θ̂* is the MLE of $\\A\\theta$.

    We define a set of generic functions for multiplying `mle` objects
    by (non-random) matrices and vectors to permit these sort of
    operations, with the asymptotic distributions of the operators
    applied to the `mle` objects exactly known if the asymptotic
    distribution of the `mle` objects are exactly known.

    See `rmap` for when the operation is non-linear.
