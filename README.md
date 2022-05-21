Algebraic maximum likelihood estimators
================

-   <a href="#r-package-algebraicmle" id="toc-r-package-algebraicmle">R
    package: <code>algebraic.mle</code></a>
    -   <a href="#installation" id="toc-installation">Installation</a>
    -   <a href="#api" id="toc-api">API</a>

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

For a full list of functions, see the [function
reference](https://queelius.github.io/algebraic.mle/reference/index.html)
for `algebraic.mle`. In what follows, we briefly define the API (generic
functions, mostly) with default implementations for objects that inherit
from `mle`:

1.  `vcov(mle)` returns the variance-covariance matrix of the model’s
    parameter estimates. `mse(mle)` computes the asymptotic mean squared
    error of the `mle` object. `sd(mle)` computes the asymptotic
    standard deviation of the `mle` object. `confint(mle)` returns the
    asymptotic confidence intervals of the components of the parameter
    value estimated by the `mle` object. `point(mle)` returns the point
    that maximizes a likelihood of the model.

2.  `sampler(mle)` maps to a function that may be used to sample from
    the sampling distribution of the `mle` object.

3.  `mle_weighted(mles)` computes the *weighted* maximum likelihood
    estimate from independent `mle` objects in the list argument `mles`.

4.  `fisher_info(mle)` returns the Fisher information matrix of the
    model’s parameters.

5.  `rmap(mle,g)` provides an approximation of the maximum likelihood
    estimator of
    ![g(\theta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g%28%5Ctheta%29 "g(\theta)")
    where `mle` is an estimator of
    ![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta").
    The returned value is an `numerical_mle` estimator of
    ![g(\theta)](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;g%28%5Ctheta%29 "g(\theta)").
    Note that `numerical_mle` inherits from `mle`.

6.  `linear_transform(A,mle)`

    The normal distribution is closed under linear combinations and
    transformations. For instance, let
    ![A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A "A")
    be a
    ![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")-by-![q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q "q")
    matrix and
    ![\hat\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%5Ctheta "\hat\theta")
    be a
    ![q](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;q "q")
    dimensional random vector (MLE of
    ![\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta "\theta")),
    then
    ![A\hat\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%5Chat%5Ctheta "A\hat\theta")
    is a
    ![p](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;p "p")
    dimensional multivariate normal with a mean
    ![A\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%5Ctheta "A\theta")
    and a variance-covariance
    ![\operatorname{trans}(A)\operatorname{vcov}(\hat\theta)A](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Coperatorname%7Btrans%7D%28A%29%5Coperatorname%7Bvcov%7D%28%5Chat%5Ctheta%29A "\operatorname{trans}(A)\operatorname{vcov}(\hat\theta)A").
    By the invariance property of the MLE,
    ![A\hat\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%5Chat%5Ctheta "A\hat\theta")
    is the MLE of
    ![A\theta](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;A%5Ctheta "A\theta").

7.  `+(mle1,mle2)`, `-(mle1,mle2)`, `neg(mle)`, and `pos(mle)`

    Assuming `mle1` and `mle2` are independent `mle` objects, we also
    define addition and subtraction of `mle` objects, e.g.,
    ![\hat\theta_1 + \hat\theta_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Chat%5Ctheta_1%20%2B%20%5Chat%5Ctheta_2 "\hat\theta_1 + \hat\theta_2")
    is asymptotically distributed as

    ![N(\mu = \theta_1+\theta_2,\sigma^2 = \sigma\_{\hat\theta_1}^2 +
        \sigma\_{\hat\theta_2}^2).](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;N%28%5Cmu%20%3D%20%5Ctheta_1%2B%5Ctheta_2%2C%5Csigma%5E2%20%3D%20%5Csigma_%7B%5Chat%5Ctheta_1%7D%5E2%20%2B%0A%20%20%20%20%5Csigma_%7B%5Chat%5Ctheta_2%7D%5E2%29. "N(\mu = \theta_1+\theta_2,\sigma^2 = \sigma_{\hat\theta_1}^2 +
        \sigma_{\hat\theta_2}^2).")

    This linear combination is the MLE of
    ![\theta_1+\theta_2](https://latex.codecogs.com/png.image?%5Cdpi%7B110%7D&space;%5Cbg_white&space;%5Ctheta_1%2B%5Ctheta_2 "\theta_1+\theta_2").
