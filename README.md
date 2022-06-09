
-   <a href="#r-package-algebraicmle" id="toc-r-package-algebraicmle">R
    package: <code>algebraic.mle</code></a>
    -   <a href="#installation" id="toc-installation">Installation</a>
    -   <a href="#todo" id="toc-todo">TODO</a>
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

## TODO

-   Use generics and standard notation for operations on MLEs, like
    addition. May need to look into multiple dispatch.

## API

The object representing a fitted model is a type of `mle` object, the
maximum likelihood estimator of the model with respect to observed data.

The API mostly consists of generic methods with implementations for
various `mle` type objects. For a full list of functions, see the
[function
reference](https://queelius.github.io/algebraic.mle/reference/index.html)
for `algebraic.mle`.
