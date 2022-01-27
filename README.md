
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package: `algebraic.mle`

<!-- badges: start -->
<!-- badges: end -->

An algebra over maximum likelihood estimators.

## Installation

You can install the development version of `algebraic.mle` from
[GitHub](https://github.com/) with:

``` r
# install.packages("devtools")
devtools::install_github("queelius/algebraic.mle")
```

## Example: MLE of rate parameter in exponential distribution

Suppose we have a sample of *n* = 50 draws from EXP (*λ*=1).

``` r
n = 50
rate = 1
x <- stats::rexp(n,rate)
head(x)
#> [1] 0.1407042 0.3839293 2.7061341 0.2310358 2.2283157 0.3381703
```

Then, we can estimate *λ* with:

``` r
library(algebraic.mle)
(rate.hat <- algebraic.mle::mle_exp(x))
#> $theta.hat
#> [1] 0.8838971
#> 
#> $info
#>          [,1]
#> [1,] 63.99802
#> 
#> $sigma
#>            [,1]
#> [1,] 0.01562548
#> 
#> $sample_size
#> [1] 50
#> 
#> attr(,"class")
#> [1] "mle_exp"   "mle"       "estimator"
```

## MLE from sum of MLEs

Consider *k* mutually independent MLE estimators of parameter *θ*,
*θ̂*<sub>1</sub>, …, *θ̂*<sub>*k*</sub>, where
*θ̂*<sub>*j*</sub> ∼ *N*(*θ*,*I*<sub>*j*</sub><sup>−1</sup>(*θ*)).

Then, the maximum likelihood estimator of *θ* that incorporates all of
the data in *θ̂*<sub>1</sub>, …, *θ̂*<sub>*k*</sub> is given by the
inverse-variance weighted mean,
*θ̂* = (∑*I*<sub>*j*</sub>(*θ*))<sup>−1</sup>(∑*I*<sub>*j*</sub>(*θ*)*θ*<sub>*j*</sub>).
