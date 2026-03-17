# `rmap.mle` computes the distribution of `f(x)` where `x` is an `mle` object.

By the invariance property of the MLE, if `x` is an `mle` object, then
under the right conditions, asymptotically, `f(x)` is normally
distributed with an approximation given by
`f(x) ~ normal(f(point(x)),sigma)` where `sigma` is the
variance-covariance of `f(x)`, e.g., the sample covariance of
`f(x1),...f(xn)` where `xj` is sampled from `sampler(x)`.

## Usage

``` r
mle_apply(x, g, method = "delta", keep_obs = F)
```

## Arguments

- x:

  an `mle` object.

- g:

  a function that accepts objects like `point(x)`.

## Details

delta method bootstrap method monte carlo method
