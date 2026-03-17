# MLE of the rate parameter when we assume the sample is i.i.d. and drawn from the weibull distribution.

MLE of the rate parameter when we assume the sample is i.i.d. and drawn
from the weibull distribution.

## Usage

``` r
mle_weibull(
  x,
  par0 = c(1, 1),
  pgtol = 1e-07,
  maxit = 10000L,
  keep_obs = FALSE,
  ...
)
```

## Arguments

- x:

  a sample of observations

- par0:

  initial estimate of the parameters, defaults to \`c(1,1)\`

- pgtol:

  the convergence tolerance for the optimization algorithm

- keep_obs:

  Boolean, specifies whether to keep observations

## Value

an \`mle\` object.
