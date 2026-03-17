# MLE of the rate parameter when we assume the sample is i.i.d. and drawn from the weibull distribution.

MLE of the rate parameter when we assume the sample is i.i.d. and drawn
from the weibull distribution.

## Usage

``` r
mle_weibull_shape_scale(x, k0 = 1, eps = 1e-07, keep_obs = F)
```

## Arguments

- x:

  a sample of observations

- eps:

  we numerically solve the MLE equation, \`\|old-new\| \<= eps\` is
  stopping condition

- keep_obs:

  Boolean, specifies whether to keep observations

- k:

  initial estimate of shape parameter k

## Value

an \`mle\` object.
