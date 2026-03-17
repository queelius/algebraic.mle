# Empirical sampling distribution of the MLE

NOTE: If you store the results of a Monte-carlo simulation in this
wrapper object, then many of the functions defined for this object be
misleading. For example, \`confint.mle_emp\` will return a probability
interval. This is actually a still useful metric to have, since we can
determine the probability that an estimate will contain the true
parameter value, but this is not the , which is typically not what we
want. However, it is good for other things, like computing the bias and
mean squared error.

## Usage

``` r
mle_emp(mles, nobs = NULL)
```

## Arguments

- mles:

  the MLEs

- nobs:

  sample size used to compute each MLE

## Value

an \`mle_empirical\` object
