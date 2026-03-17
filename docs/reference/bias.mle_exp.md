# Computes the Bias of an Exponential MLE

This function computes the bias of the Maximum Likelihood Estimator
(MLE) of the rate parameter for an Exponential distribution. The bias is
\`1/(nobs(x)-1)\*par\` where \`x\` is an \`mle_exp\` object and \`par\`
is the true rate parameter. If \`par\` is not provided, the point
estimate from \`x\` is used to estimate the bias.

## Usage

``` r
# S3 method for mle_exp
bias(x, par = NULL, ...)
```

## Arguments

- x:

  an \`mle_exp\` object from which to compute the bias.

- par:

  a numeric value representing the true rate parameter. If NULL, the
  point estimate from \`x\` is used.

- ...:

  additional arguments (not used).

## Value

Returns a named numeric value representing the estimated bias of the
rate parameter.

## Examples

``` r
if (FALSE) {
  obs <- rexp(100, rate = 0.5)
  rate.hat <- mle_exp(obs)
  print(bias(rate.hat, par = 0.5))
}
```
