# Expectation operator applied to \`x\` of type \`mle\` with respect to a function \`g\`. That is, \`E(g(x))\`.

Optionally, we use the CLT to construct a CI(\`alpha\`) for the estimate
of the expectation. That is, we estimate \`E(g(x))\` with the sample
mean and Var(g(x)) with the sigma^2/n, where sigma^2 is the sample
variance of g(x) and n is the number of samples. From these, we
construct the CI.

## Usage

``` r
# S3 method for class 'mle'
expectation(x, g = function(t) t, ..., control = list())
```

## Arguments

- x:

  \`mle\` object

- g:

  characteristic function of interest, defaults to identity

- ...:

  additional arguments to pass to \`g\`

- control:

  a list of control parameters: compute_stats - Whether to compute CIs
  for the expectations, defaults to FALSE n - The number of samples to
  use for the MC estimate, defaults to 10000 alpha - The significance
  level for the confidence interval, defaults to 0.05

## Value

If \`compute_stats\` is FALSE, then the estimate of the expectation,
otherwise a list with the following components: value - The estimate of
the expectation ci - The confidence intervals for each component of the
expectation n - The number of samples
