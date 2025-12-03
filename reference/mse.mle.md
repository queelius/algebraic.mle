# Computes the MSE of an \`mle\` object.

The MSE of an estimator is just the expected sum of squared differences,
e.g., if the true parameter value is \`x\` and we have an estimator
\`x.hat\`, then the MSE is “\` mse(x.hat) = E\[(x.hat-x) vcov(x.hat) +
bias(x.hat, x) “\`

## Usage

``` r
# S3 method for class 'mle'
mse(x, theta = NULL)
```

## Arguments

- x:

  the \`mle\` object to compute the MSE of.

- theta:

  true parameter value, defaults to \`NULL\` for unknown. If \`NULL\`,
  then we let the bias method deal with it. Maybe it has a nice way of
  estimating the bias.

## Details

Since \`x\` is not typically known, we normally must estimate the bias.
Asymptotically, assuming the regularity conditions, the bias of an MLE
is zero, so we can estimate the MSE as \`mse(x.hat) = vcov(x.hat)\`, but
for small samples, this is not generally the case. If we can estimate
the bias, then we can replace the bias with an estimate of the bias.

Sometimes, we can estimate the bias analytically, but if not, we can use
something like the bootstrap. For example, if we have a sample of size
\`n\`, we can bootstrap the bias by sampling \`n\` observations with
replacement, computing the MLE, and then computing the difference
between the bootstrapped MLE and the MLE. We can repeat this process
\`B\` times, and then average the differences to get an estimate of the
bias.
