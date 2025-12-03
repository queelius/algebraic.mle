# Bootstrapped MLE

Sometimes, the large sample asymptotic theory of MLEs is not applicable.
In such cases, we can use the bootstrap to estimate the sampling
distribution of the MLE.

## Usage

``` r
mle_boot(x)
```

## Arguments

- x:

  the \`boot\` return value

## Value

an \`mle_boot\` object (wrapper for \`boot\` object)

## Details

This takes an approach similiar to the \`mle_numerical\` object, which
is a wrapper for a \`stats::optim\` return value, or something that is
compatible with the \`optim\` return value. Here, we take a \`boot\`
object, which is the sampling distribution of an MLE, and wrap it in an
\`mle_boot\` object and then provide a number of methods for the
\`mle_boot\` object that satisfies the concept of an \`mle\` object.

Look up the \`boot\` package for more information on the bootstrap.
