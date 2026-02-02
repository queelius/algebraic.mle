# Computes the estimate of the MSE of a \`boot\` object.

Computes the estimate of the MSE of a \`boot\` object.

## Usage

``` r
# S3 method for class 'mle_boot'
mse(x, theta = NULL, ...)
```

## Arguments

- x:

  the \`boot\` object to compute the MSE of.

- theta:

  true parameter value (not used for \`mle_boot\`)

- ...:

  pass additional arguments into \`vcov\`

## Value

The MSE matrix estimated from bootstrap variance and bias.
