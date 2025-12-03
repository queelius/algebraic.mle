# Computes the variance-covariance matrix of \`boot\` object. Note: This impelements the \`vcov\` method defined in \`stats\`.

Computes the variance-covariance matrix of \`boot\` object. Note: This
impelements the \`vcov\` method defined in \`stats\`.

## Usage

``` r
# S3 method for class 'mle_boot'
vcov(object, ...)
```

## Arguments

- object:

  the \`boot\` object to obtain the variance-covariance of

- ...:

  additional arguments to pass into \`stats::cov\`
