# Computes the score of an \`mle_fit\` object (score evaluated at the MLE).

If reguarlity conditions are satisfied, it should be zero (or
approximately, if rounding errors occur).

## Usage

``` r
# S3 method for class 'mle_fit'
score_val(x, ...)
```

## Arguments

- x:

  the \`mle_fit\` object to compute the score of.

- ...:

  additional arguments to pass (not used)

## Value

The score vector evaluated at the MLE, or NULL if not available.
