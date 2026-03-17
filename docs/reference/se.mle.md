# Function for obtaining an estimate of the standard error of the MLE object \`x\`.

Function for obtaining an estimate of the standard error of the MLE
object \`x\`.

## Usage

``` r
# S3 method for class 'mle'
se(x, se.matrix = FALSE, ...)
```

## Arguments

- x:

  the MLE object

- se.matrix:

  if \`TRUE\`, return the square root of the variance-covariance

- ...:

  additional arguments to pass (not used)

## Value

Vector of standard errors, or matrix if `se.matrix = TRUE`, or NULL if
variance-covariance is not available.
