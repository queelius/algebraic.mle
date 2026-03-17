# Function for obtaining an empirical sampling distribution from the asymptotic distribution of the `mle` object `x` when applied to a function `g`.

Function for obtaining an empirical sampling distribution from the
asymptotic distribution of the `mle` object `x` when applied to a
function `g`.

## Usage

``` r
pred.interval(x, g, n = 1000, alpha = 0.05, ...)
```

## Arguments

- x:

  the `mle` object

- g:

  the function to predict the value of with respect to input `x`

- n:

  the sample size

- alpha:

  confidence-level, (1-alpha)

- ...:

  pass additional arguments
