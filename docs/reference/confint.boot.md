# Function to compute the confidence intervals of `mle` objects.

Function to compute the confidence intervals of `mle` objects.

## Usage

``` r
# S3 method for boot
confint(object, parm = NULL, level = 0.95, type = "perc", ...)
```

## Arguments

- object:

  the `boot` object to compute the confidence intervals for

- parm:

  parameter indexes to compute the confidence intervals for, defaults to
  all. NOTE: not implemented so ignored for now.

- level:

  confidence level, defaults to 0.95 (alpha=.05)

- type:

  A vector of character strings representing the type of intervals
  required. The value should be any subset of the values
  `c("norm","basic", "stud", "perc", "bca")` or simply "all".

- ...:

  additional arguments to pass
