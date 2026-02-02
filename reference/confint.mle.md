# Function to compute the confidence intervals of \`mle\` objects.

Function to compute the confidence intervals of \`mle\` objects.

## Usage

``` r
# S3 method for class 'mle'
confint(object, parm = NULL, level = 0.95, use_t_dist = FALSE, ...)
```

## Arguments

- object:

  the \`mle\` object to compute the confidence intervals for

- parm:

  the parameters to compute the confidence intervals for (not used)

- level:

  confidence level, defaults to 0.95 (alpha=.05)

- use_t_dist:

  logical, whether to use the t-distribution to compute the confidence
  intervals.

- ...:

  additional arguments to pass

## Value

Matrix of confidence intervals with columns for lower and upper bounds.
