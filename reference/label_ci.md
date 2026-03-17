# Label a confidence interval matrix with column and row names.

Label a confidence interval matrix with column and row names.

## Usage

``` r
label_ci(ci, theta, alpha)
```

## Arguments

- ci:

  A two-column matrix of lower and upper bounds.

- theta:

  Named numeric vector of parameter estimates.

- alpha:

  The lower-tail probability (e.g., 0.025 for a 95 percent CI).

## Value

The `ci` matrix with column names set to percentile labels and row names
set to parameter names (or `"param1"`, `"param2"`, ...).
