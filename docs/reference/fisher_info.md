# Generic method for obtaining the fisher information matrix of an `mle` object.

Fisher information is a way of measuring the amount of information that
an observable random variable \`X\` carries about an unknown parameter
`theta` upon which the probability of \`X\` depends.

## Usage

``` r
fisher_info(x, ...)
```

## Arguments

- x:

  the object to obtain the fisher information of

- ...:

  additional arguments to pass

## Details

The inverse of the Fisher information matrix is the variance-covariance
of the MLE for `theta`.
