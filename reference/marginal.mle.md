# Method for obtaining the marginal distribution of an MLE that is based on asymptotic assumptions:

\`x ~ MVN(params(x), inv(H)(x))\`

## Usage

``` r
# S3 method for class 'mle'
marginal(x, indices)
```

## Arguments

- x:

  The distribution object.

- indices:

  The indices of the marginal distribution to obtain.

## Details

where H is the (observed or expecation) Fisher information matrix.
