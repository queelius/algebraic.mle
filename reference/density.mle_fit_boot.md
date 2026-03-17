# PDF of the empirical distribution of bootstrap replicates.

PDF of the empirical distribution of bootstrap replicates.

## Usage

``` r
# S3 method for class 'mle_fit_boot'
density(x, ...)
```

## Arguments

- x:

  An `mle_fit_boot` object.

- ...:

  Additional arguments (not used).

## Value

A function computing the empirical PMF at given points.

## Examples

``` r
set.seed(1)
b <- boot::boot(rexp(50, 2), function(d, i) 1/mean(d[i]), R = 99)
fit <- mle_boot(b)
density(fit)
#> function (t, log = FALSE) 
#> {
#>     stopifnot(is.numeric(t))
#>     if (is.matrix(t)) {
#>         stopifnot(ncol(t) == p)
#>         counts <- apply(t, 1, function(t_ob) {
#>             sum(apply(xobs, 1, function(xob) {
#>                 all(xob == t_ob)
#>             }))
#>         })
#>     }
#>     else {
#>         stopifnot(length(t) == p)
#>         counts <- sum(apply(xobs, 1, function(xob) {
#>             all(xob == t)
#>         }))
#>     }
#>     if (log) {
#>         return(log(counts) - log(n))
#>     }
#>     else {
#>         return(counts/n)
#>     }
#> }
#> <bytecode: 0x558fe64f0fe8>
#> <environment: 0x558fe64f24c0>
```
