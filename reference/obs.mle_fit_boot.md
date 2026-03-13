# Method for obtaining the observations used by the \`mle_fit_boot\`.

Method for obtaining the observations used by the \`mle_fit_boot\`.

## Usage

``` r
# S3 method for class 'mle_fit_boot'
obs(x)
```

## Arguments

- x:

  the \`mle_fit_boot\` object to obtain the number of observations for

## Value

The original data used for bootstrapping.

## Examples

``` r
set.seed(1)
b <- boot::boot(rexp(50, 2), function(d, i) 1/mean(d[i]), R = 99)
fit <- mle_boot(b)
head(obs(fit))
#> [1] 0.37759092 0.59082139 0.07285336 0.06989763 0.21803431 1.44748427
```
