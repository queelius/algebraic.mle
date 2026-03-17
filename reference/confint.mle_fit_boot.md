# Method for obtained the confidence interval of an \`mle_fit_boot\` object. Note: This impelements the \`vcov\` method defined in \`stats\`.

Method for obtained the confidence interval of an \`mle_fit_boot\`
object. Note: This impelements the \`vcov\` method defined in \`stats\`.

## Usage

``` r
# S3 method for class 'mle_fit_boot'
confint(
  object,
  parm = NULL,
  level = 0.95,
  type = c("norm", "basic", "perc", "bca"),
  ...
)
```

## Arguments

- object:

  the \`mle_fit_boot\` object to obtain the confidence interval of

- parm:

  character or integer vector of parameters to return intervals for. If
  NULL (default), all parameters are returned.

- level:

  the confidence level

- type:

  the type of confidence interval to compute

- ...:

  additional arguments to pass into \`boot.ci\`

## Value

Matrix of bootstrap confidence intervals with columns for lower and
upper bounds.

## Examples

``` r
set.seed(1)
b <- boot::boot(rexp(50, 2), function(d, i) 1/mean(d[i]), R = 99)
fit <- mle_boot(b)
confint(fit)
#>            2.5%    97.5%
#> param1 1.445856 2.571918
```
