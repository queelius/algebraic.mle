# Method for obtained the confidence interval of an \`mle_boot\` object. Note: This impelements the \`vcov\` method defined in \`stats\`.

Method for obtained the confidence interval of an \`mle_boot\` object.
Note: This impelements the \`vcov\` method defined in \`stats\`.

## Usage

``` r
# S3 method for class 'mle_boot'
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

  the \`mle_boot\` object to obtain the confidence interval of

- parm:

  the parameter to obtain the confidence interval of (not used)

- level:

  the confidence level

- type:

  the type of confidence interval to compute

- ...:

  additional arguments to pass into \`boot.ci\`

## Value

Matrix of bootstrap confidence intervals with columns for lower and
upper bounds.
