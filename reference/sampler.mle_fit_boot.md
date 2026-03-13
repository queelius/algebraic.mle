# Method for sampling from an \`mle_fit_boot\` object.

It creates a sampler for the \`mle_fit_boot\` object. It returns a
function that accepts a single parameter \`n\` denoting the number of
samples to draw from the \`mle_fit_boot\` object.

## Usage

``` r
# S3 method for class 'mle_fit_boot'
sampler(x, ...)
```

## Arguments

- x:

  the \`mle_fit_boot\` object to create sampler for

- ...:

  additional arguments to pass (not used)

## Value

A function that takes parameter `n` and returns `n` samples drawn from
the bootstrap replicates.

## Details

Unlike the \`sampler\` method for the more general \`mle_fit\` objects,
for \`mle_fit_boot\` objects, we sample from the bootstrap replicates,
which are more representative of the sampling distribution, particularly
for small samples.

## Examples

``` r
set.seed(1)
b <- boot::boot(rexp(50, 2), function(d, i) 1/mean(d[i]), R = 99)
fit <- mle_boot(b)
sampler(fit)(5)
#> [1] 1.869070 2.257127 2.497161 2.295957 1.925246
```
