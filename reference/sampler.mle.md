# Method for sampling from an \`mle\` object.

It creates a sampler for the \`mle\` object. It returns a function that
accepts a single parameter \`n\` denoting the number of samples to draw
from the \`mle\` object.

## Usage

``` r
# S3 method for class 'mle'
sampler(x, ...)
```

## Arguments

- x:

  the \`mle\` object to create sampler for

- ...:

  additional arguments to pass

## Value

A function that takes parameter `n` and returns `n` samples from the
asymptotic distribution of the MLE.
