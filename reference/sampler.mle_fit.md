# Method for sampling from an \`mle_fit\` object.

It creates a sampler for the \`mle_fit\` object. It returns a function
that accepts a single parameter \`n\` denoting the number of samples to
draw from the \`mle_fit\` object.

## Usage

``` r
# S3 method for class 'mle_fit'
sampler(x, ...)
```

## Arguments

- x:

  the \`mle_fit\` object to create sampler for

- ...:

  additional arguments to pass

## Value

A function that takes parameter `n` and returns `n` samples from the
asymptotic distribution of the MLE.
