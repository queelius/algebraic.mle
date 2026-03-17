# Method for sampling from an \`mle_boot\` object.

It creates a sampler for the \`mle_boot\` object. It returns a function
that accepts a single parameter \`n\` denoting the number of samples to
draw from the \`mle_boot\` object.

## Usage

``` r
# S3 method for class 'mle_boot'
sampler(x, ...)
```

## Arguments

- x:

  the \`mle_boot\` object to create sampler for

- ...:

  additional arguments to pass (not used)

## Value

A function that takes parameter `n` and returns `n` samples drawn from
the bootstrap replicates.

## Details

Unlike the \`sampler\` method for the more general \`mle\` objects, for
\`mle_boot\` objects, we sample from the bootstrap replicates, which are
more representative of the sampling distribution, particularly for small
samples.
