# Method for sampling from an `boot` object.

It creates a sampler for the `boot` object. It returns a function that
accepts a single parameter `n` denoting the number of samples to draw
from the `boot` object.

## Usage

``` r
# S3 method for boot
sampler(x, ...)
```

## Arguments

- x:

  the `boot` object to create sampler for

- ...:

  additional arguments to pass
