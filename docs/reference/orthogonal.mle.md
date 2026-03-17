# Method for determining the orthogonal components of an \`mle\` object \`x\`.

Method for determining the orthogonal components of an \`mle\` object
\`x\`.

## Usage

``` r
# S3 method for class 'mle'
orthogonal(x, tol = sqrt(.Machine$double.eps), ...)
```

## Arguments

- x:

  the \`mle\` object

- tol:

  the tolerance for determining if a number is close enough to zero

- ...:

  pass additional arguments

## Value

Logical matrix indicating which off-diagonal FIM elements are
approximately zero (orthogonal parameters), or NULL if FIM unavailable.
