# Linearly transform an `mle` object `x` by multiplying it on the LHS by a matrix `A`.

Linearly transform an `mle` object `x` by multiplying it on the LHS by a
matrix `A`.

## Usage

``` r
mle_linear_transform(A, x, b = NULL)
```

## Arguments

- A:

  a (non-random) matrix

- x:

  an `mle` object to linearly transform

- b:

  a (non-random) vector, defaults to \`NULL\` (zero vector)

## Value

a \`mle_linear_transform\` object
