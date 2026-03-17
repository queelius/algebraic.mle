# mle_optim

This function takes the output of \`optim\` and turns it into an \`mle\`
object.

## Usage

``` r
mle_optim(sol)
```

## Arguments

- sol:

  the output of \`optim\`

## Value

a \`numerical_mle\` object, specialized for \`optim\` (stats package)
solutions.
