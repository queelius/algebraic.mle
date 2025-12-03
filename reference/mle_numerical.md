# This function takes the output of \`optim\`, \`newton_raphson\`, or \`sim_anneal\` and turns it into an \`mle_numerical\` (subclass of \`mle\`) object.

This function takes the output of \`optim\`, \`newton_raphson\`, or
\`sim_anneal\` and turns it into an \`mle_numerical\` (subclass of
\`mle\`) object.

## Usage

``` r
mle_numerical(sol, options = list(), superclasses = NULL)
```

## Arguments

- sol:

  the output of \`optim\` or \`newton_raphson\`

- options:

  list, options for things like sigma and FIM

- superclasses:

  list, superclasses to add to the \`mle_numerical\` object

## Value

a \`numerical_mle\` object.
