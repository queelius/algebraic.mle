# mle_local_search

This assumes the MLE is an interior point and that you have provided an
initial guess \`theta0\` that is near it. Use a global search method
like \`sim_anneal\` to find a good initial guess.

## Usage

``` r
mle_local_search(ll, dir, theta0, options = list())
```

## Arguments

- ll:

  function, log-likelihood function

- dir:

  function, promising direction function

- theta0:

  numeric, initial guess

- options:

  list, options for the local search, see function description.

## Value

an \`mle\` object with additional attributes \`iter\` and \`converged\`
and optionally \`path\` if \`trace\` is TRUE.

## Functions

- `mle_local_search()`: Optional Arguments

## Fields

- `sup`:

  function, domain of support for log-likelihood

- `eta`:

  numeric, learning rate, defaults to 1

- `max_iter`:

  integer, maximum number of iterations, defaults to 1000

- `max_iter_ls`:

  integer, maximum number of iterations for the line search, defaults to
  1000

- `abs_tol`:

  numeric, tolerance for convergence, defaults to NULL (use rel_tol
  instead)

- `rel_tol`:

  numeric, relative tolerance for convergence, defaults to 1e-5

- `r`:

  numeric, backtracking line search parameter, defaults to 0.8

- `proj`:

  function, projection function to enforce domain of support

- `norm`:

  function, we pass the difference of successive theta updates to the
  norm to check for convergence, defaults to the infinity norm.

- `debug`:

  logical, output debugging information if TRUE; default FALSE

- `trace`:

  logical, if TRUE store the path of the search in the \`path\`
  attribute of the output; default FALSE
