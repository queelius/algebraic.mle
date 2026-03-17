# General iterative MLE method

General iterative MLE method

## Usage

``` r
mle_iterative(
  l,
  theta0,
  dir = NULL,
  eps = 1e-05,
  sup = NULL,
  eta = 1,
  r = 0.5,
  max_iter = 0L,
  debug = F
)
```

## Arguments

- l:

  log-likelihood function of type \\R^p \mapsto R\\

- theta0:

  initial guess of \\\theta\\ with \\p\\ components

- dir:

  promising direction function of type \\R^p \mapsto R^p\\

- sup:

  domain of support function of type \\R^p \mapsto \\T,F\\\\ for the
  log-likelihood function \\l\\

- eta:

  learning rate, defaults to 1

- r:

  backing tracking parameter

- max_iter:

  maximum number of iterations

- debug:

  Boolean, output debugging information if TRUE; defaults to FALSE

- stop_cond:

  stopping condition function of type \\R^p \times R^p \mapsto \\T,F\\\\
