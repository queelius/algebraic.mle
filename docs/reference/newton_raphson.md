# newton_raphson

Performs Newton-Raphson to find a solution that finds a local maxima of
the objective function \`fn\`. We assume the maxima is an interior point
of the support and that an initial guess par is near the local maxima.

## Usage

``` r
newton_raphson(par, fn, gr, hess, inverted = FALSE, control = list(), ...)
```

## Arguments

- par:

  numeric, initial guess

- fn:

  function, objective function to maximized

- gr:

  function, gradient function

- hess:

  function, hessian function

- control:

  list, control for the local search, see function description.

## Value

list, designed to be mostly consistent with \`optim\` in stats. list
elements: \`par\`: best solution \`value\`: function value \`fn(par)\`
\`gr\`: gradient at \`par\` \`hess\`: hessian at \`par\`
\`convergence\`: convergence status, 0 if converged, 1 if not
(consistent with \`optim\`) \`count\`: a count of \`fn\` and \`gr\`
evaluations (consistent with \`optim\`) \`trace_info\`: a matrix of
trace information (optional)

## Details

Use a global search method like \`sim_anneal\` to find a good initial
guess that is near a global maximum of \`fn\`.

## Functions

- `newton_raphson()`: control

## Fields

- `eta`:

  numeric, learning rate

- `fnscale`:

  numeric, scaling factor for \`fn\`. If negative, then turns the
  problem into a maximization problem.

- `maxit`:

  integer, maximum number of iterations, defaults to 1000

- `convergence`:

  function, check for convergence, defaults to \`\|\|gr(par)\|\|\` and
  \`\|\|(par1 - par0)\|\|\` being approximately zero (gradient test uses
  absolute tolerance and parameter difference test uses relative
  tolerance, see \`abs_tol\` and \`rel_tol\`).

- `proj`:

  function, projection function to enforce domain of support

- `trace`:

  logical, keep track of trace information

- `abs_tol`:

  numeric, absolute tolerance for convergence, which is used by the
  default \`convergence\` function

- `rel_tol`:

  numeric, absolute tolerance for convergence, which is used by the
  default \`convergence\` function

- `REPORT`:

  integer, frequency of tracer reports, defautls to every 10 iterations
