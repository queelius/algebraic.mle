# This function implements the simulated annealing algorithm, which is a global optimization algorithm that is useful for finding a good starting point for a local optimization algorithm. We do not return this as an MLE object because, to be a good estimate of the MLE, the gradient of \`f\` evaluated at its solution should be close to zero, assuming the MLE is interior to the domain of \`f\`. However, since this algorithm is not guided by gradient information, it is not sensitive to the gradient of \`f\` and instead only seeks to maximize \`f\`.

This function implements the simulated annealing algorithm, which is a
global optimization algorithm that is useful for finding a good starting
point for a local optimization algorithm. We do not return this as an
MLE object because, to be a good estimate of the MLE, the gradient of
\`f\` evaluated at its solution should be close to zero, assuming the
MLE is interior to the domain of \`f\`. However, since this algorithm is
not guided by gradient information, it is not sensitive to the gradient
of \`f\` and instead only seeks to maximize \`f\`.

## Usage

``` r
sim_anneal(par, fn, control = list(), ...)
```

## Arguments

- par:

  Initial guess

- fn:

  Objective function to maximize

- control:

  List of optional arguments

- ...:

  Additional arguments that may be passed; loads into options, and is
  also passed into \`neigh\`

## Value

list, members include: \`par\`: best solution \`value\`: function value
\`fn(par)\` \`fn_count\`: a count of \`fn\` invocations \`accepted\`: a
count of accepted moves \`trace_info\`: a matrix of trace information
(optional)

## Functions

- `sim_anneal()`: control

## Fields

- `t_init`:

  Initial temperature, defaults to 100

- `t_end`:

  Final temperature, defaults to 0

- `alpha`:

  Numberic, cooling factor, defaults to .95

- `REPORT`:

  The frequency of reports if control\$debug \> 0. Defaults to every 100
  iterations.

- `it_per_temp`:

  Ineger, iterations per temperature, defaults to 100

- `maxit`:

  Integer, maximum number of iterations, defaults to 100000

- `accept_p`:

  Acceptance probability function, defaults to \`runif(1) \< exp((val0 -
  val1) / temp)\`, where \`val0\` is the function value at the current
  position and \`val1\` is the function value at the proposed position.
  We pass \`temp\` (temperature), \`val1\` (new value at candidate
  position \`par1\`), \`val0\` (old value at current position), \`it\`
  (iteration), and \`...\` as arguments to \`accept_p\`.

- `fnscale`:

  Scaling factor for \`fn\`, defaults to 1. If negative, then turns the
  problem into a maximization problem.

- `sup`:

  Support function, returns TRUE if \`par\` is in the domain of \`fn\`

- `proj`:

  Projection function, returns a vector in the domain of \`fn\`

- `neigh`:

  Neighborhood function, returns a random neighbor of \`par\`, defaults
  to \`par + rnorm(length(par))\`. We pass \`par\` \`temp\`, \`it\`,
  \`value\`, and \`...\` as arguments to \`neigh\`.

- `trace`:

  logical, whether to store current changes in position and associated
  other values in a \`trace_info\` matrix.
