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

An object of class `mle_numerical` (subclass of `mle`).

## Examples

``` r
# Fit exponential distribution using optim
set.seed(123)
x <- rexp(100, rate = 2)

# Log-likelihood for exponential distribution
loglik <- function(rate) {
  if (rate <= 0) return(-Inf)
  sum(dexp(x, rate = rate, log = TRUE))
}

# Optimize (maximize by setting fnscale = -1)
result <- optim(
  par = 1,
  fn = loglik,
  method = "Brent",
  lower = 0.01, upper = 10,
  hessian = TRUE,
  control = list(fnscale = -1)
)

# Wrap in mle_numerical
fit <- mle_numerical(result, options = list(nobs = length(x)))
params(fit)
#> [1] 1.91256
se(fit)
#> [1] 0.191256
```
