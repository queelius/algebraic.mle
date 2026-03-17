# mle_newton_raphson

MLE method using the Newton-Raphson method

## Usage

``` r
mle_newton_raphson(ll, score, info, theta0, inverted = FALSE, options = list())
```

## Arguments

- ll:

  log-likelihood function

- score:

  score Function, the gradient of log-likelihood

- info:

  FIM function

- theta0:

  initial guess

- inverted:

  logical, if TRUE \`info\` is covariance instead of FIM

- options:

  list, options for the local search, see \`mle_local_search\`

## Value

an object of class \`mle_newton_raphson\`, which is an \`mle\` object
