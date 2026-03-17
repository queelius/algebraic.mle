# mle_gradient_ascent

MLE method using gradient ascent

## Usage

``` r
mle_gradient_ascent(ll, theta0, score, options)
```

## Arguments

- ll:

  log-likelihood function

- theta0:

  initial guess of theta with \`p\` components

- score:

  score function of type \`R^p -\> R^p\`, the gradient of \`ll\`
