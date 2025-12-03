# Generic method for computing the observed FIM of an \`mle\` object.

Fisher information is a way of measuring the amount of information that
an observable random variable \`X\` carries about an unknown parameter
\`theta\` upon which the probability of \`X\` depends.

## Usage

``` r
observed_fim(x, ...)
```

## Arguments

- x:

  the object to obtain the fisher information of

- ...:

  additional arguments to pass

## Details

The inverse of the Fisher information matrix is the variance-covariance
of the MLE for \`theta\`.

Some MLE objects do not have an observed FIM, e.g., if the MLE's
sampling distribution was bootstrapped.
