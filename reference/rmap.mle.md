# Computes the distribution of \`g(x)\` where \`x\` is an \`mle\` object.

By the invariance property of the MLE, if \`x\` is an \`mle\` object,
then under the right conditions, asymptotically, \`g(x)\` is normally
distributed, g(x) ~ normal(g(point(x)),sigma) where \`sigma\` is the
variance-covariance of \`f(x)\`

## Usage

``` r
# S3 method for class 'mle'
rmap(x, g, ..., n = 1000L, method = c("mc", "delta"))
```

## Arguments

- x:

  an \`mle\` object

- g:

  a function

- ...:

  additional arguments to pass to the \`g\` function

- n:

  number of samples to take to estimate distribution of \`g(x)\` if
  \`method == "mc"\`.

- method:

  method to use to estimate distribution of \`g(x)\`, "delta" or "mc".

## Details

We provide two different methods for estimating the variance-covariance
of \`f(x)\`: method = "delta" -\> delta method method = "mc" -\> monte
carlo method
