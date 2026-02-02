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

## Value

An `mle` object of class `rmap_mle` representing the transformed MLE
with variance estimated by the specified method.

## Details

We provide two different methods for estimating the variance-covariance
of \`f(x)\`: method = "delta" -\> delta method method = "mc" -\> monte
carlo method

## Examples

``` r
# MLE for normal distribution
set.seed(123)
x <- rnorm(100, mean = 5, sd = 2)
n <- length(x)
fit <- mle(
  theta.hat = c(mu = mean(x), var = var(x)),
  sigma = diag(c(var(x)/n, 2*var(x)^2/n)),
  nobs = n
)

# Transform: compute MLE of standard deviation (sqrt of variance)
# Using delta method
g <- function(theta) sqrt(theta[2])
sd_mle <- rmap(fit, g, method = "delta")
params(sd_mle)
#>      var 
#> 1.825632 
se(sd_mle)
#> [1] 0.1290917
```
