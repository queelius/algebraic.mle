
<!-- README.md is generated from README.Rmd. Please edit that file -->

# R package: `algebraic.mle`

<!-- badges: start -->
<!-- badges: end -->

An algebra over maximum likelihood estimators (MLE).

MLEs have many desirable, well-defined statistical properties. We define
an algebra over MLEs.

## Installation

You can install the development version of `algebraic.mle` from
[GitHub](https://github.com/) with:

``` r
#install.packages("devtools")
devtools::install_github("queelius/algebraic.mle")
#> Downloading GitHub repo queelius/algebraic.mle@HEAD
#> 
#>      checking for file ‘/tmp/RtmpoJ0zsT/remotes6d18d50f6e26a/queelius-algebraic.mle-3e1c165/DESCRIPTION’ ...  ✓  checking for file ‘/tmp/RtmpoJ0zsT/remotes6d18d50f6e26a/queelius-algebraic.mle-3e1c165/DESCRIPTION’
#>   ─  preparing ‘algebraic.mle’:
#>    checking DESCRIPTION meta-information ...  ✓  checking DESCRIPTION meta-information
#>   ─  checking for LF line-endings in source and make files and shell scripts
#>   ─  checking for empty or unneeded directories
#>    Omitted ‘LazyData’ from DESCRIPTION
#>   ─  building ‘algebraic.mle_0.1.0.tar.gz’
#>      
#> 
#> Installing package into '/usr/local/lib/R/site-library'
#> (as 'lib' is unspecified)
```

## MLE of exponential distribution’s rate parameter

In what follows, to demonstrate the `algebraic.mle` R package, we
consider a random sample from exponentially distributed random variable
*X*<sub>*i*</sub> ∼ EXP (*λ*=1) for *i* = 1, …, *n*.

We generate an observation from this random sample with:

``` r
library(stats)
library(tidyverse)
n <- 1000
rate <- 1
x <- tibble(x=stats::rexp(n,rate))
```

We have observed a sample of size *n* = 1000. We show some observations
from this sample (data frame) with:

``` r
print(x)
#> # A tibble: 1,000 × 1
#>         x
#>     <dbl>
#>  1 0.146 
#>  2 0.182 
#>  3 1.35  
#>  4 0.0481
#>  5 1.25  
#>  6 0.189 
#>  7 1.43  
#>  8 0.111 
#>  9 0.0494
#> 10 0.617 
#> # … with 990 more rows
```

We show a histogram of the sample, and a plot of the exponential
function’s pdf, with:

``` r
library(ggplot2)
ggplot(x, aes(x=x)) + geom_histogram(aes(y=..density..),alpha=.2) +
    xlim(0,6) + 
    geom_function(fun=dexp)
```

<img src="man/figures/README-unnamed-chunk-5-1.png" width="100%" />

If we would like to estimate *λ*, we can do so using maximum likelihood
estimation as implemented by the `algebraic.mle` package:

``` r
library(algebraic.mle)
#> 
#> Attaching package: 'algebraic.mle'
#> The following object is masked from 'package:stats':
#> 
#>     nobs
rate.hat <- mle_exp(x$x)
summary(rate.hat)
#> Maximum likelihood estimator, of type mle_exp ,
#> is normally distributed with mean
#> [1] 0.9509951
#> and variance-covariance
#>              [,1]
#> [1,] 0.0009043917
#> ---
#> The asymptotic mean squared error 0.0009043917 
#> The asymptotic 95% confidence interval is
#>       2.5 %   97.5 %
#> 1 0.9015292 1.000461
#> The log-likelihood is -1050.246 
#> The AIC is 2102.493
```

We can show the point estimator with:

``` r
point(rate.hat)
#> [1] 0.9509951
```

We can show the Fisher information and variance-covariance matrices
with:

``` r
fisher_info(rate.hat)
#>          [,1]
#> [1,] 1105.716
vcov(rate.hat)
#>              [,1]
#> [1,] 0.0009043917
```

(If `rate.hat` had been a vector, `vcov` would have output a
variance-covariance matrix. We may consider the above outputs 1 × 1
matrices.)

We can show the confidence interval with:

``` r
confint(rate.hat)
#>       2.5 %   97.5 %
#> 1 0.9015292 1.000461
```

## Sampling distribution of the MLE

In general, to estimate the sampling distribution, we generate
*B* = 10000 samples (of size 1000) and their corresponding estimators,
*θ̂*<sup>(1)</sup>, …, *θ̂*<sup>(*B*)</sup>.

Normally, we do not have *B* samples, and if we did, we would gather all
*B* samples into one sample (or used a weighted MLE), which would
contain more (Fisher) information about *θ*.

However, a nice property of MLEs is that, asymptotically, they converge
to a normal distribution with a mean given by the true parameter, in
this case *λ*, and a variance-covariance given by the inverse of the
Fisher information matrix, i.e., *λ̂* ∼ *N*(*λ*,*I*<sup>−1</sup>(*λ*))
where *I* is the Fisher information matrix, in this case
*I*(*λ*) = *n*/*λ*<sup>2</sup>.

We observe the empirical sampling distribution of *θ̂* overlaid with the
theoretical asymptotic distribution with:

``` r
B <- 1000
data0 <- numeric(length=B)
for (i in 1:B)
{
    x <- stats::rexp(n,rate)
    data0[i] <- point(mle_exp(x))
}
ggplot(tibble(rate.hat=data0), aes(x=rate.hat)) +
    geom_histogram(aes(y=..density..),alpha=.3,bins=50) +
    geom_function(fun=function(x) { dnorm(x,mean=rate,sd=rate/sqrt(n)) })
```

<img src="man/figures/README-unnamed-chunk-10-1.png" width="100%" />

We do not know *λ*, but we may estimate it from a sample, and thus we
may approximate the sampling distribution of *λ̂* with
*N*(*λ̂*,*I*<sup>−1</sup>(*λ̂*)).

Since we are only giving one sample, we cannot do as we did before to
provide *B* estimates of *λ*. However, we can sample from the
approximation of the asymptotic distribution of *λ̂* with:

``` r
data1 <- sampler(rate.hat)(B)
```

We visually compare the two MLE samples, `data0` and `data1`, with:

``` r
data <- data.frame(values=c(data0,data1), group=c(rep("empirical",B),rep("approximation",B)))

ggplot(data,aes(x=values,fill=group)) +
    geom_histogram(aes(y=..density..),position="identity", alpha=0.2, bins=50)
```

<img src="man/figures/README-unnamed-chunk-12-1.png" width="100%" />

Due to sampling error, we see that the approximation, the estimate of
the asymptotic sampling distribution *N*(*λ̂*,*λ̂*/*n*<sup>1/2</sup>), is
shifted to the left of the sample from the true distribution, but they
appear to be quite similar otherwise.

## Invariance property of the MLE

An interesting property of an MLE *λ̂* is that the MLE of *f*(*λ*) is
given by *f*(*λ̂*).

The method `rmap` applied to an object `x` for which
`is(x,"mle") == TRUE` and a function `f` compatible with `point(x)` (and
optionally a simulation sample size) computes the MLE of `f(x)`.

### Example

We know that the MLE *λ̂* ∼ *N*(*λ*,*λ*<sup>2</sup>/*n*). We seek a
transformation *g*(*λ̂*) such that its expectation is 2*λ*, i.e.,
*g*(*λ*) = 2*λ*:

``` r
f <- function(lambda) 2*lambda
```

We compute the MLE of 2*λ̂* with:

``` r
g.hat <- rmap(rate.hat,f)
```

Then, var (*f*(*λ̂*)) = 4*λ*<sup>2</sup>/*n*. Letting *λ* = 1, we see
that, asymptotically, *g*(*λ̂*) ∼ *N*(2,4/*n*).

We expected a mean of 2 and a variance of 4/*n* = 0.004, which is pretty
close to what we obtained:

``` r
point(g.hat)
#> [1] 1.90199
vcov(g.hat)
#>             [,1]
#> [1,] 0.003469069
```

## Weighted MLE: a weighted sum of maximum likelihood estimators

Since the variance-covariance of an MLE is inversely proportional to the
Fisher information that the MLE is defined with respect to, we can
combine multiple MLEs of *θ*, each of which may be defined with respect
to a different kind of sample, to arrive at the MLE that incorporates
the Fisher information in all of those samples.

Consider *k* mutually independent MLE estimators of parameter *θ*,
*θ̂*<sub>1</sub>, …, *θ̂*<sub>*k*</sub>, where
*θ̂*<sub>*j*</sub> ∼ *N*(*θ*,*I*<sub>*j*</sub><sup>−1</sup>(*θ*)).

Then, the maximum likelihood estimator of *θ* that incorporates all of
the data in *θ̂*<sub>1</sub>, …, *θ̂*<sub>*k*</sub> is given by the
inverse-variance weighted mean,
*θ̂* = (∑*I*<sub>*j*</sub>(*θ*))<sup>−1</sup>(∑*I*<sub>*j*</sub>(*θ*)*θ*<sub>*j*</sub>).

To evaluate the performance of this estimator, we generate a sample of
*N* = 1000 from EXP (*λ*=1). We compute the MLE for the observed sample,
denoted by *θ*. We then divide the observed sample into *r* = 5
sub-samples, each of size *N*/*r* = 100, and compute the MLE for each
sub-sampled, denoted by *θ*<sup>(1)</sup>, …, *θ*<sup>(*r*)</sup>.
Finally, we come these MLEs into the weighted MLE, denoted by
*θ*<sub>*w*</sub>.

``` r
N <- 1000
r <- 5
data3 <- rexp(N)
data3.mat <- matrix(data3,nrow=r)
mles <- list(length=r)
for (i in 1:r)
    mles[[i]] <- mle_exp(data3.mat[i,])

mle.wt <- mle_weighted(mles)
summary(mle.wt)
#> Maximum likelihood estimator, of type mle_weighted ,
#> is normally distributed with mean
#>           [,1]
#> [1,] 0.9592375
#> and variance-covariance
#>              [,1]
#> [1,] 0.0009211142
#> ---
#> The asymptotic mean squared error 0.0009211142 
#> The asymptotic 95% confidence interval is
#>       2.5 %   97.5 %
#> 1 0.9093164 1.009159
#> The log-likelihood is -1040.018 
#> The AIC is 2082.035

mle.tot <- mle_exp(data3)
summary(mle.tot)
#> Maximum likelihood estimator, of type mle_exp ,
#> is normally distributed with mean
#> [1] 0.9602567
#> and variance-covariance
#>              [,1]
#> [1,] 0.0009220928
#> ---
#> The asymptotic mean squared error 0.0009220928 
#> The asymptotic 95% confidence interval is
#>       2.5 %   97.5 %
#> 1 0.9103091 1.010204
#> The log-likelihood is -1040.555 
#> The AIC is 2083.109
```

We see that *θ̂* and *θ̂*<sub>*w*</sub> model approximately the same
sampling distribution for *θ*.
