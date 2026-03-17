# Fitting models to unknown DGPs

## Introduction

We are interested in the generative process that gave rise to the data
we observed. In the real world, systems are quite complex, but we settle
for simpler models for analytical tractability. So, we usually assume:

1.  The sample is i.i.d.

2.  The observations may be reasonably modeled by some relatively simple
    process.

3.  Other factors, like ambient temperature, are either negligible or
    are more or less constant in the sample.

Since this is a simulation, we know the underlying DGP (data generating
process). It’s just \\ T_i = W_i + \epsilon_i \\ where \\ W_i \sim
\operatorname{weibull}(k,\lambda) \\ and \\ \epsilon_i \sim
\operatorname{normal}(0,\sigma). \\

In the real world, we do not know the DGP. In this study, we will assume
that either \\T_1,\ldots,T_n\\ comes from Weibull or Normal. Clearly,
the true DGF is a bit more complicated but still very simple compared to
a more realistic DGP.

Then, the process of parametrically modeling the observed data may take
the following steps:

1.  Visualize the data, e.g., plot a histogram of the data.

2.  Guess which parametric distribution (for the components) might fit
    the observed data for the system lifetime.

3.  Use a statistical test for goodness-of-fit.

4.  Repeat steps 2 and 3 if the measure of goodness of fit is not
    satisfactory.

## Simulation parameters and generation

The simulation parameters are given by:

``` r
library(tibble)
library(stats)

sim.n <- 27
sim.err.sd <- 0.1
sim.shape <- 20
sim.scale <- 3
sim.theta = c(sim.shape,sim.scale)
set.seed(142334)
```

We generate the data with the following R code:

``` r
sim.df <- tibble(lifetime=
  rweibull(n=sim.n, shape=sim.shape, scale=sim.scale) +
  rnorm(n=sim.n, mean=0, sd=sim.err.sd))
```

A few elements from the sample are given by:

    #> # A tibble: 6 × 1
    #>   lifetime
    #>      <dbl>
    #> 1     2.94
    #> 2     2.80
    #> 3     3.10
    #> 4     3.01
    #> 5     3.23
    #> 6     2.95

## Visualizing the data

Visualizing the data is a good first step in the analysis of the data.
If the data is univariate or bivariate, we can plot a histogram of the
data pretty easily (if it’s multivariate, we can plot the marginal
distributions of the data).

We show a histogram of the simulated data below:
![](unknow_dgp_files/figure-html/histo-1.png)

## Parametrically modeling the data

If we only had this sample, what might we conclude? This can be a very
difficult problem. In our case, we know that the simulated data is drawn
from the distribution \\T_i = W_i + \epsilon_i\\ where \\ W_i \sim
\operatorname{weibull}(\lambda = 20, k = 3) \\ and \\ \epsilon_i \sim
\operatorname{normal}(\mu=0,\sigma=0.1). \\ However, in real-world data
sets, we do not know the distribution. So, let us suppose that we do not
know the true distribution of the data.

If we were only interested in, say, *prediction*, and we had a
sufficiently large sample, we could use a non-parametric methods and
“let the data speak for itself.” However, if we are interested in
inference (e.g., explaining the data) or the sample was small, then we
usually need to make some assumptions about the data.

In this case, we will assume that the data is drawn from a parametric
distribution. There are many well-known, named parametric distributions,
e.g., Pareto, Weibull, and Normal, to name a few. From experience, it
seems like the Weibull and the normal might be good fits to the data.
However, note that since the normal distribution permits negative values
to be realized, it may not be an appropriate choice. Still, since these
are only approximations anyway, this may not be a big deal.

## Maximum likelihood estimation

First, let us fit the Weibull distribution by choosing appropriate shape
\\\lambda\\ and scale \\k\\ parameters using the maximum likelihood
estimator.

To find the MLE of \\\theta = (\lambda,k)'\\, we need the log-likelihood
function, which is given by the following R code:

``` r
loglike <- function(theta) sum(dweibull(
  sim.df$lifetime, shape=theta[1], scale=theta[2], log=TRUE))
```

This is the definition of the log-likelihood, where we have used the
built-in R function `dweibull` to compute the sum of the log-density
function for each observation. However, we can also use the
`algebraic.mle` package to compute the log-likelihood function in a more
efficient manner.

``` r
library(algebraic.mle)
#> 
#> Attaching package: 'algebraic.mle'
#> The following object is masked from 'package:stats':
#> 
#>     confint
ll.wei <- weibull_loglike(sim.df$lifetime)
```

An MLE is a point \\(\hat k,\hat\lambda)'\\ that is a maximum of the
log-likelihood function over the support of the parameters. In cases
where a closed solution isn’t possible, we use iterative techniques.

A popular choice is Newton-Raphson, which is a local search method that
relies upon both the gradient and the Hessian of the log-likelihood
(Jocobian of the gradient) to iteratively improve the solution to the
MLE equation. Here is an example of Newton-raphson code:

``` r
# f is the function we want to find the root of
# Jf is the jocabian of f
# x0 is the starting point
newton_raphson <- function(f, Jf, x0, iter=100)
{
  for (i in 1:iter)
  {
    fx <- f(x0)              # new function value
    d <- solve(Jf(x0), fx)   # newton-raphson direction (pointing uphill)
    x0 <- x0 - d             # newton-raphson update
  }
  x0
}
```

Newton-Raphson and other local search methods need a good starting point
(that is already close to the MLE) to converge to the MLE. To find a
good starting point, we use a global search method, Simulated Annealing,
implemented by `sim_anneal` function. Here is the code for finding a
good starting point:

``` r
# find a good starting position
start <- sim_anneal(
    fn = ll.wei,
    par = sim.theta,
    control = list(
        t_init=100,
        t_end=.1,
        fnscale=-1,
        it_per_temp=20,
        maxit=10000L,
        sup=function(theta) all(theta > 0),
        trace=TRUE))
cat("initial guess (k0,lambda0)' = (", start$par, ")'.\n")
#> initial guess (k0,lambda0)' = ( 18.10665 2.956376 )'.
```

To get an idea about what the algorithm is doing, we can plot the trace
of the algorithm. Let’s take a look at some plots.
![](unknow_dgp_files/figure-html/unnamed-chunk-8-1.png)![](unknow_dgp_files/figure-html/unnamed-chunk-8-2.png)

In the convergence plot, we see the history of log-likelihood values as
the algorithm progresses. The left parameter trace plot shows the best
value for the shape parameter, \\k\\, as the algorithm progresses, and
the right plot does the same for the \\\lambda\\ parameter.

With this starting point in hand, we find an MLE with:

``` r
library(algebraic.mle)
mle.wei <- mle_numerical(optim(
    par=start$par,
    fn=ll.wei,
    gr=weibull_score(sim.df$lifetime),
    method="BFGS",
    control=list(fnscale=-1)))
```

The `mle_numerical` function takes anything that is like an `optim`
object and returns an `mle` object for it.

Here’s a summary of it:

``` r
summary(mle.wei)
#> Maximum likelihood estimator of type mle_numerical is normally distributed.
#> The estimates of the parameters are given by:
#> [1] 18.884786  2.953124
#> The MSE of the individual componetns in a multivariate estimator is:
#> <0 x 0 matrix>
#> The MSE of the estimator is   .
#> The log-likelihood is  9.141514 .
#> The AIC is  -14.28303 .
```

A similar result can be obtained with the `mle_weibull` function from
the `algebraic.mle` package:

``` r
summary(mle_weibull(sim.df$lifetime))
#> Maximum likelihood estimator of type mle_weibull is normally distributed.
#> The estimates of the parameters are given by:
#> [1] 18.884786  2.953124
#> The standard error is  2.749064 0.03182087 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>             2.5%     97.5%
#> param1 13.234004 24.535568
#> param2  2.887715  3.018533
#> The MSE of the individual componetns in a multivariate estimator is:
#> [1] 7.557352664 0.001012568
#> The MSE of the estimator is  7.557353 0.02842253 0.02842253 0.001012568 .
#> The log-likelihood is  9.141514 .
#> The AIC is  -14.28303 .
```

In either case, the `mle` object exposes an API that provides a number
of conventient methods, such as estimating the variance-covariance
matrix, confidence intervals, bias, and so on. Here’s the code to print
a summary of the MLE object:

Let’s do the same for the normal distribution. We will use the
`mle_normal` function from the `algebraic.mle` package:

``` r
mle.norm <- mle_normal(sim.df$lifetime)
summary(mle.norm)
#> Maximum likelihood estimator of type mle_normal is normally distributed.
#> The estimates of the parameters are given by:
#>         mu        var 
#> 2.87408422 0.02896918 
#> The standard error is  0.03275565 0.007884412 .
#> The asymptotic 95% confidence interval of the parameters are given by:
#>           2.5%      97.5%
#> mu  2.80675402 2.94141443
#> var 0.01276254 0.04517582
#> The MSE of the individual componetns in a multivariate estimator is:
#> [1] 1.072933e-03 6.331514e-05
#> The MSE of the estimator is  0.001072933 0 0 6.331514e-05 .
#> The log-likelihood is  9.499217 .
#> The AIC is  -14.99843 .
```

Let’s plot the pdfs of the Weibull and normal distributions:
![](unknow_dgp_files/figure-html/unnamed-chunk-13-1.png)

In purple, we have the true density (DGP). In red, we have the Weibull
density. In green, we have the normal density. From the plot, it’s hard
to tell which distribution is a better fit to the DGP.

Interestingly, the tails of the true distribution seem a bit heavier
than the tails of the Weibull and Normal. This may suggest that a
heavier-tailed model may be a better fit, such as the lognormal
distribution, but we will not pursue this.

How do we choose between the Weibull and Normal distributions? We will
discuss this in the next section.

## Goodness of fit

We are fitting a model to the data that does not precisely capture the
generative model \\W\\. So, how good of a fit is it?

We will conduct a goodness of fit test, \\\begin{align} H_0 &: \text{the
data is compatible with the Weibull distribution}\\ H_A &: \text{the
data is not compatible with the Weibull distribution}. \end{align}\\

To perform this test, we will use the Cramer-von Mises test. This test
is based on the Cramer-von Mises statistic, which is a measure of the
distance between the empirical distribution function of the data and the
distribution function of the model. The Cramer-von Mises statistic is
given by \\ \hat D_n^2 = \frac{1}{n}\sum\_{i=1}^n \left(\hat F_n(x_i) -
F(x_i)\right)^2 \\ where \\\hat F_n\\ is the empirical distribution
function of the data and \\F\\ is the distribution function of the
model.

``` r
cramer.test <- function(obs.dat,ref.dat)
{
  stat <- CDFt::CramerVonMisesTwoSamples(obs.dat,ref.dat)
  list(p.value=exp(-stat)/6.0,
       cramer.stat=stat,
       obs.size=length(obs.dat),
       ref.size=length(ref.dat))
}

wei.shape <- point(mle.wei)[1]
wei.scale <- point(mle.wei)[2]
ref.dat <- rweibull(1000000,shape=wei.shape,scale=wei.scale)
cramer.test(sim.df$lifetime,ref.dat)
#> $p.value
#> [1] 0.1566145
#> 
#> $cramer.stat
#> [1] 0.0622083
#> 
#> $obs.size
#> [1] 27
#> 
#> $ref.size
#> [1] 1000000
```

Looking at the \\p\\-value, we see that the data is compatible with the
Weibull distribution. Now, let’s do the same for the normal
distribution:

``` r
norm.mu <- point(mle.norm)[1]
norm.var <- point(mle.norm)[2]
ref.dat <- rnorm(1000000,mean=norm.mu,sd=sqrt(norm.var))
cramer.test(sim.df$lifetime,ref.dat)
#> $p.value
#> [1] 0.1534112
#> 
#> $cramer.stat
#> [1] 0.08287383
#> 
#> $obs.size
#> [1] 27
#> 
#> $ref.size
#> [1] 1000000
```

They are both compatible with the data. However, the Weibull
distribution has a larger \\p\\-value, which may suggest it is a better
fit. We also have the AIC measure of goodness of fit. The AIC is given
by \\ \text{AIC} = -2\log L + 2k, \\ where \\L\\ is the likelihood of
the model and \\k\\ is the number of parameters in the model. The AIC is
a measure of the tradeoff between the goodness of fit and the complexity
of the model.

``` r
aic(mle.wei)
#> [1] -14.28303
aic(mle.norm)
#>       var 
#> -14.99843
```

A lower AIC value indicates a better fit. Thus, according to this
measure, the Weibull distribution is the better fit.

## Conclusion

In this post, we have shown how to fit a Weibull and Normal
distributions to a simulated dataset whose true distribution, while
known, does not have a common name.

We have shown how to compare the two models using the Cramer-von Mises
test and the AIC measure of goodness of fit.

We came to no definitive conclusion about which model is better, but the
Weibull distribution has a larger \\p\\-value from the Cramer-von Mises
test, and a lower AIC value, which serves as some evidence that it is a
better fit. We saw the true DGP is visually different from both the
Weibull and the normal distributions. Notably, the DGP has longer tails
than both, suggesting that an even better fit may be a long-tail
distribution like the log-normal or the Pareto distribution.
