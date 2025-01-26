## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>"
)

options(digits=3)
options(scipen=999)

## ----setup, warning=F, message=F----------------------------------------------
library(algebraic.mle)
library(algebraic.dist)
library(boot)

## ----sim----------------------------------------------------------------------
n <- 100
err <- 0.1
shape <- 2
scale <- 10
theta = c(shape, scale)
set.seed(142334)

## -----------------------------------------------------------------------------
x <- rweibull(n = n, shape = shape, scale = scale) +
  rnorm(n = n, mean = 0, sd = err)

## -----------------------------------------------------------------------------
head(x, n = 4)

## ----histo, fig.align='center', echo=F----------------------------------------
library(ggplot2)
ggplot(data.frame(x = x), aes(x = x)) +
    geom_histogram(color = "dark blue",
                   fill = "light blue",
                   bins = 25) +
    labs(title = "Simulated Data",
         subtitle = "Histogram")

## -----------------------------------------------------------------------------
fit_normal <- function(data) {
    loglik <- function(theta) {
        sum(dnorm(data, mean = theta[1], sd = sqrt(theta[2]), log = TRUE))
    }
    mu.hat <- mean(data)
    sigma2.hat <- mean((data - mu.hat)^2)
    H <- -numDeriv::hessian(loglik, c(mu.hat, sigma2.hat))
    mle(theta.hat = c(mu.hat, sigma2.hat),
        loglike = loglik(c(mu.hat, sigma2.hat)),
        score = numDeriv::grad(loglik, c(mu.hat, sigma2.hat)),
        sigma = MASS::ginv(H),
        info = H,
        obs = data,
        nobs = length(data),
        superclasses = c("mle_normal"))
}

fit_weibull <- function(data) {
    loglik <- function(theta) {
        sum(dweibull(data, shape = theta[1], scale = theta[2], log = TRUE))
    }

    sol <- stats::optim(
        par = c(shape, scale),
        fn = loglik,
        hessian = TRUE,
        method = "L-BFGS-B",
        lower = c(0, 0),
        #method = "Nelder-Mead",
        control = list(maxit = 10000, fnscale = -1))

    mle(theta.hat = sol$par,
        loglike = sol$value,
        sigma = MASS::ginv(-sol$hessian),
        info = -sol$hessian,
        obs = data,
        nobs = length(data),
        superclasses = c("mle_weibull"))
}

bias.mle_normal <- function(x, theta = NULL) {
    if (is.null(theta))
        theta <- params(x)
    c(0, -theta[2] / nobs(x))
}

theta.hat <- fit_normal(x)
summary(theta.hat)

theta.weibull <- fit_weibull(x)
summary(theta.weibull)

## ----fig.align='center', fig.height=4, fig.width=4, echo=F--------------------
df1 <- data.frame(Value=rweibull(n=100000,
    shape=params(theta.weibull)[1],
    scale=params(theta.weibull)[2]), Group = "Weibull")

df2 <- data.frame(Value=rnorm(n=100000,
    mean=params(theta.hat)[1],
    sd=sqrt(params(theta.hat)[2])), Group = "Normal")

df3 <- data.frame(Value=rweibull(n = 100000,shape = shape, scale = scale) +
                        rnorm(n = 100000, mean = 0, sd = err), Group = "True")
df <- rbind(df1, df2, df3)
ggplot(df, aes(x = Value, fill = Group)) +
  geom_density(alpha = 0.3, adjust = 1.5) +
  scale_fill_manual(values = c(
    Weibull="red",
    Normal="green",
    True="purple")) +
  labs(title = "Density plot", x = "Value", y = "Density") +
  theme_minimal()

## -----------------------------------------------------------------------------
bias(theta.hat,theta)

## -----------------------------------------------------------------------------
bias(theta.hat)

## -----------------------------------------------------------------------------
round(mse(theta.hat), digits=3)
round(mse(theta.weibull), digits=3)  # true MSE

## -----------------------------------------------------------------------------
A <- matrix(c(2,3),nrow=2)
b <- c(1,0)
g <- function(x) A %*% x + b

## -----------------------------------------------------------------------------
g.mc <- rmap(theta.weibull,g,n=100000L)
g.delta <- rmap(theta.weibull,g,method="delta")
round(vcov(g.mc), digits=3)
round(vcov(g.delta), digits=3)

## -----------------------------------------------------------------------------
N <- 500
r <- 5
samp <- rnorm(N, mean = theta[1], sd = sqrt(theta[2]))
samp.sub <- matrix(samp, nrow = r)
mles.sub <- list(length = r)
for (i in 1:r)
    mles.sub[[i]] <- fit_normal(samp.sub[i,])

mle.wt <- mle_weighted(mles.sub)
mle <- fit_normal(samp)

## -----------------------------------------------------------------------------
params(mle.wt)
round(mse(mle.wt), digits=3)

## -----------------------------------------------------------------------------
params(mle)
round(mse(mle), digits=3)

## -----------------------------------------------------------------------------
theta.boot <- mle_boot(
    boot(data = x,
         statistic = function(x, i) params(fit_weibull(x[i])),
         R = 1000))

## -----------------------------------------------------------------------------
print(theta.boot)
bias(theta.boot)

## -----------------------------------------------------------------------------
cramer.test <- function(obs.dat,ref.dat)
{
  stat <- CDFt::CramerVonMisesTwoSamples(obs.dat,ref.dat)
  list(p.value=exp(-stat)/6.0,
       cramer.stat=stat,
       obs.size=length(obs.dat),
       ref.size=length(ref.dat))
}

wei.shape <- params(theta.weibull)[1]
wei.scale <- params(theta.weibull)[2]
ref.dat <- rweibull(1000000, shape = wei.shape, scale = wei.scale)
cramer.test(x, ref.dat)

## -----------------------------------------------------------------------------
norm.mu <- params(theta.hat)[1]
norm.var <- params(theta.hat)[2]
ref.dat <- rnorm(1000000, mean = norm.mu, sd = sqrt(norm.var))
cramer.test(x, ref.dat)

## -----------------------------------------------------------------------------
aic(theta.weibull)
aic(theta.hat)

## -----------------------------------------------------------------------------
pred(x=theta.hat, samp=function(n=1, theta) rnorm(n,theta[1],theta[2]))

## -----------------------------------------------------------------------------
mu <- params(theta.hat)[1]
sd <- sqrt(params(theta.hat)[2])
c(mean=mu,lower=qnorm(.025,mean=mu, sd=sd),upper=qnorm(.975,mean=mu, sd=sd))

