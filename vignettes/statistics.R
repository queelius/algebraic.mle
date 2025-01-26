## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")
options(digits = 3)
options(scipen = 999)
print_num <- function(x) print(sprintf("%.3f", x))

## ----install, eval = FALSE----------------------------------------------------
#  install.packages("devtools")
#  devtools::install_github("queelius/algebraic.mle")

## ----setup, include = FALSE---------------------------------------------------
library(algebraic.dist)
library(algebraic.mle)
library(numDeriv)
library(MASS)
library(ggplot2)
library(tibble)
library(mvtnorm)
set.seed(1234)

## -----------------------------------------------------------------------------
fit_normal <- function(data) {
    sigma <- function(data) {
        mean((data - mean(data))^2)
    }
    loglik <- function(par, data) {
        n <- length(data)
        -n / 2 * log(2 * pi * par[2]) - 1 / (2 * par[2]) *
            (sum(data^2) - 2 * par[1] * sum(data) + n * par[1]^2)
    }
    par.hat <- c(mu = mean(data), var = sigma(data))
    H <- numDeriv::hessian(func = loglik, x = par.hat, data = data)
    algebraic.mle::mle(
        theta.hat = par.hat,
        loglike = loglik(par.hat, data),
        score = numDeriv::grad(func = loglik, x = par.hat, data = data),
        sigma = MASS::ginv(-H),
        info = -H,
        obs = NULL,
        nobs = length(data),
        superclasses = c("mle_normal"))
}

## ----theta-samp-mc------------------------------------------------------------
theta_samp_mc <- function(n, theta, B = 10000) {
    mu <- theta[1]
    var <- theta[2]
    mles <- matrix(NA, nrow = B, ncol = 2)
    for (i in 1:B) {
        d <- rnorm(n, mean = mu, sd = sqrt(var))
        mles[i, ] <- params(fit_normal(d))
    }
    colnames(mles) <- c("mu", "var")
    mles
}

## -----------------------------------------------------------------------------
theta.mc <- algebraic.dist::empirical_dist(mles)

## -----------------------------------------------------------------------------
(mu.mc <- mean(theta.mc))

## -----------------------------------------------------------------------------
bias.mle_normal <- function(x, par = NULL, ...) {
    if (is.null(par)) {
        par <- params(x)
    }
    c(mu = 0, var = -(1 / nobs(x)) * par[2])
}

## ----bias-1-------------------------------------------------------------------
# first, we sample some data from the true distribution
data <- rnorm(n = n, mean = mu, sd = sqrt(var))

# now we fit it to the normal distribution
theta.hat <- fit_normal(data)

# now we compute the bias, first using the asymptotic theory
bias(theta.hat, theta)
# now using the empirical sampling distribution
expectation(theta.mc, function(x) x - theta) # mean(theta.mc) - theta

## ----bias-plot, echo=FALSE, fig.height=4, fig.width=6-------------------------
plot(ns, bias_var, type = "l",
    xlab = "Sample size", ylab = "Bias(Var)")
# plot the asymptotic bias
lines(ns, -(1 / ns) * var, lty = 2)

## -----------------------------------------------------------------------------
round(vcov(theta.hat), digits=3)
round(vcov(theta.mc), digits=3)

## ----confint-1----------------------------------------------------------------
confint(theta.hat)

## ----mse-1--------------------------------------------------------------------
mse.hat <- mse(theta.hat, theta)
mse.mc <- matrix(expectation(theta.mc,
    function(x) (x - theta) %*% t(x - theta)), nrow = 2)

round(mse.hat, digits = 3)
round(mse.mc, digits = 3)

## -----------------------------------------------------------------------------
# temporarily show more digits in the numbers/outputs for this code block
options(digits = 12)
# mse(mu)
expectation(theta.mc, function(x) (x[1] - mu)^2)
# variance(mu)
(mu.var <- expectation(theta.mc, function(x) (x[1] - mean(theta.mc)[1])^2))
b <- expectation(theta.mc, function(x) x[1] - mu)
# mse = bias^2 + variance
b^2 + mu.var
options(digits = 3)

## ----mse-plots, echo = FALSE, fig.height = 4, fig.width = 6-------------------
plot(ns, mses.mc[, 1], type = "l", col = rgb(0,0,1,0.5), ylab = "MSE(mu)", xlab = "Sample size") # blue with transparency
lines(ns, mses.hat[, 1], col = rgb(0,1,0,0.5)) # green with transparency
lines(ns, mses.hat.hat[, 1], col = rgb(1,0,0,0.5)) # red with transparency
legend("topright", legend = c("MC", "Asymptotic", "Asymptotic (estimate)"), 
    col = c(rgb(0,0,1,0.5), rgb(0,1,0,0.5), rgb(1,0,0,0.5)), lty = 1)

plot(ns, mses.mc[, 2], type = "l", col = rgb(0,0,1,0.5), ylab = "MSE(var)", xlab = "Sample size") # blue with transparency
lines(ns, mses.hat[, 2], col = rgb(0,1,0,0.5)) # green with transparency
lines(ns, mses.hat.hat[, 2], col = rgb(1,0,0,0.5)) # red with transparency
legend("topright", legend = c("MC", "Asymptotic", "Asymptotic (estimate)"), 
    col = c(rgb(0,0,1,0.5), rgb(0,1,0,0.5), rgb(1,0,0,0.5)), lty = 1)

## -----------------------------------------------------------------------------
params(theta.boot)
confint(theta.boot)

## -----------------------------------------------------------------------------
theta.b <- empirical_dist(theta.boot$t)

## ----bias-2-------------------------------------------------------------------
bias(theta.boot)
expectation(theta.mc, function(x) x - theta)
bias(theta.hat, theta)

## -----------------------------------------------------------------------------
round(vcov(theta.b), digits = 3)
round(vcov(theta.mc), digits = 3)
round(vcov(theta.hat), digits = 3)

## -----------------------------------------------------------------------------
samp <- function(n, par) rnorm(n = n, mean = par[1], sd = sqrt(par[2]))
pred(x = theta.hat, samp = samp)

## -----------------------------------------------------------------------------
par <- params(theta.hat)
mu.hat <- par[1]
var.hat <- par[2]
c(mu.hat, qnorm(c(.025,.975), mean = mu.hat, sd = sqrt(var.hat)))

