## ----include = FALSE----------------------------------------------------------
knitr::opts_chunk$set(
  collapse = TRUE,
  comment = "#>")

## ----dgp-1--------------------------------------------------------------------
library(algebraic.dist)
b0 <- -.1
b1 <- 0.5

dgp <- function(n, x) {
    rate <- exp(b0 + b1 * x)
    X <- rexp(n, rate)
    W <- rnorm(n, 0, 1e-3)
    return(X + W)
}

## ----data-1-------------------------------------------------------------------
n <- 75 # number of observations
set.seed(1231) # for reproducibility
df <- data.frame(x = rep(NA, n), y = rep(NA, n))
for (i in 1:n) {
    x <- runif(1, -10, 10)
    y <- dgp(n = 1, x = x)
    df[i, ] <- c(x, y)
}

## ----model-1------------------------------------------------------------------
resp <- function(df) df$y
rate <- function(df, beta) exp(beta[1] + beta[2] * df$x)
loglik <- function(df, resp, rate) {
  function(beta) sum(dexp(x = resp(df), rate = rate(df, beta), log = TRUE))
}

## ----mle-1--------------------------------------------------------------------
library(algebraic.mle)

# initial guess for the parameters
par0 <- c(0, 0)
names(par0) <- c("b0", "b1")

sol <- mle_numerical(optim(par = par0,
    fn = loglik(df, resp, rate),
    control = list(fnscale = -1),
    hessian = TRUE))
summary(sol)

## ----plot-model-1, fig.width = 6, fig.height = 6, fig.align = "center"--------
# plot the x-y points from the data frame
plot(df$x,df$y)

# now overlay a plot of the conditional mean
x <- seq(-10, 10, .1)
b0.hat <- params(sol)[1]
b1.hat <- params(sol)[2]
y.hat <- 1/exp(b0.hat + b1.hat*x)
y <- 1/exp(b0 + b1*x)
lines(x, y, col = "green", lwd = 10)
lines(x, y.hat, col = "blue", lwd = 10)

## ----null-1, warning = FALSE, message = FALSE---------------------------------
# construct null model where b1 = 0
rate_b0_zero <- function(df, b1) exp(b1 * df$x)

# initial guess for the parameters
# fit the model under the null hypothesis
sol2 <- mle_numerical(optim(par = 0,
    fn = loglik(df, resp, rate_b0_zero),
    control = list(fnscale = -1),
    hessian = TRUE,
    method = "BFGS"))
summary(sol2)

## -----------------------------------------------------------------------------
(lrt.sol2 <- -2 * (loglik_val(sol2) - loglik_val(sol)))
pchisq(lrt.sol2, df = 1, lower.tail = FALSE) # compute the p-value

## -----------------------------------------------------------------------------
aic(sol)
aic(sol2)

## ----null-b1-1, warning = FALSE, message = FALSE------------------------------
rate_b1_zero <- function(df, b0) exp(b0)
# fit the model under the null hypothesis
sol3 <- mle_numerical(optim(par = 0,
    fn = loglik(df, resp, rate_b1_zero),
    control = list(fnscale = -1),
    hessian = TRUE,
    method = "BFGS"))
(lrt.sol3 <- -2 * (loglik_val(sol2) - loglik_val(sol)))
pchisq(lrt.sol3, df = 1, lower.tail = FALSE) # compute the p-value

## -----------------------------------------------------------------------------
CI <- confint(sol)
# print the confidence interval
print(CI)

