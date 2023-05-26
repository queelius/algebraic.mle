#' MLE of the rate parameter when we assume the sample is i.i.d. and drawn
#' from the weibull distribution.
#'
#' @param x a sample of observations
#' @param par0 initial estimate of the parameters, defaults to `c(1,1)`
#' @param pgtol the convergence tolerance for the optimization algorithm
#' @param keep_obs Boolean, specifies whether to keep observations
#' @return an `mle` object.
#' @importFrom MASS ginv
#' @importFrom stats optim
#' @export
mle_weibull <- function(x, par0 = c(1,1),
    pgtol = 1e-7, maxit = 10000L, keep_obs = FALSE, ...) {

    n <- length(x)
    stopifnot(n > 0, all(par0 > 0), pgtol > 0, maxit >= 1)

    ll <- weibull_loglike(x)
    scr <- weibull_score(x)
    sol <- optim(par = par0, fn = weibull_loglike(x),
        gr = weibull_score(x), method = "L-BFGS-B",
        lower = 0, hessian = TRUE,
        control = list(fnscale = -1,
                       pgtol = 1e-16,
                       maxit = 1000L,
                       ...))
    mle_sol <- mle_numerical(
        sol = sol,
        superclasses = c("mle_weibull", "mle_numerical"))
    mle_sol$obs <- if (keep_obs) x else NULL
    mle_sol$nobs <- n
    mle_sol
}

#' A log-likelihood function generator given data `x` for the weibull
#' distribution.
#'
#' The returned log-likelihood function takes a single vector `par` of
#' size `2` (at least) where the first component is the shape parameter
#' `k` and the second component is the scale parameter `lambda` such
#' that E[X] = lambda * GAMMA(1 + 1/k).
#'
#' @param x data
#' @return log-likelihood function for the weibull distribution
#'         with shape parameter `k` and scale parameter `lambda`
#'         given data `x`. (Accepts a parameter vector `par` of size 2.)
#' @export
weibull_loglike <- function(x) {
    n <- length(x)
    stopifnot(n > 0, all(x > 0)) # we need at least one observation
                                # and all observations must be positive
    sum_log_x <- sum(log(x))
    function (par) {
        k <- par[1] # shape
        l <- par[2] # scale
        n * log(k / l^k) +
            (k - 1) * sum_log_x - sum((x / l)^k)
    }
}

#' score function generator given data `x` for the weibull
#' distributison given a simple random sample.
#'
#' @param x data
#' @export
weibull_score <- function(x) {
    n <- length(x)
    stopifnot(n > 0, all(x > 0)) # we need at least one observation
                                 # and all observations must be positive
    s <- sum(log(x))
    function(par) {
        k <- par[1] # shape
        l <- par[2] # scale
        x_l <- x / l
        x_l_k <- x_l^k
        c(s + n / k - n * log(l) - sum(x_l_k * log(x_l)),
          k * sum(x_l_k) / l - n * k / l)
    }
}

#' weibull_fim
#' 
#' Observed Fisher information generator given data `x` for the weibull
#' distribution
#'
#' @param x observed data
#' @export
weibull_fim <- function(x) {
    n <- length(x)
    function(par) {
        k <- par[1]
        l <- par[2]
        x_l_k <- (x / l)^k
        log_x_l <- log(x / l)
        p_k_l <- n / l - sum(x_l_k * (1 + k * log_x_l)) / l
        matrix(c(
            # d^2/dk^2
            n / k^2 + sum(x_l_k * log_x_l^2),
            # d^2/dkdl
            p_k_l,
            # d^2/dkdl
            p_k_l,
            # d^2/dl^2
            -n * k / l^2 + k * (k + 1) / l^2 * sum(x_l_k)
        ), nrow = 2)
    }
}


#' Log-Likelihood Function Generator for the Weibull Distribution
#' with a shape and scale parameters that are a function of predictors in `df`,
#'
#'    `resp(df) ~ weibull(shape(df, beta), scale(df, beta))`
#'
#' where `beta` is a vector of parameters.
#' 
#' The function returned by this function is suitable for use with
#' `optim` or `nlm` to find the maximum likelihood estimate of `beta`.
#' 
#' Note that the expected value given the data is
#' `scale(...) * gamma(1 + 1/shape(...))`.
#' 
#' @param df a numeric vector representing a sample of observations.
#' @param resp a function that returns the response variable given a data frame.
#' @param shape a function that returns the shape given a data frame
#'             and a parameter vector `beta`
#' @param scale a function that returns the scale given a data frame and a
#'              parameter vector `beta`
#' @return Returns a function that computes the conditional log-likelihood
#' @export
weibull_conditional_shape_scale_loglike <- function(df, resp, shape, scale) {
    function(beta) {
        sum(dweibull(x = resp(df),
            shape = shape(df, beta),
            scale = scale(df, beta), log = TRUE))
    }
}
