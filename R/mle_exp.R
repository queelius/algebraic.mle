#' Maximum Likelihood Estimation (MLE) for the Exponential Distribution
#' 
#' This function computes the Maximum Likelihood Estimation (MLE) of the rate
#' parameter for a given sample that is assumed to be independently and
#' identically distributed (i.i.d) and drawn from the Exponential distribution.
#'
#' @param x a numeric vector representing a sample of observations.
#' @param keep_obs logical. If TRUE, the observations are stored with the `mle`
#'                 object. Default is FALSE.
#' @return Returns an `mle` object that contains:
#' 
#' @examples
#' \dontrun{
#'   obs <- rexp(100, rate = 0.5)
#'   rate.hat <- mle_exp(obs)
#'   summary(rate.hat)
#' }
#' 
#' @export
mle_exp <- function(x, keep_obs = FALSE) {
    n <- length(x)
    stopifnot(n > 0)
    rate.hat <- 1 / mean(x)

    names(rate.hat) <- c("rate")
    mle(theta.hat = rate.hat,
        loglike = n * log(rate.hat) - n,
        score = 0, # n/rate.hat-sum(x)) == 0
        sigma = rate.hat^2 / n,
        info = n / rate.hat^2,
        obs = if (keep_obs) x else NULL,
        nobs = n,
        superclasses = c("mle_exp"))
}

#' Log-Likelihood Function Generator for the Exponential Distribution
#' 
#' This function returns a function of the rate parameter that computes the
#' log-likelihood of the given data assuming an Exponential distribution.
#' 
#' @param x a numeric vector representing a sample of observations.
#' 
#' @return Returns a function that computes the log-likelihood of `x` given a
#'         rate parameter.
#' 
#' @examples
#' \dontrun{
#'   obs <- rexp(100, rate = 0.5)
#'   loglike <- exp_loglike(obs)
#'   plot(loglike, from = 0, to = 2)
#' }
#'
#' @export
exp_loglike <- function(x) {
    n <- length(x)
    s <- sum(x)
    function(rate) n * log(rate) - rate * s
}

#' Score Function Generator for the Exponential Distribution
#' 
#' This function returns a function of the rate parameter that computes the
#' score (derivative of log-likelihood) of the given data assuming an
#' Exponential distribution.
#' 
#' @param x a numeric vector representing a sample of observations.
#' 
#' @return Returns a function that computes the score of `x` given a rate
#'         parameter.
#'
#' @export
exp_score <- function(x) {
    n <- length(x)
    s <- sum(x)
    function(rate) n / rate - s
}

#' Fisher Information Function Generator for the Exponential Distribution
#' 
#' This function returns a function of the rate parameter that computes the
#' Fisher Information of the given data assuming an Exponential distribution.
#' 
#' @param x a numeric vector representing a sample of observations.
#' 
#' @return Returns a function that computes the Fisher Information of `x` given
#'         a rate parameter.
#'
#' @export
exp_fim <- function(x) {
    n <- length(x)
    function(rate) n / rate^2
}

#' Computes the Bias of an Exponential MLE
#'
#' This function computes the bias of the Maximum Likelihood Estimator (MLE) of
#' the rate parameter for an Exponential distribution. The bias is
#' `1/(nobs(x)-1)*par` where `x` is an `mle_exp` object and `par` is the true
#' rate parameter. If `par` is not provided, the point estimate from `x` is
#' used to estimate the bias.
#' 
#' @param x an `mle_exp` object from which to compute the bias.
#' @param par a numeric value representing the true rate parameter. If NULL, the
#'            point estimate from `x` is used.
#' @param ... additional arguments (not used).
#' 
#' @return Returns a named numeric value representing the estimated bias of the rate parameter.
#' 
#' @examples
#' \dontrun{
#'   obs <- rexp(100, rate = 0.5)
#'   rate.hat <- mle_exp(obs)
#'   print(bias(rate.hat, par = 0.5))
#' }
#' 
#' @export
bias.mle_exp <- function(x, par = NULL, ...) {
    if (is.null(par)) {
        par <- point(x)
    }
    stopifnot(length(par) == nparams(x))
    b <- c(1 / (nobs(x) - 1) * par)
    names(b) <- c("bias(rate)")
    b
}
