#' MLE of the rate parameter when we assume the sample is i.i.d. and drawn
#' from the exponential distribution.
#'
#' Of course, the draws are unlikely to be exponential, but the exponential may
#' be a sufficiently good model. Hypothesis testing, such as relative
#' likelihoods, can be used to assess the appropriateness of the exponential
#' model to the data.
#'
#' @param x a sample of observations
#' @param keep_obs store the observations with the `mle` object, default is `F`
#' @return an `mle` object.
#' @export
mle_exp_rate <- function(x, keep_obs = F) {
    n <- length(x)
    stopifnot(n > 0)
    rate.hat <- 1 / mean(x)
    l <- n * log(rate.hat) - n
    names(rate.hat) <- c("rate")
    mle(
        theta.hat = rate.hat,
        loglike = l,
        score = 0, # n/rate.hat-sum(x)) == 0
        sigma = rate.hat^2 / n,
        info = n / rate.hat^2,
        obs = if (keep_obs) x else NULL,
        nobs = n,
        superclasses = c("mle_exp_rate")
    )
}

#' log-likelihood function generator given data `x` for the exponential
#' distribution
#'
#' @param x data
#' @export
exp_rate_loglike <- function(x) {
    n <- length(x)
    s <- sum(x)
    function(rate) matrix(n * log(rate) - rate * s)
}

#' score (derivative of log-likelihood) function generator given data `x`
#' for the exponential distribution
#'
#' @param x data
#'
#' @export
exp_rate_score <- function(x) {
    n <- length(x)
    s <- sum(x)
    function(rate) matrix(n / rate - s)
}

#' Fisher information function generator for the exponential
#' distribution
#'
#' @param x data
#'
#' @export
exp_rate_fisher_info <- function(x) {
    n <- length(x)
    function(rate) matrix(n / rate^2)
}

#' Computes the bias of an `exp_mle` object (exponential mle).
#'
#' An unbiased estimator of the rate parameter of the exponential distribution
#' is given by: `1/(nobs(x)-1)*bias(x)`, where `x` is an
#' `mle_exp_rate` object.
#'
#' @param x the `mle_exp_rate` object to compute the bias of.
#' @param par the true rate parameter value. Usually, rate is not known,
#'             and so we estimate the bias
#' @param ... pass additional arguments
#' @export
bias.mle_exp_rate <- function(x, par = NULL, ...) {
    if (is.null(par)) {
        par <- point(x)
    }
    stopifnot(length(par) == nparams(x))
    b <- c(1 / (nobs(x) - 1) * par)
    names(b) <- c("bias(rate)")
    b
}
