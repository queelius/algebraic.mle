#' MLE of the rate parameter when we assume the sample is i.i.d. and drawn
#' from the exponential distribution.
#'
#' Of course, the draws are unlikely to be exponential, but the exponential may
#' be a sufficiently good model. Hypothesis testing, such as relative
#' likelihoods, can be used to assess the appropriateness of the exponential
#' model to the data.
#'
#' @param x a sample of observations
#' @param keep_obs store the observations with the \code{mle} object, default is \code{F}
#' @return an \code{mle} object.
#' @export
mle_exp <- function(x,keep_obs=F)
{
    n <- length(x)
    stopifnot(n > 0)
    rate.hat <- mle_exp_rate(x)
    mle <- make_mle(
        theta.hat=rate.hat,
        loglike=n*log(rate.hat)-n,
        score=matrix(0), # n/rate.hat-sum(x)) == 0
        sigma=matrix(rate.hat^2/n),
        info=matrix(n/rate.hat^2),
        obs=ifelse(keep_obs,x,NULL),
        sample_size=n,
        superclasses=c("mle_exp"))
}

mle_exp_rate <- function(x)
{
    matrix(1/mean(x))
}

#' log-likelihood function generator given data \code{x} for the exponential
#' distribution
#'
#' @param x data
#' @export
exp_rate_loglike <- function(x)
{
    n <- length(x)
    s <- sum(x)
    function(rate) matrix(n*log(rate) - rate*s)
}

#' score (derivative of log-likelihood) function generator given data \code{x}
#' for the exponential distribution
#'
#' @param x data
#'
#' @export
exp_rate_score <- function(x)
{
    n <- length(x)
    s <- sum(x)
    function(rate) matrix(n/rate - s)
}

#' Fisher information function generator for the exponential
#' distribution
#'
#' @param x data
#'
#' @export
exp_rate_fisher_info <- function(x)
{
    n <- length(x)
    function(rate) matrix(n/rate^2)
}

#' Computes the bias of an \code{exp_mle} object (exponential mle).
#'
#' An unbiased estimator of the rate parameter of the exponential distribution
#' is given by: \code{1/(nobs(x)-1)*bias(x)}, where \code{x} is an
#' \code{mle_exp} object.
#'
#' @param x the \code{mle_exp} object to compute the bias of.
#' @param rate the true rate parameter value
#'
#' @export
bias.mle_exp <- function(x,rate)
{
    matrix(1/(nobs(x)-1)*rate)
}
