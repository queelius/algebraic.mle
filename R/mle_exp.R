#' MLE of the rate parameter when we assume the sample is i.i.d. and drawn
#' from the exponential distribution.
#'
#' Of course, the draws are unlikely to be exponential, but the exponential may
#' be a sufficiently good model. Hypothesis testing, such as relative
#' likelihoods, can be used to assess the appropriateness of the exponential
#' model to the data.
#'
#' @param x a sample of observations
#' @return an \code{mle} object.
#' @export
mle_exp <- function(x)
{
    n <- length(x)
    stopifnot(n > 0)
    rate.hat <- 1/mean(x)
    make_mle(
        theta.hat=rate.hat,
        loglike=n*log(rate.hat)-n,
        score=matrix(0),
        sigma=matrix(rate.hat^2/n),
        info=matrix(n/rate.hat^2),
        obs=x,
        sample_size=n)
}

#' log-likelihood function generator given data \code{x} for the exponential
#' distribution
#'
#' @param x data
#' @export
exp_loglike <- function(x)
{
    n <- length(x)
    s <- sum(x)
    function(rate) n*log(rate) - rate*s
}

