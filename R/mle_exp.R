
#' Method for sampling from an \code{mle_exp} object.
#'
#' It creates a sampler for the \code{mle_exp} object. It returns a function
#' that accepts a single parameter \code{n} denoting the number of samples
#' to draw from the \code{mle_exp} object.
#'
#' @param x the \code{mle_exp} object to create sampler for
#' @param ... additional arguments to pass to \code{mle_exp} objects.
#' @export
mle_exp_sampler <- function(x,...)
{
    function(n=1)
    {
        stats::rnorm(n,point(x),sqrt(vcov(x)),...)
    }
}
.S3method("sampler", "mle_exp", mle_exp_sampler)

#' MLE of the rate parameter of the exponential distribution given a
#' random sample of observations drawn from it.
#'
#' Of course, the draws are unlikely to be exponential, but the exponential may
#' be a sufficiently good model. Hypothesis testing, such as relative
#' likelihoods, can be used to assess the appropriateness of the exponential
#' model to the data.
#'
#' @param x vector of draws from exponential distribution.
#' @return an \code{mle_exp} object.
#' @export
mle_exp <- function(x)
{
    rate.hat <- 1/mean(x)
    n <- length(x)

    structure(list(
        theta.hat=c(rate.hat),
        info=matrix(n/rate.hat^2),
        sigma=matrix(rate.hat^2/n),
        sample_size=n),
        class=c("mle_exp","mle","estimator"))
}
