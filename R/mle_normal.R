#' MLE of the (mu,sigma) parameter vector when we assume the sample is i.i.d. and
#' drawn from the normal distribution.
#'
#' Of course, the draws are unlikely to be normal, but the normal distribution
#' is often a good model. Hypothesis testing, such as relative
#' likelihoods, can be used to assess the appropriateness of the normal
#' model to the data.
#'
#' @param x a sample of observations
#' @return an \code{mle} object.
#' @export
mle_normal <- function(x)
{
    n <- length(x)
    stopifnot(n > 0)

    theta.hat <- c(mean(x),(n-1)/n*var(x))
    l <- normal_loglike(x)
    info.hat <- -hessian(l,theta.hat)

    make_mle(
        theta.hat=theta.hat,
        loglike=l(theta.hat),
        score=grad(l,theta.hat),
        sigma=ginv(info.hat),
        info=info.hat,
        obs=x,
        sample_size=n)
}

#' log-likelihood function generator given data \code{x} for the normal
#' distribution
#'
#' @param x data
#' @export
normal_loglike <- function(x)
{
    # more efficiently computes:
    #    sum(dnorm(x,mean=theta[1],sd=sqrt(theta[2]),log=T))
    n <- length(x)
    S <- sum(x)
    S2 <- sum(x^2)
    function(theta) -n/2*log(2*pi*theta[2]) -
        1/(2*theta[2])*(S2-2*theta[1]*S+n*theta[1]^2)
}

