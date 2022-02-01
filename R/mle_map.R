#' \code{rmap.mle} computes the distribution of \code{f(x)} as a function of
#' \code{mle} object \code{x}.
#'
#' When we apply \code{f} to an argument \code{x}, denoted by \code{f(x)}, we
#' are evaluating \code{f} by substituting the formal parameters of \code{f}
#' with the argument \code{x}. However, if the argument is a random variable,
#' then \code{f(x)} is also a random variable.
#'
#' By the invariance property of the MLE, if \code{x} is an \code{mle} object,
#' asymptotically, \code{f(x)} is normally distributed with an approximation
#' given by \code{f(x) ~ normal(f(point(x)),sigma)} where \code{sigma} is an
#' estimate of the variance-covariance of \code{f(x)}, e.g., the sample
#' covariance of \code{f(x1),...f(xn)} where \code{xj} is sampled from
#' \code{sampler(x)}.
#'
#' @param x a list of \code{mle} objects.
#' @param f a function that accepts objects like x (e.g., a vector).
#' @param n number of samples to take from \code{x} to estimate distribution
#'        of \code{f(x)}.
#' @param ... additional arguments to pass to the \code{mle} sampler.
#' @export
rmap.mle <- function(x,f,n=1000,...)
{
    theta.hat <- f(point(x))
    data <- sampler(x,...)(n)
    p <- length(f(data[1]))
    f.data <- matrix(nrow=n,ncol=p)
    i <- 1
    for (x in data)
    {
        fx <- f(x)
        f.data[i,] <- fx
        i <- i + 1
    }

    sigma <- NULL
    info <- NULL
    if (p == 1)
    {
        sigma <- stats::var(f.data)
        info <- 1/sigma
    }
    else
    {
        sigma <- stats::cov(f.data)
        info <- MASS::ginv(sigma)
    }

    structure(list(
        sigma=sigma,
        info=info,
        theta.hat=theta.hat),
        class=c("mle_func","mle"))
}
