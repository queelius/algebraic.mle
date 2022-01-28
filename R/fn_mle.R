#' Computes the distribution of \code{f(x)} as a function of \code{x}, where
#' \code{x} models a \code{mle} object.
#'
#' Since \code{x} is an MLE, asymptotically,
#' \code{f(x) ~ N(mean(f(x)),vcov(f(x)))}.
#'
#' @param x an \code{mle} object.
#' @param f a function that accepts objects like x (e.g., a vector).
#' @param n number of samples to take from \code{x} to estimate distribution
#'        of \code{f(x)}.
#' @param ... additional arguments to pass to the \code{mle} sampler.
#' @export
fn_distr.mle <- function(x,f,n=1000,...)
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
        class=c("mle_func","mle","estimator"))
}
