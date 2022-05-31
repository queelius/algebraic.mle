#' \code{rmap.mle} computes the distribution of \code{f(x)} where \code{x}
#' is an \code{mle} object.
#'
#' By the invariance property of the MLE, if \code{x} is an \code{mle} object,
#' then under the right conditions, asymptotically, \code{f(x)} is normally
#' distributed with an approximation given by
#' \code{f(x) ~ normal(f(point(x)),sigma)} where \code{sigma} is
#' the variance-covariance of \code{f(x)}, e.g., the sample
#' covariance of \code{f(x1),...f(xn)} where \code{xj} is sampled from
#' \code{sampler(x)}.
#'
#' @param x an \code{mle} object.
#' @param g a function that accepts objects like \code{point(x)}.
#' @param n number of samples to take to estimate distribution
#'        of \code{f(x)}.
#' @param ... additional arguments to pass to the \code{mle} sampler.
#' @importFrom stats var
#' @importFrom stats cov
#' @importFrom MASS ginv
#' @export
rmap.mle <- function(x,g,n=1000,...)
{
    # by the invariance property of the MLE, the MLE of f(theta) is f(theta.hat)
    g.theta.hat <- g(point(x))
    mle.samp <- as.matrix(sampler(x,...)(n),nrow=n)
    p <- length(g(mle.samp[1,]))
    g.theta.hat <- as.matrix(g.theta.hat,ncol=p)
    g.mle.samp <- matrix(nrow=n,ncol=p)
    for (i in 1:n)
        g.mle.samp[i,] <- g(mle.samp[i,])
    sigma <- stats::cov(g.mle.samp)
    info <- MASS::ginv(sigma)
    os <- obs(x)

    make_mle(g.theta.hat,
             ifelse(is.null(os),NULL,loglike(x)-sum(log(grad(g,os)))),
             sigma,
             info,
             NULL,
             nobs(x))
}


