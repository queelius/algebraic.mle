#' Computes the distribution of `g(x)` where `x` is an `mle` object.
#'
#' By the invariance property of the MLE, if `x` is an `mle` object,
#' then under the right conditions, asymptotically, `g(x)` is normally
#' distributed,
#'     g(x) ~ normal(g(point(x)),sigma)
#' where `sigma` is the variance-covariance of `f(x)`
#'
#' We provide two different methods for estimating the
#' variance-covariance of `f(x)`:
#'     method = "delta" -> delta method
#'     method = "mc" -> monte carlo method
#'
#' @param x an `mle` object
#' @param g a function 
#' @param ... additional arguments to pass to the `g` function
#' @param n number of samples to take to estimate distribution of `g(x)` if
#'         `method == "mc"`.
#' @param method method to use to estimate distribution of `g(x)`,
#'               "delta" or "mc".
#' @importFrom stats cov vcov nobs
#' @importFrom numDeriv jacobian
#' @importFrom algebraic.dist rmap sampler params
#' @importFrom MASS ginv
#' @export
rmap.mle <- function(x, g, ...,
    n = 1000L, method = c("mc", "delta"))
{
    stopifnot(is.integer(n), n > 0, is_mle(x), is.function(g))

    method <- match.arg(method)
    if (method == "mc") {
        mle.samp <- as.matrix(sampler(x)(n), nrow = n)
        p <- length(g(mle.samp[1, ], ...))

        g.mle.samp <- matrix(nrow=n,ncol=p)
        for (i in 1:n)
            g.mle.samp[i,] <- g(mle.samp[i,], ...)
        g.sigma <- cov(g.mle.samp)
    }
    else { # method == "delta"
        J <- jacobian(func = g, x = params(x), ...)
        g.sigma <- J %*% vcov(x) %*% t(J)
    }

    mle(theta.hat=g(params(x)),
        loglike=NULL,
        score=NULL,
        sigma=g.sigma,
        info=ginv(g.sigma),
        obs=NULL,
        nobs=nobs(x),
        superclasses=c("rmap_mle"))
}
