#' rmap.mle
#' 
#' Computes the distribution of `g(x)` where `x` is an `mle` object.
#'
#' @details
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

#' @inheritParams rmap
#' @param n number of samples to take to estimate distribution of `g(x)` if
#'         `method=="mc"`.
#' @param method method to use to estimate distribution of `g(x)`,
#'               "delta" or "mc".
#' @param ... additional arguments to pass to the `mle` sampler if
#'            `method=="mc"` or to the `grad` function if `method=="delta"`.
#' 
#' @method rmap mle
#' @importFrom stats cov
#' @importFrom numDeriv jacobian
#' @importFrom MASS ginv
#' @export
rmap.mle <- function(x,g,n=1000L,method="mc",...)
{
    stopifnot(is.integer(n), n > 0, is_mle(x), is.function(g))

    # let g(theta) = A %*% theta + b for instance to do a linear transformation
    # of the parameters, where A is a matrix and b is a vector

    g.sigma <- NULL
    if (method=="mc") {
        mle.samp <- as.matrix(sampler(x,...)(n),nrow=n)
        p <- length(g(mle.samp[1,]))

        g.mle.samp <- matrix(nrow=n,ncol=p)
        for (i in 1:n)
            g.mle.samp[i,] <- g(mle.samp[i,])
        g.sigma <- cov(g.mle.samp)
    }
    else if (method=="delta") {
        J <- jacobian(g,point(x))
        g.sigma <- J %*% vcov(x) %*% t(J)
    }
    else {
        stop("method must be 'mc' or 'delta'")
    }

    mle(theta.hat=g(point(x)),
        loglike=NULL,
        score=NULL,
        sigma=g.sigma,
        info=ginv(g.sigma),
        obs=NULL,
        nobs=NULL,
        superclasses=c("rmap_mle"))
}
