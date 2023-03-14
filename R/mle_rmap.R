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
#' delta method
#' bootstrap method
#' monte carlo method
#'
#' @param x an \code{mle} object.
#' @param g a function that accepts objects like \code{point(x)}.
#' @param n number of samples to take to estimate distribution
#'        of \code{f(x)}.
#' @param method method to use, defaults to "mc" for monte carlo, but "delta"
#'               for delta method and "bootstrap" for bootstrap are other options
#' @param ... additional arguments to pass to the \code{mle} sampler.
#' @importFrom stats var
#' @importFrom stats cov
#' @importFrom MASS ginv
#' @export
rmap.mle <- function(x,g,n=1000L,method="mc",...)
{
    stopifnot(is.integer(n) && n >= 1)
    stopifnot(is_mle(x))
    stopifnot(is.function(g))

    mle.samp <- as.matrix(sampler(x,...)(n),nrow=n)
    p <- length(g(mle.samp[1,]))

    # by the invariance property of the MLE,
    # the MLE of g(theta) is g(theta.hat)
    g.theta.hat <- as.matrix(g(point(x)),ncol=p)
    g.mle.samp <- matrix(nrow=n,ncol=p)
    for (i in 1:n)
        g.mle.samp[i,] <- g(mle.samp[i,])
    sigma <- stats::cov(g.mle.samp)
    info <- MASS::ginv(sigma)

    #data <- g(obs(x))
    #loglik <- ifelse(is.null(data),NULL,loglike(x)-sum(log(grad(g,data))))

    mle(theta.hat=g.theta.hat,
        loglike=NULL,
        score=NULL,
        sigma=sigma,
        info=info,
        obs=g(obs(x)),
        nobs=nobs(x))
}



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
#' delta method
#' bootstrap method
#' monte carlo method
#'
#' @param x an \code{mle} object.
#' @param g a function that accepts objects like \code{point(x)}.
#' @export
mle_apply <- function(x,g,method="delta",keep_obs=F)
{
    stopifnot(is_mle(x))
    stopifnot(is.function(g))

    if (method=="delta")
    {
        dgdv <- deriv(g)
        sigma <- dgdv %*% vcov(x) %*% t(dgdv)
        mle(theta.hat=g(point(x)),
            loglike=NULL,
            score=NULL,
            sigma=sigma,
            info=ginv(sigma),
            obs=if (keep_obs) g(obs) else NULL,
            nobs=nobs(x))
    }
    else if (method=="mc")
    {
        NULL
    }
    else if (method=="bootstrap")
    {
        NULL
    }
}
