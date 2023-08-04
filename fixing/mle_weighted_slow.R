#' Accepts a list of \code{mle} objects for some parameter, say \code{theta},
#' and combines them into a single estimator \code{mle_weighted}.
#'
#' It combines the \code{mle} objects by adding them together, weighted by
#' the inverse of their respective variance-covariance matrix. Intuitively,
#' the higher the variance, the less weight an \code{mle} is given in the
#' summation.
#'
#' @param mles A list of \code{mle} objects, all for the same parameter.
#' @return an object of type \code{mle_weighted} which inherits from \code{mle}.
#'
#' @importFrom MASS ginv
##' @export
mle_weighted_slow <- function(mles)
{
    stopifnot(is.list(mles))
    if (!all(sapply(mles, is_mle)))
        stop("Invalid input: not all elements are 'mle' objects.")

    A <- fim(mles[[1]])
    B <- A %*% point(mles[[1]])
    info.wt <- A
    loglik <- loglike(mles[[1]])
    n <- nobs(mles[[1]])
    total_obs <- obs(mles[[1]])

    for (i in 2:length(mles))
    {
        n <- n + nobs(mles[[i]])
        total_obs <- append(total_obs,obs(mles[[i]]))
        loglik <- loglik + loglike(mles[[i]])
        A <- fim(mles[[i]])
        info.wt <- info.wt + A
        B <- B + A %*% point(mles[[i]])
    }
    info.wt <- (info.wt + t(info.wt))/2
    cov.wt <- ginv(info.wt)
    cov.wt <- (cov.wt + t(cov.wt))/2
    theta.wt <- cov.wt %*% B

    mle(theta.hat=theta.wt,
        loglike=loglik,
        score=NULL,
        sigma=cov.wt,
        info=info.wt,
        obs=total_obs,
        nobs=n,
        superclasses=c("mle_weighted"))
}
