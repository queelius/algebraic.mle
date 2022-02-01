#' Accepts a list of \code{mle} objects for some parameter, say \code{theta},
#' and combines them into a single estimator \code{mle_weighted}
#'
#' It combines the \code{mle} objects by adding them together, weighted by
#' the inverse of their respective variance-covariance matrix. Intuitively,
#' the higher the variance, the less weight an \code{mle} is given in the
#' summation.
#'
#' @param mles A list of \code{mle} objects, all for the same parameter.
#' @return an object of type \code{mle_weighted} which inherits from \code{mle}.
#'
#' @export
mle_weighted <- function(mles)
{
    A <- fisher_info(mles[[1]])
    B <- A %*% point(mles[[1]])
    info.wt <- A
    loglike <- loglike(mles[[1]])
    n <- nobs(mles[[1]])

    for (i in 2:length(mles))
    {
        n <- n + nobs(mles[[i]])
        loglike <- loglike + loglike(mles[[i]])
        A <- fisher_info(mles[[i]])
        info.wt <- info.wt + A
        B <- B + A %*% point(mles[[i]])
    }

    cov.wt <- MASS::ginv(info.wt)
    theta.wt <- cov.wt %*% B

    structure(list(
        theta.hat=matrix(theta.wt,nrow=nrow(info.wt)),
        loglike=loglike,
        info=info.wt,
        sigma=cov.wt),
        class=c("mle_weighted","mle"))
}
