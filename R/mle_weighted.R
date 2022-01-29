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
#' @examples
#' data1 <- rexp(100,2)
#' data2 <- rexp(200,2)
#' mle1 <- mle_exp(data1)
#' cat(point(mle1))
#' mle2 <- mle_exp(data2)
#' cat(point(mle2))
#' mle.wt <- mle_weighted(list(mle1,mle2))
#' summary(mle.wt)
#' mle <- mle_exp(c(data1,data2))
#' point(mle)
#' aic(mle)
#' loglike(mle)
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

    #cov.wt <- solve(info.wt,diag(1,nrow(info.wt)))
    cov.wt <- MASS::ginv(info.wt)
    theta.wt <- cov.wt %*% B

    structure(list(
        theta.hat=matrix(theta.wt,nrow=nrow(info.wt)),
        loglike=loglike,
        info=info.wt,
        sigma=cov.wt),
        class=c("mle_weighted","mle","estimate","normal","dist"))
}
