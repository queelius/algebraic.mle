
#' Accepts a list of \code{mle} objects for some parameter, say \code{theta},
#' and combines them into a single estimator \code{mle_weighted}.
#' 
#' It combines the \code{mle} objects by adding them together, weighted by
#' the inverse of their respective variance-covariance matrix (information
#' matrix). Intuitively, the higher the variance, the less weight an \code{mle}
#' is given in the summation, or alternatively, the more information it
#' has about the parameter, the more weight it is given in the summation.
#' 
#' Each \code{mle} object should have a \code{fim} method, which returns
#' the information matrix (inverse of the variance-covariance matrix).
#' 
#' @param mles A list of \code{mle} objects, all for the same parameter.
#' @return an object of type \code{mle_weighted} which inherits from \code{mle}.
#'
#' @importFrom MASS ginv
#' @export
mle_weighted <- function(mles)
{
    if (is.null(mles))
        stop("Invalid input: null list of 'mle' objects.")

    if (is_mle(mles))
        return(mles)

    if (!is.list(mles))
        stop("Invalid input: not a list of 'mle' objects.")

    if (!all(sapply(mles, is_mle)))
        stop("Invalid input: not all elements are 'mle' objects.")

    if (length(mles) == 1)
        return(mles[[1]])

    A <- lapply(mles, fim)
    info.wt <- Reduce(`+`, A)
    cov.wt <- ginv(info.wt)
    theta.wt <- cov.wt %*% Reduce(`+`, Map(`%*%`, A, lapply(mles, point)))

    mle(theta.hat=theta.wt,
        loglike=sum(sapply(mles, loglike)),
        score=mean(sapply(mles, score)),
        sigma=cov.wt,
        info=info.wt,
        obs=lapply(mles, obs),
        nobs=sum(sapply(mles, nobs)),
        superclasses=c("mle_weighted"))
}
