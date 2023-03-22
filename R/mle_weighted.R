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
#' We assume that the observations used to estimate the MLE objects are
#' independent, which might not always be the case. If this is not the case,
#' then the computation of the log-likelihood is necessarily meaningful.
#'
#' @param mles A list of \code{mle} objects, all for the same parameter.
#'
#'
#' @return An object of type \code{mle_weighted} (which inherits from
#'         \code{mle}, which is the weighted sum of the \code{mle} objects.
#'         The \code{mle_weighted} object has the following attributes:
#'
#' @importFrom MASS ginv
#' @export
mle_weighted <- function(mles) {
    if (is.null(mles)) {
        stop("Invalid input: null list of 'mle' objects.")
    }

    # If only one mle object, return it.
    if (is_mle(mles)) {
        return(mles)
    }

    if (!is.list(mles)) {
        stop("Invalid input: not a list of 'mle' objects.")
    }

    if (!all(sapply(mles, is_mle))) {
        stop("Invalid input: not all elements are 'mle' objects.")
    }

    if (length(mles) == 1) {
        return(mles[[1]])
    }

    fims <- lapply(mles, fim)
    info.wt <- Reduce(`+`, fims)
    cov.wt <- ginv(info.wt)
    theta.wt <- as.vector(cov.wt %*%
        Reduce(`+`, Map(`%*%`, fims, lapply(mles, point))))

    names(theta.wt) <- get_first_attr(mles, point, c(names, colnames))

    mle(
        theta.hat = theta.wt,
        loglike = sum(sapply(mles, loglike)),
        score = mean(sapply(mles, score)),
        sigma = cov.wt,
        info = info.wt,
        obs = lapply(mles, obs),
        nobs = sum(sapply(mles, nobs)),
        superclasses = c("mle_weighted")
    )
}
