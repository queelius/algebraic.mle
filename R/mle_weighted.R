#' Accepts a list of `mle` objects for some parameter, say `theta`,
#' and combines them into a single estimator `mle_weighted`.
#'
#' It combines the `mle` objects by adding them together, weighted by
#' the inverse of their respective variance-covariance matrix (information
#' matrix). Intuitively, the higher the variance, the less weight an `mle`
#' is given in the summation, or alternatively, the more information it
#' has about the parameter, the more weight it is given in the summation.
#'
#' Each `mle` object should have an `observed_fim` method, which returns
#' the Fisher information matrix (FIM) for the parameter. The FIM is
#' assumed to be the negative of the expected value of the Hessian of the
#' log-likelihood function. The `mle` objects should also have a `params`
#' method, which returns the parameter vector.
#'
#' We assume that the observations used to estimate each of the MLE objects
#' in `mles` are independent.
#'
#' @param mles A list of `mle` objects, all for the same parameter.
#' @return An object of type `mle_weighted` (which inherits from
#'         `mle`) that is the weighted sum of the `mle` objects.
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

    fims <- lapply(mles, observed_fim)
    info.wt <- Reduce(`+`, fims)
    cov.wt <- ginv(info.wt)
    theta.wt <- as.vector(cov.wt %*%
        Reduce(`+`, Map(`%*%`, fims, lapply(mles, params))))

    mle(theta.hat = theta.wt,
        loglike = NULL,
        score = NULL,
        sigma = cov.wt,
        info = info.wt,
        obs = NULL,
        nobs = sum(sapply(mles, nobs)),
        superclasses = c("mle_weighted")
    )
}
