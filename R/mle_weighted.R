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
#' @return An object of type \code{mle_weighted} (which inherits from
#'         \code{mle}) that is the weighted sum of the \code{mle} objects.
#' @examples
#' # Combine three independent estimates of mean
#' set.seed(123)
#'
#' # Three independent samples
#' x1 <- rnorm(50, mean = 10, sd = 2)
#' x2 <- rnorm(30, mean = 10, sd = 2)
#' x3 <- rnorm(70, mean = 10, sd = 2)
#'
#' # Create MLE objects for each sample
#' make_mean_mle <- function(x) {
#'   n <- length(x)
#'   s2 <- var(x)
#'   mle(theta.hat = mean(x),
#'       sigma = matrix(s2/n),
#'       info = matrix(n/s2),
#'       nobs = n)
#' }
#'
#' fit1 <- make_mean_mle(x1)
#' fit2 <- make_mean_mle(x2)
#' fit3 <- make_mean_mle(x3)
#'
#' # Combine using inverse-variance weighting
#' combined <- mle_weighted(list(fit1, fit2, fit3))
#' params(combined)
#' se(combined)
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
