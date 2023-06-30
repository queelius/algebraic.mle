#' Empirical sampling distribution of the MLE
#' 
#' NOTE: If you store the results of a Monte-carlo simulation
#' in this wrapper object, then many of the functions defined for this
#' object be misleading. For example, `confint.mle_emp` will
#' return a probability interval. This is actually a still useful
#' metric to have, since we can determine the probability that
#' an estimate will contain the true parameter value, but this is
#' not the , which is
#' typically not what we want. However, it is good for other things, like
#' computing the bias and mean squared error.
#'
#' @param mles the MLEs
#' @param nobs sample size used to compute each MLE
#' @return an `mle_empirical` object
#' @export
mle_emp <- function(
    mles, nobs = NULL) {

    stopifnot(!is.null(mles), !is.list(mles))

    if (!is.matrix(mles)) {
        mles <- as.matrix(mles, ncol = 1)
    }

    object <- list(mles = mles, nobs = nobs)
    class(object) <- c("mle_emp", "mle")
    object
}

#' Determine if an object is an `mle_emp` object.
#' @param x the object to test
#' @export
is_mle_emp <- function(x) inherits(x, "mle_emp")

#' Method for obtaining the parameters of an `mle_emp` object.
#'
#' @param x the `mle_emp` object to obtain the parameters of.
#' @export
params.mle_emp <- function(x) point(x)

#' Method for obtaining the number of parameters of an `mle_emp` object.
#'
#' @param x the `mle_emp` object to obtain the number of parameters of
#'
#' @export
nparams.mle_emp <- function(x) ncol(x$mles)

#' Method for obtaining the number of observations in the sample used by
#' an `mle`.
#'
#' @param x the `mle` object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @importFrom stats nobs
#' @export
nobs.mle_emp <- function(x, ...) x$nobs

#' Computes the variance-covariance matrix of `mle_emp` object.
#'
#' @param x the `mle_emp` object to obtain the variance-covariance of
#' @param ... additional arguments to pass
#'
#' @importFrom stats cov
#' @export
vcov.mle_emp <- function(x, ...) {
    cov(x$mles, ...)
}

#' Computes the estimate of the bias of a `mle_emp` object.
#'
#' @param x the `mle_emp` object to compute the bias of.
#' @param par The true parameter value.
#' @param ... pass additional arguments
#' @export
bias.mle_emp <- function(x, par, ...) {
    stopifnot(!is.null(par), length(par) == nparams(x))
    params(x) - par
}

#' Computes the point estimate of an `mle` object.
#'
#' @param x the `mle` object to compute the point estimate of
#' @param ... pass additional arguments
#' @export
point.mle_emp <- function(x, ...) {
    stopifnot(!is.null(x), nrow(x$mles) > 0)
    colMeans(x$mles)
}

#' Method for sampling from an `mle_emp` object.
#'
#' @param x the `mle_emp` object to create sampler for
#' @param ... additional arguments to pass
#' @export
sampler.mle_emp <- function(x, ...) {
    stopifnot(!is.null(x), nrow(x$mles) > 0)
    function(n) {
        x$mles[sample.int(nrow(x$mles), n, replace = TRUE), ]
    }
}

#' Method for obtained the confidence interval of an `mle_emp` object.
#' @param x the `mle_emp` object to obtain the confidence interval of
#' @param level the confidence level
#' @param ... additional arguments to pass
#' @importFrom stats quantile
#' @export
confint.mle_emp <- function(x, level = 0.95, ...) {

    stopifnot(is.numeric(level), level >= 0, level <= 1)

    alpha <- (1 - level) / 2
    # for each column, return the desired lower and upper percentile
    # of that column
    t(apply(x$mles, 2, quantile, probs = c(alpha, 1 - alpha)))
}

