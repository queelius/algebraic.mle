#' @title Bootstrap MLE
#'
#' @description
#' Sometimes, the large sample asymptotic theory of MLEs is not applicable.
#' In such cases, we can use the bootstrap to estimate the sampling distribution
#' of the MLE.
#'
#' Look up the `boot` package for more information on the bootstrap and
#' how to use it. You can pass additional arguments to the `boot` function
#' using the `...` argument to `mle_boot`.
#'
#' @param mle_solver given a data, find the MLE.
#' @param data data for resampling, where for each resample we generate an MLE
#' @param R bootstrap replicates, defaults to 999
#' @param ... additional arguments to pass.
#' @return an `mle_boot` object, which is an `mle` object with a
#' `boot` object as its parent.
#' @importFrom boot boot
#' @export
#' @examples
#' data <- rnorm(15)
#' solver <- function(data, ind) {
#'     point(mle_normal(data[ind]))
#' }
#' theta.bt <- mle_boot(solver, data, 1000)
#' 
mle_boot <- function(mle_solver, data, R=999, ...) {
    stopifnot(is.function(mle_solver))
    theta.bt <- boot(
        data = data,
        statistic = mle_solver,
        R = R, ...
    )
    class(theta.b) <- c("mle_boot", "mle", class(theta.bt))
    theta.bt
}

#' Method for obtaining the parameters of an `boot` object.
#'
#' @param x the `boot` object to obtain the parameters of.
#' @export
params.mle_boot <- function(x) point(x)

#' Method for obtaining the number of parameters of an `boot` object.
#'
#' @param x the `boot` object to obtain the number of parameters of
#'
#' @export
nparams.mle_boot <- function(x) length(x$t0)

#' Method for obtaining the number of observations in the sample used by
#' an `mle`.
#'
#' @param object the `mle` object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @importFrom stats nobs
#' @export
nobs.mle_boot <- function(object, ...) length(object$data)

#' Method for obtaining the observations used by the `mle`.
#'
#' @param object the `mle` object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @export
obs.mle_boot <- function(object, ...) object$data

#' Computes the variance-covariance matrix of `boot` object.
#'
#' @param object the `boot` object to obtain the variance-covariance of
#' @param ... additional arguments to pass
#'
#' @importFrom stats cov
#' @export
vcov.mle_boot <- function(object, ...) {
    # remove the Bessel correction, since this is an MLE
    cov(object$t, ...) * (nobs(object) - 1) / nobs(object
}

#' Computes the estimate of the MSE of a `boot` object.
#'
#' @param x the `boot` object to compute the MSE of.
#' @param par if the true parameter value is known, you may provide it;
#'            otherwise we use the MLE of `par`.
#' @param ... pass additional arguments
#' @export
mse.mle_boot <- function(x, par = NULL, ...) {
    if (is.null(par)) {
        par <- point(x)
    }
    mean(rowSums(t(t(x$t) - as.vector(par))^2))
}

#' Computes the estimate of the bias of a `mle_boot` object.
#'
#' Generally, we do not trust this to be a good estimate of the bias
#' of the MLE estimator, but it is still useful for comparing different
#' estimators.
#'
#' @param x the `mle_boot` object to compute the bias of.
#' @param par if the true parameter value is known, you may provide it;
#'            otherwise we use the MLE of `par`.
#' @param ... pass additional arguments
#' @export
bias.mle_boot <- function(x, par = NULL, ...) {
    if (is.null(par)) {
        par <- point(x)
    }
    stopifnot(length(par) == nparams(x))

    if (length(par) == 1) {
        return(mean(x$t) - par)
    } else {
        return(colMeans(x$t) - par)
    }
}

#' Computes the point estimate of an `mle` object.
#'
#' @param x the `boot` object.
#' @param ... pass additional arguments
#' @export
point.mle_boot <- function(x, ...) {
    x$t0
}

#' Function for obtaining an estimate of the standard error of the bootstrap of
#' the MLE `object`.
#'
#' @param object the bootstrapped MLE object
#' @export
se.mle_boot <- function(object) {
    sqrt(diag(vcov(object)))
}

#' Method for sampling from an `mle_boot` object.
#'
#' It creates a sampler for the `mle_boot` object. It returns a function
#' that accepts a single parameter `n` denoting the number of samples
#' to draw from the `mle_boot` object.
#'
#' Unlike the `sampler` method for the more general `mle` objects,
#' for `mle_boot` objects, we sample from the bootstrap replicates, which
#' are more representative of the sampling distribution, particularly for
#' small samples.
#'
#' @param x the `mle_boot` object to create sampler for
#' @param ... additional arguments to pass
#' @export
sampler.mle_boot <- function(x, ...) {
    function(n) {
        x$t[sample.int(nrow(x$t), n, replace = TRUE), ]
    }
}
