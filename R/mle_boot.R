#' Bootstrapped MLE
#'
#' Sometimes, the large sample asymptotic theory of MLEs is not applicable.
#' In such cases, we can use the bootstrap to estimate the sampling distribution
#' of the MLE.
#' 
#' This takes an approach similiar to the `mle_numerical` object, which is
#' a wrapper for a `stats::optim` return value, or something that is compatible
#' with the `optim` return value. Here, we take a `boot` object, which is
#' the sampling distribution of an MLE, and wrap it in an `mle_boot` object
#' and then provide a number of methods for the `mle_boot` object that satisfies
#' the concept of an `mle` object.
#'
#' Look up the `boot` package for more information on the bootstrap.
#'
#' @param object the `boot` return value
#' @return an `mle_boot` object (wrapper for `boot` object)
#' @export
mle_boot <- function(object) {
    class(object) <- c("mle_boot", "mle", class(object))
    object
}

#' Determine if an object is an `mle_boot` object.
#' @param x the object to test
#' @export
is_mle_boot <- function(x) inherits(x, "mle_boot")

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
    cov(object$t, ...) * (nobs(object) - 1) / nobs(object)
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
    if (length(x$t) == 1) {
        function(n) {
            x$t[sample.int(nrow(x$t), n, replace = TRUE)]
        }
    } else {
        function(n) {
            x$t[sample.int(nrow(x$t), n, replace = TRUE), ]
        }
    }
}

#' Method for obtained the confidence interval of an `mle_boot` object.
#' @param object the `mle_boot` object to obtain the confidence interval of
#' @param parm the parameter to obtain the confidence interval of
#' @param level the confidence level
#' @param ... additional arguments to pass
#' @importFrom boot boot.ci
#' @export
confint.mle_boot <- function(object, parm, level = 0.95, ...) {
    boot.ci(object, conf = level, ...)
}
