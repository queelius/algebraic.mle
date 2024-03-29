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
#' @param x the `boot` return value
#' @return an `mle_boot` object (wrapper for `boot` object)
#' @export
mle_boot <- function(x) {
    class(x) <- c("mle_boot", "mle", class(x))
    x
}

#' Determine if an object is an `mle_boot` object.
#' @param x the object to test
#' @export
is_mle_boot <- function(x) inherits(x, "mle_boot")

#' Method for obtaining the parameters of an `boot` object.
#'
#' @param x the `boot` object to obtain the parameters of.
#' @importFrom algebraic.dist params
#' @export
params.mle_boot <- function(x) {
    x$t0
}

#' Method for obtaining the number of parameters of an `boot` object.
#'
#' @param x the `boot` object to obtain the number of parameters of
#' @importFrom algebraic.dist nparams
#' @export
nparams.mle_boot <- function(x) {
    length(x$t0)
}

#' Method for obtaining the number of observations in the sample used by
#' an `mle`.
#'
#' @param object the `mle` object to obtain the number of observations for
#' @param ... additional arguments to pass (not used)
#' @importFrom stats nobs
#' @export
nobs.mle_boot <- function(object, ...) {
    if (is.matrix(object$data) || is.data.frame(object$data)) {
        return(nrow(object$data))
    } else {
        length(object$data)
    }
}

#' Method for obtaining the observations used by the `mle`.
#'
#' @param x the `mle` object to obtain the number of observations for
#' @importFrom algebraic.dist obs
#' @export
obs.mle_boot <- function(x) {
    x$data
}

#' Computes the variance-covariance matrix of `boot` object.
#' Note: This impelements the `vcov` method defined in `stats`.
#'
#' @param object the `boot` object to obtain the variance-covariance of
#' @param ... additional arguments to pass into `stats::cov`
#' @importFrom stats vcov cov
#' @export
vcov.mle_boot <- function(object, ...) {
    cov(object$t, ...)
}

#' Computes the estimate of the MSE of a `boot` object.
#'
#' @param x the `boot` object to compute the MSE of.
#' @param theta true parameter value (not used for `mle_boot`)
#' @param ... pass additional arguments into `vcov`
#' @importFrom algebraic.dist params
#' @importFrom stats vcov
#' @export
mse.mle_boot <- function(x, theta = NULL, ...) {
    vcov(x, ...) + bias(x, theta) %*% t(bias(x, theta))
}

#' Computes the estimate of the bias of a `mle_boot` object.
#'
#' @param x the `mle_boot` object to compute the bias of.
#' @param theta true parameter value (not used for `mle_boot`).
#' @param ... pass additional arguments (not used)
#' @importFrom algebraic.dist params nparams
#' @export
bias.mle_boot <- function(x, theta = NULL, ...) {
    if (!is.null(theta)) {
        warning("We ignore the `theta` argument for `mle_boot` objects.")
    }
    if (length(params(x)) == 1) {
        return(mean(x$t) - params(x))
    } else {
        return(colMeans(x$t) - params(x))
    }
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
#' @param ... additional arguments to pass (not used)
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
#' Note: This impelements the `vcov` method defined in `stats`.
#'
#' @param object the `mle_boot` object to obtain the confidence interval of
#' @param parm the parameter to obtain the confidence interval of (not used)
#' @param level the confidence level
#' @param type the type of confidence interval to compute
#' @param ... additional arguments to pass into `boot.ci`
#' @importFrom boot boot.ci
#' @importFrom stats confint
#' @importFrom algebraic.dist nparams
#' @importFrom utils tail
#' @export
confint.mle_boot <- function(object, parm = NULL, level = 0.95,
    type = c("norm", "basic", "perc", "bca"), ...) {

    type <- match.arg(type)
    type_long <- switch(type,
        norm = "normal",
        basic = "basic",
        stud = "student",
        perc = "percent",
        bca = "bca")

    p <- nparams(object)
    alpha <- 1 - level
    CI <- matrix(NA, nrow = p, ncol = 2)
    colnames(CI) <- c(paste0(alpha / 2 * 100, "%"),
                      paste0((1 - alpha / 2) * 100, "%"))

    for (j in 1:p) {
        ci <- boot.ci(object, conf = level, type = type, index = j, ...)
        CI[j, ] <- tail(ci[[type_long]][1,], 2)
    }

    theta <- params(object)
    if (is.null(names(theta))) {
        rownames(CI) <- paste0("param", 1:p)
    } else {
        rownames(CI) <- names(theta)[1:p]
    }
    CI
}
