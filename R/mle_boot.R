#' @title Bootstrap MLE
#'
#' @description
#' Sometimes, the large sample asymptotic theory of MLEs is not applicable.
#' In such cases, we can use the bootstrap to estimate the sampling distribution
#' of the MLE.
#'
#' Look up the \code{boot} package for more information on the bootstrap and
#' how to use it. You can pass additional arguments to the \code{boot} function
#' using the \code{...} argument to \code{mle_boot}.
#'
#' @param mle_solver given a data, find the MLE.
#' @param data data for resampling, where for each resample we generate an MLE
#' @param R bootstrap replicates
#' @param ... additional arguments to pass.
#' @return an \code{mle_boot} object, which is an \code{mle} object with a
#' \code{boot} object as its parent.
#' @importFrom boot boot
#' @export
#' @examples
#' data <- rnorm(15)
#' solver <- function(data, ind) {
#'     point(mle_normal_mu_var(data[ind]))
#' }
#' mle_boot(solver, data, 1000)
mle_boot <- function(mle_solver, data, R, ...) {
    stopifnot(is.function(mle_solver))
    theta.boot <- boot::boot(
        data = data,
        statistic = mle_solver,
        R = R, ...
    )
    class(theta.boot) <- c("mle_boot", "mle", class(theta.boot))
    theta.boot
}

#' Method for obtaining the parameters of an \code{boot} object.
#'
#' @param x the \code{boot} object to obtain the parameters of.
#' @export
params.mle_boot <- function(x) point(x)

#' Method for obtaining the number of parameters of an \code{boot} object.
#'
#' @param x the \code{boot} object to obtain the number of parameters of
#'
#' @export
nparams.mle_boot <- function(x) length(x$t0)

#' Method for obtaining the number of observations in the sample used by
#' an \code{mle}.
#'
#' @param object the \code{mle} object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @importFrom stats nobs
#' @export
nobs.mle_boot <- function(object, ...) length(object$data)

#' Method for obtaining the observations used by the \code{mle}.
#'
#' @param object the \code{mle} object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @export
obs.mle_boot <- function(object, ...) object$data

#' Computes the variance-covariance matrix of \code{boot} object.
#'
#' @param object the \code{boot} object to obtain the variance-covariance of
#' @param ... additional arguments to pass
#'
#' @importFrom stats cov
#' @export
vcov.mle_boot <- function(object, ...) {
    cov(object$t, ...)
}

#' Computes the variance-covariance matrix of \code{boot} object.
#'
#' @param object the \code{boot} object to obtain the variance-covariance of
#' @param ... additional arguments to pass
#'
#' @importFrom MASS ginv
#' @export
fim.mle_boot <- function(object, ...) {
    ginv(vcov(object, ...))
}


#' Computes the estimate of the MSE of a \code{boot} object.
#'
#' @param x the \code{boot} object to compute the MSE of.
#' @param par if the true parameter value is known, you may provide it;
#'            otherwise we use the MLE of \code{par}.
#' @param ... pass additional arguments
#' @export
mse.mle_boot <- function(x, par = NULL, ...) {
    if (is.null(par)) {
        par <- point(x)
    }
    mean(rowSums(t(t(x$t) - as.vector(par))^2))
}

#' Computes the estimate of the bias of a \code{boot} object.
#'
#' Normally, the bias is defined as \eqn{\mathbb{E}(\hat{\theta} - \theta)}.
#' However, in this case, we have a bootstrap sample, so we compute the
#' bias as \eqn{\mathbb{E}(\hat{\theta} - \hat{\theta}_0)}. This is
#' equivalent to the bias of the MLE estimator.
#'
#' Generally, we do not trust this to be a good estimate of the bias
#' of the MLE estimator, but it is still useful for comparing different
#' estimators.
#'
#' @param x the \code{boot} object to compute the bias of.
#' @param par if the true parameter value is known, you may provide it;
#'            otherwise we use the MLE of \code{par}.
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

#' Computes the point estimate of an \code{mle} object.
#'
#' @param x the \code{boot} object.
#' @param ... pass additional arguments
#' @export
point.mle_boot <- function(x, ...) {
    x$t0
}

#' Function for obtaining an estimate of the standard error of the bootstrap of
#' the MLE \code{object}.
#'
#' @param object the bootstrapped MLE object
#' @export
se.mle_boot <- function(object) {
    sqrt(diag(vcov(object)))
}

#' Method for sampling from an \code{mle_boot} object.
#'
#' It creates a sampler for the \code{mle_boot} object. It returns a function
#' that accepts a single parameter \code{n} denoting the number of samples
#' to draw from the \code{mle_boot} object.
#'
#' Unlike the \code{sampler} method for the more general \code{mle} objects,
#' for \code{mle_boot} objects, we sample from the bootstrap replicates, which
#' are more representative of the sampling distribution, particularly for
#' small samples.
#'
#' @param x the \code{mle_boot} object to create sampler for
#' @param ... additional arguments to pass
#' @export
sampler.mle_boot <- function(x, ...) {
    function(n) {
        x$t[sample.int(nrow(x$t), n, replace = TRUE), ]
    }
}
