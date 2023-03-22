#' A function for estimating the empirical sampling distribution the
#' MLE given a MLE solver and a sample using the Bootstrap method.
#'
#' @param mle_solver given a data, find the MLE.
#' @param data data for resampling, where for each resample we generate an MLE
#' @param R bootstrap replicates
#' @param ... additional arguments to pass.
#' @importFrom boot boot
#' @export
mle_boot <- function(mle_solver,data,R,...)
{
    stopifnot(!is.function(mle_solver))
    boot::boot(data=data,
         statistic=function(x,idx) point(mle_solver(x[idx,])),
         R=R,...)
}

#' Method for obtaining the parameters of an \code{boot} object.
#'
#' @param x the \code{boot} object to obtain the parameters of.
#' @export
params.boot <- function(x) x$t0

#' Method for obtaining the number of parameters of an \code{boot} object.
#'
#' @param x the \code{boot} object to obtain the number of parameters of
#'
#' @export
nparams.boot <- function(x) length(x$t0)

#' Method for obtaining the number of observations in the sample used by
#' an \code{mle}.
#'
#' @param object the \code{mle} object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @importFrom stats nobs
#' @export
nobs.boot <- function(object,...) length(object$data)

#' Method for obtaining the observations used by the \code{mle}.
#'
#' @param object the \code{mle} object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @export
obs.boot <- function(object,...) object$data

#' Computes the variance-covariance matrix of \code{boot} object.
#'
#' @param object the \code{boot} object to obtain the variance-covariance of
#' @param ... additional arguments to pass
#'
#' @importFrom stats cov
#' @export
vcov.boot <- function(object,...)
{
    cov(object$t)
}

#' Computes the estimate of the MSE of a \code{boot} object.
#'
#' @param x the \code{boot} object to compute the MSE of.
#' @param par if the true parameter value is known, you may provide it;
#'            otherwise we use the MLE of \code{par}.
#' @param ... pass additional arguments
#' @export
mse.boot <- function(x,par=NULL,...)
{
    if (is.null(par))
        par <- point(x)
    mean(rowSums(t(t(x$t)-as.vector(par))^2))
}

#' Computes the estimate of the bias of a \code{boot} object.
#'
#' @param x the \code{boot} object to compute the bias of.
#' @param par if the true parameter value is known, you may provide it;
#'            otherwise we use the MLE of \code{par}.
#' @param ... pass additional arguments
#' @export
bias.boot <- function(x,par=NULL,...)
{
    if (is.null(par))
        par <- point(x)

    mean(x$t)-par
}

#' Computes the point estimate of an \code{mle} object.
#'
#' @param x the \code{boot} object.
#' @param ... pass additional arguments
#' @export
point.boot <- function(x,...)
{
    x$t0
}

#' Function for obtaining an estimate of the standard error of the bootstrap of
#' the MLE \code{object}.
#'
#' @param object the bootstrapped MLE object
#' @export
se.boot <- function(object)
{
    #summary(object)$bootSE
    sqrt(diag(vcov(object)))
}
