#' Generic method for obtaining the log-likelihood of fitted distribution
#' object, assuming it has a likelihood function, e.g., a fitted `mle` object.
#'
#' @param x the object to obtain the log-likelihood of
#' @param ... additional arguments to pass
#'
#' @export
loglike <- function(x,...)
{
    UseMethod("loglike",x)
}

#' Generic method for obtaining the AIC of a fitted distribution object fit.
#'
#' @param x the object to obtain the AIC of
#'
#' @export
aic <- function(x)
{
    UseMethod("aic",x)
}

#' Method for obtaining the parameters of a fitted distribution object.
#'
#' @param x the fitted object to obtain the parameters of
#'
#' @export
params <- function(x)
{
    UseMethod("params",x)
}


#' Generic method for obtaining the number of parameters of a fitted
#' distribution object.
#'
#' @param x the fitted object to obtain the number of parameters for
#'
#' @export
nparams <- function(x)
{
    UseMethod("nparams",x)
}

#' Generic method for obtaining the best point estimate from an estimator.
#'
#' @param x the object to obtain the point estimate of
#' @param ... additional arguments to pass
#'
#' @export
point <- function(x, ...)
{
    UseMethod("point",x)
}

#' Generic function for sampling from distribution objects.
#'
#' It creates a sampler for the \code{x} object. It returns a function
#' that accepts a single parameter \code{n} denoting the number of samples
#' to draw from the \code{x} object.
#'
#' @param x the \code{x} object to create a sampler for
#' @param ... additional arguments to pass
#'
#' @export
sampler <- function(x, ...)
{
    UseMethod("sampler",x)
}

#' Generic method for obtaining the fisher information
#' matrix of an \code{mle} object.
#'
#' Fisher information is a way of measuring the amount of
#' information that an observable random variable `X`
#' carries about an unknown parameter \code{theta}
#' upon which the probability of `X` depends.
#'
#' The inverse of the Fisher information matrix
#' is the variance-covariance of the MLE for
#' \code{theta}.
#'
#' @param x the object to obtain the fisher information of
#' @param ... additional arguments to pass
#'
#' @export
fisher_info <- function(x, ...)
{
    UseMethod("fisher_info",x)
}

#' Generic function for obtaining the mean squared error (MSE) of an estimator,
#' \code{mse(x) = E[(x-mu)^2]} where \code{mu} is the (likely unknown) true
#' parameter value.
#'
#' @param x the object to compute the MSE of
#' @param ... additional arguments to pass
#'
#' @export
mse <- function(x, ...)
{
    UseMethod("mse",x)
}

#' Computes the bias of an estimator object.
#'
#' @param x the object to compute the bias of.
#' @param ... pass additional arguments
#' @export
bias <- function(x,...)
{
    UseMethod("bias",x)
}


#' Generic function for computing the distribution of \code{f(x)} where \code{x}
#' models a distribution (random element) object.
#'
#' When we apply a function \code{f} to a distribution \code{x}, the result
#' is a distribution object \code{f(x)}.
#'
#' @param x a list of distribution (random element) objects.
#' @param g a function that accepts arguments sampled from the distribution
#'          objects in \code{x}.
#' @param n number of samples for the distribution to draw from to estimate \code{f(x)}
#' @param ... additional arguments to pass.
#' @export
rmap <- function(x,g,n,...)
{
    UseMethod("rmap",x)
}


#' Method for obtaining the observations used by a fitted model \code{object}.
#'
#' @param object the fitted object to obtain the number of observations used by the fit
#' @param ... additional arguments to pass
#' @export
obs <- function(object,...)
{
    UseMethod("obs",object)
}

#' Method for obtaining an estimate of the standard error of an estimator.
#'
#' @param object the estimator
#' @export
se <- function(object)
{
    UseMethod("se",object)
}



