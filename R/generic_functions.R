#' Generic method for obtaining the log-likelihood of an MLE fit.
#'
#' @param x the object to obtain the log-likelihood of
#' @param ... additional arguments to pass
#'
#' @export
loglike <- function(x,...)
{
    UseMethod("loglike",x)
}

#' Generic method for obtaining the AIC of an MLE fit.
#'
#' @param x the object to obtain the AIC of
#'
#' @export
aic <- function(x)
{
    UseMethod("aic",x)
}

#' Method for obtaining the parameters of a fitted object.
#'
#' @param x the fitted object to obtain the parameters of
#'
#' @export
params <- function(x)
{
    UseMethod("params",x)
}


#' Generic method for obtaining the number of parameters of an object.
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

#' Generic function for computing the distribution of \code{f(x)} where \code{x}
#' models a \code{dist} (distribution) object.
#'
#' When we apply a function \code{f} to an argument \code{x}, we are applying
#' the function to the argument, or evaluating the function by substituting the
#' formal parameters of the function with the arguments. However, if the
#' argument is a random variable, then \code{f(x)} is also a random variable.
#'
#' @param x a list of \code{dist} (distrubtion) objects.
#' @param f a function that accepts arguments sampled from the \code{dist} items
#'          in \code{x}.
#' @param ... additional arguments to pass.
#' @export
rmap <- function(x,f,...)
{
    UseMethod("rmap",x)
}
