#' Generic method for obtaining the best point estimate from an estimator.
#'
#' @param x the object to obtain the point estimate of.
#' @param ... additional arguments to pass.
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
#' @param ... additional arguments to pass.
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
#' @param ... additional arguments to pass.
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
#' @param x the object to compute the MSE of.
#' @param ... additional arguments to pass.
#'
#' @export
mse <- function(x, ...)
{
    UseMethod("mse",x)
}

#' Generic function for computing the distribution of \code{f(x)} where \code{x}
#' models another distribution object.
#'
#' @param x a \code{distr} object.
#' @param f a function of the \code{x} object.
#' @param n number of samples to take to estimate \code{f(x)}  with.
#' @param ... additional arguments to pass.
#'
#' @export
fn_distr <- function(x, f, n=1000, ...)
{
    UseMethod("fn_distr",x)
}

