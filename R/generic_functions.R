#' Generic method for obtaining the log-likelihood of fitted distribution
#' object, assuming it has a likelihood function, e.g., a fitted `mle` object.
#'
#' @param x the object to obtain the log-likelihood of
#' @param ... additional arguments to pass
#'
#' @export
loglike <- function(x, ...) {
    UseMethod("loglike", x)
}

#' Generic method for obtaining the AIC of a fitted distribution object fit.
#'
#' @param x the object to obtain the AIC of
#'
#' @export
aic <- function(x) {
    UseMethod("aic", x)
}

#' Method for obtaining the parameters of a fitted distribution object.
#'
#' @param x the fitted object to obtain the parameters of
#'
#' @export
params <- function(x) {
    UseMethod("params", x)
}


#' Generic method for obtaining the number of parameters of a fitted
#' distribution object.
#'
#' @param x the fitted object to obtain the number of parameters for
#'
#' @export
nparams <- function(x) {
    UseMethod("nparams", x)
}

#' Generic method for obtaining the best point estimate from an estimator.
#'
#' @param x the object to obtain the point estimate of
#' @param ... additional arguments to pass
#'
#' @export
point <- function(x, ...) {
    UseMethod("point", x)
}

#' Generic function for sampling from distribution objects.
#'
#' It creates a sampler for the `x` object. It returns a function
#' that accepts a single parameter `n` denoting the number of samples
#' to draw from the `x` object.
#'
#' @param x the `x` object to create a sampler for
#' @param ... additional arguments to pass
#'
#' @export
sampler <- function(x, ...) {
    UseMethod("sampler", x)
}

#' Generic method for obtaining the fisher information
#' matrix of an `mle` object.
#'
#' Fisher information is a way of measuring the amount of
#' information that an observable random variable `X`
#' carries about an unknown parameter `theta`
#' upon which the probability of `X` depends.
#'
#' The inverse of the Fisher information matrix
#' is the variance-covariance of the MLE for
#' `theta`.
#'
#' @param x the object to obtain the fisher information of
#' @param ... additional arguments to pass
#'
#' @export
fim <- function(x, ...) {
    UseMethod("fim", x)
}

#' Generic function for obtaining the mean squared error (MSE) of an estimator,
#' `mse(x) = E[(x-mu)^2]` where `mu` is the true parameter value.
#'
#' @param x the object to compute the MSE of
#' @param theta the true parameter value
#' @export
mse <- function(x, theta) {
    UseMethod("mse", x)
}

#' Computes the bias of an estimator object.
#'
#' @param x the object to compute the bias of.
#' @param theta true parameter value. usually, this is unknown (NULL), in which
#'              case we estimate the bias
#' @param ... pass additional arguments
#' @return The bias of the estimator. The return type depends on the specific
#'         method.
#' @export
bias <- function(x, theta, ...) {
    UseMethod("bias", x)
}

#' score
#' 
#' Generic function for computing the score of an estimator
#' object with respect to its log-likelihood.
#'
#' @param x the object to compute the score of.
#' @param ... pass additional arguments
#'
#' @export
score <- function(x, ...) {
    UseMethod("score", x)
}

#' Generic function for computing the distribution of `f(x)` where `x`
#' models a distribution (random element) object.
#'
#' When we apply a function `f` to a distribution `x`, the result
#' is a distribution object `f(x)`.
#'
#' @param x a list of distribution (random element) objects.
#' @param g a function that accepts arguments sampled from the distribution
#'          objects in `x`.
#' @param ... additional arguments to pass.
#' @export
rmap <- function(x, g, ...) {
    UseMethod("rmap", x)
}

#' Generic function for obtaining the observations used by a fitted model `object`.
#'
#' @param object the fitted object to obtain the number of observations used by the fit
#' @param ... additional arguments to pass
#' @export
obs <- function(object, ...) {
    UseMethod("obs", object)
}

#' Method for obtaining the standard error of an estimator.
#'
#' @param object the estimator
#' @export
se <- function(object) {
    UseMethod("se", object)
}

#' Method for determining the orthogonal parameters of an estimator.
#'
#' @param x the estimator
#' @param tol the tolerance for determining if a number is close enough to zero
#' @param ... additional arguments to pass
#' @export
orthogonal <- function(x, tol, ...) {
    UseMethod("orthogonal", x)
}

#' Compute the predictive confidence interval given an estimator object `x`.
#'
#' @param x the estimator object
#' @param alpha (1-alpha)/2 confidence interval
#' @param samp a sampler for random variable that is parameterized by mle `x`
#' @param ... additional arguments to pass
#' @export
pred <- function(x, samp = NULL, alpha = .05, ...) {
    UseMethod("pred", x)
}

