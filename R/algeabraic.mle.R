#' `algebraic.mle`: A package for algebraically operating on and generating
#' maximum likelihood estimators from existing maximum likelihood estimators.
#'
#' The object representing a fitted model is a type of `mle` object, the maximum
#' likelihood estimator of the model with respect to observed data.
#'
#' It has a relatively rich API for working with these objects to help you
#' understand your MLE estimator.#'
#' 
#' @docType package
#' @name algebraic.mle
NULL
#> NULL


#' Generic method for obtaining the log-likelihood value of a fitted MLE
#' object.
#'
#' @param x the object to obtain the log-likelihood of
#' @param ... additional arguments to pass
#'
#' @export
loglik_val <- function(x, ...) {
    UseMethod("loglik_val", x)
}

#' Generic method for obtaining the AIC of a fitted distribution object fit.
#'
#' @param x the object to obtain the AIC of
#'
#' @export
aic <- function(x) {
    UseMethod("aic", x)
}

#' Generic method for obtaining the observed FIM
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
#' @export
observed_fim <- function(x, ...) {
    UseMethod("observed_fim", x)
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
score_val <- function(x, ...) {
    UseMethod("score_val", x)
}

#' Method for obtaining the standard error of an estimator.
#'
#' @param x the estimator
#' @param ... additional arguments to pass
#' @export
se <- function(x, ...) {
    UseMethod("se", x)
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


#' Generic method for obtaining the parameters of a fitted MLE object.
#' 
#' @param x the object to obtain the parameters of
#' @param ... additional arguments to pass
#' @export
params <- function(x, ...) {
    UseMethod("params", x)
}