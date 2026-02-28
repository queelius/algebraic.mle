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

#' @importFrom algebraic.dist params
#' @export
algebraic.dist::params

#' @importFrom algebraic.dist nparams
#' @export
algebraic.dist::nparams

#' @importFrom algebraic.dist obs
#' @export
algebraic.dist::obs

#' @importFrom algebraic.dist sampler
#' @export
algebraic.dist::sampler

#' @importFrom algebraic.dist rmap
#' @export
algebraic.dist::rmap

#' @importFrom algebraic.dist expectation
#' @export
algebraic.dist::expectation

#' @importFrom algebraic.dist marginal
#' @export
algebraic.dist::marginal

#' @importFrom algebraic.dist as_dist
#' @export
algebraic.dist::as_dist

#' @importFrom algebraic.dist cdf
#' @export
algebraic.dist::cdf

#' @importFrom algebraic.dist conditional
#' @export
algebraic.dist::conditional

#' @importFrom algebraic.dist inv_cdf
#' @export
algebraic.dist::inv_cdf

#' @importFrom algebraic.dist sup
#' @export
algebraic.dist::sup

#' @importFrom algebraic.dist normal mvn empirical_dist
NULL

#' Generic method for obtaining the log-likelihood value of a fitted MLE
#' object.
#'
#' @param x the object to obtain the log-likelihood of
#' @param ... additional arguments to pass
#' @return The log-likelihood value (numeric).
#' @export
loglik_val <- function(x, ...) {
    UseMethod("loglik_val", x)
}

#' Generic method for obtaining the AIC of a fitted distribution object fit.
#'
#' @param x the object to obtain the AIC of
#' @return The Akaike Information Criterion value (numeric).
#' @export
aic <- function(x) {
    UseMethod("aic", x)
}

#' Generic method for computing the observed FIM
#' of an `mle` object.
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
#' Some MLE objects do not have an observed FIM,
#' e.g., if the MLE's sampling distribution was
#' bootstrapped.
#'
#' @param x the object to obtain the fisher information of
#' @param ... additional arguments to pass
#' @return The observed Fisher Information Matrix.
#' @export
observed_fim <- function(x, ...) {
    UseMethod("observed_fim", x)
}

#' Generic method for computing the mean squared error (MSE) of an estimator,
#' `mse(x) = E[(x-mu)^2]` where `mu` is the true parameter value.
#'
#' @param x the object to compute the MSE of
#' @param theta the true parameter value
#' @return The mean squared error (matrix or scalar).
#' @export
mse <- function(x, theta) {
    UseMethod("mse", x)
}

#' Generic method for computing the bias of an estimator object.
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

#' Generic method for computing the score of an estimator
#' object (gradient of its log-likelihood function evaluated
#' at the MLE).
#'
#' @param x the object to compute the score of.
#' @param ... pass additional arguments
#' @return The score vector evaluated at the MLE.
#' @export
score_val <- function(x, ...) {
    UseMethod("score_val", x)
}

#' Generic method for obtaining the standard errors of an estimator.
#'
#' @param x the estimator
#' @param ... additional arguments to pass
#' @return Vector of standard errors for each parameter.
#' @export
se <- function(x, ...) {
    UseMethod("se", x)
}

#' Generic method for determining the orthogonal parameters of an estimator.
#'
#' @param x the estimator
#' @param tol the tolerance for determining if a number is close enough to zero
#' @param ... additional arguments to pass
#' @return Logical vector or matrix indicating which parameters are orthogonal.
#' @export
orthogonal <- function(x, tol, ...) {
    UseMethod("orthogonal", x)
}

#' Generic method for computing the predictive confidence interval given an estimator object `x`.
#'
#' @param x the estimator object
#' @param alpha (1-alpha)/2 confidence interval
#' @param samp a sampler for random variable that is parameterized by mle `x`
#' @param ... additional arguments to pass
#' @return Matrix of predictive confidence intervals.
#' @export
pred <- function(x, samp = NULL, alpha = .05, ...) {
    UseMethod("pred", x)
}

