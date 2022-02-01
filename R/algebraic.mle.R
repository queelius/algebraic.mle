#' \code{algebraic.mle}: A package for algebraically operating on maximum
#' likelihood estimators.
#'
#' The object representing a fitted model is a type of \code{mle} object, the maximum
#' likelihood estimator of the model with respect to observed data.
#' In what follows, we briefly define the API (generic functions, mostly) with
#' default implementations for objects that inherit from \code{mle}:
#'
#' \code{vcov(mle)} returns the variance-covariance matrix of the model's
#' parameter estimates.
#'
#' \code{point(mle)} returns the point that maximizes the likelihood
#'   of the model.
#'
#'\code{confint(mle)} returns confidence intervals of the true parameter value
#'estimated by the \code{mle} object.
#'
#' \code{sampler(mle)} maps to a function that may be used to sample from the
#' sampling distribution of the \code{mle} object.
#'
#' \code{mse(mle)} computes the asymptotic mean squared error of the \code{mle} object.
#'
#' \code{mle_weighted(mles)} computes the weighted maximum likelihood estimate
#' from independent \code{mle} objects in the list argument \code{mles}.
#'
#' \code{fisher_info(mle)} returns the Fisher information matrix of the model's
#'parameters.
#'
#' \code{rmap(mle,g) provides an approximation of the maximum likelihood estimator
#' that is a function of another \code{mle} object. Note that code{numerical_mle}
#'inherits from \code{mle}.
#'
#' Finally, since normal distributions are closed under linear transformations,
#' then letting \code{A} be a \code{p}-by-\code{q} matrix and \code{x} be a
#' \code{q} dimensional \code{mle} object, then \code{A %*% x} is a \code{p}
#' dimensional multivariate normal with a mean \code{A %*% x} and a
#' variance-covariance \code{tr(A) %*% vcov(x) %*% A}.
#' Furthermore, by the invariance property of the MLE, \code{A %*% x} is the
#' MLE of \code{A %*% true(x)} where \code{true(x)} is the parameter being
#' estimated.
#'
#' We define a set of generic functions for multiplying \code{mle} objects by
#' (non-random) matrices and vectors to permit these sort of operations, with
#' the asymptotic distributions of the operators applied to the \code{mle} objects
#' exactly known if the asymptotic distribution of the \code{mle} objects are
#' exactly known.
#'
#' See \code{rmap} for when the operation is non-linear.
#'
#' @docType package
#' @name algebraic.mle
NULL
#> NULL
