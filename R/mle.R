#' @export
nparams.mle <- function(x) length(x$theta.hat)

#' @export
aic.mle <- function(x) -2 * loglike(x) + 2 * nparams(x)

#' @export
nobs.mle <- function(x) x$sample_size

#' @export
loglike.mle <- function(x,...) x$loglike

#' Function to compute the confidence intervals of \code{mle} objects.
#'
#' @param object the \code{mle} object to compute the confidence intervals for
#' @param parm parameter indexes to compute the confidence intervals for,
#'             defaults to all
#' @param level confidence level, defaults to 0.95 (alpha=.05)
#' @param ... additional arguments to pass
#'
#' @importFrom stats confint
#' @export
confint.mle <- function(object, parm=NULL, level=0.95, ...)
{
    sigma <- diag(vcov(object,...))
    theta <- point(object,...)
    p <- length(theta)
    q <- stats::qnorm(level)
    if (is.null(parm))
        parm <- 1:p

    parm <- parm[parm >= 1 & parm <= p]
    ci <- matrix(nrow=length(parm),ncol=2)
    colnames(ci) <- c(paste((1-level)/2*100,"%"),
                      paste((1-(1-level)/2)*100,"%"))

    i <- 1
    for (j in parm)
    {
        ci[i,] <- c(theta[j] - q * sqrt(sigma[j]),
                    theta[j] + q * sqrt(sigma[j]))
        i <- i + 1
    }
    rownames(ci) <- parm
    ci
}

#' Method for sampling from an \code{mle} object.
#'
#' It creates a sampler for the \code{mle} object. It returns a function
#' that accepts a single parameter \code{n} denoting the number of samples
#' to draw from the \code{mle} object.
#'
#' @param x the \code{mle} object to create sampler for
#' @param ... additional arguments to pass to \code{mle} objects, like sampling
#'            method, \code{?mvtnorm}
#' @export
sampler.mle <- function(x,...)
{
    sigma <- vcov(x)
    theta <- point(x)

    function(n=1)
    {
        mvtnorm::rmvnorm(n,theta,sigma,...)
    }
}

#' @importFrom stats vcov
#' @export
vcov.mle <- function(object,...)
{
    object$sigma
}

#' Computes the MSE of an MLE
#'
#' The MSE of an estimator \code{theta.hat} is defined as
#'     \code{mse(theta.hat) = E[(theta.hat-theta)^2] = vcov(theta.hat) + (bias(theta.hat))^2}.
#' where bias = E(theta.hat-theta), i.e., the MSE depends on theta.
#' The definition of the MSE of \code{theta.hat} is given by
#' the difference between the
#' $$
#'     f(theta.hat)
#' $$
#'
#' @param x the \code{mle} object to compute the MSE of.
#' @param ... this is not used for object type \code{mle}.
#' @export
mse.mle <- function(x,...)
{
    sum(diag(vcov(x)))
}

#' Computes the point estimate of an mle object.
#'
#' @param x the mle object(s).
#' @param ... unused by \code{mle} objects. particular specializations of
#'            \code{mle} objects may use it, however.
#' @export
point.mle <- function(x,...)
{
    x$theta.hat
}

#' Function for obtaining the fisher information matrix of an \code{mle} object.
#'
#' @param x the \code{mle} object to obtain the fisher information of.
#' @param ... unused by \code{mle} objects. particular specializations of
#'            \code{mle} objects may use it, however.
#' @export
fisher_info.mle <- function(x, ...)
{
    x$info
}

#' @export
summary.mle <- function(object,...)
{
    cat("Maximum likelihood estimator, of type",class(object)[1],",\n")
    cat("is normally distributed with mean\n")
    print(point(object))
    cat("and variance-covariance\n")
    print(vcov(object))
    cat("---\n")
    cat("The asymptotic mean squared error",mse(object),"\n")
    cat("The asymptotic 95% confidence interval is\n")
    print(confint(object))
    cat("The log-likelihood is",loglike(object),"\n")
    cat("The AIC is",aic(object),"\n")
}
