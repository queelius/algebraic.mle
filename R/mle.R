#' Method for obtaining the parameters of an \code{mle} object.
#'
#' @param x the \code{mle} object to obtain the parameters of
#'
#' @export
params.mle <- function(x) point(x)

#' Method for obtaining the number of parameters of an \code{mle} object.
#'
#' @param x the \code{mle} object to obtain the number of parameters of
#'
#' @export
nparams.mle <- function(x) length(params(x))

#' Method for obtaining the AIC of an \code{mle} object.
#'
#' @param x the \code{mle} object to obtain the AIC of
#'
#' @export
aic.mle <- function(x) -2 * loglike(x) + 2 * nparams(x)

#' Method for obtaining the number of observations in the sample used by
#' an \code{mle}.
#'
#' @param object the \code{mle} object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @importFrom stats nobs
#' @export
nobs.mle <- function(object,...) object$sample_size

#' Method for obtaining the log-likelihood of an \code{mle} object.
#'
#' @param x the log-likelihood \code{l} evaluated at \code{x}, \code{l(x)}.
#' @param ... additional arguments to pass
#'
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
#' @param ... additional arguments to pass
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

#' Computes the variance-covariance matrix of \code{mle} objects.
#'
#' @param object the \code{mle} object to obtain the variance-covariance of
#' @param ... additional arguments to pass
#'
#' @importFrom stats vcov
#' @export
vcov.mle <- function(object,...)
{
    object$sigma
}

#' Computes the asymptotic MSE of an \code{mle} object.
#'
#' The MSE of an estimator is just the expected sum of squared differences,
#' e.g., if the true parameter value is \code{x} and we have an estimator \code{x.hat},
#' then the MSE is
#' \code{mse(x.hat) = E[(x.hat-x)^2] = trace(vcov(x.hat)) + (bias(x.hat))^2}.
#'
#' Since \code{x} is not typically known, the bias is not a statistic, and thus the
#' MSE for biased estimators is not typically a statistic that can be computed.
#' However, since the MLE is asymptotically unbiased, asymptotically,
#' \code{mse(x.hat) = trace(vcov(theta.hat))}.
#'
#' @param x the \code{mle} object to compute the MSE of.
#' @param ... pass additional arguments
#' @export
mse.mle <- function(x,...)
{
    sum(diag(vcov(x)))
}

#' Computes the point estimate of an \code{mle} object.
#'
#' @param x the \code{mle} object.
#' @param ... pass additional arguments
#' @export
point.mle <- function(x,...)
{
    x$theta.hat
}

#' Function for obtaining the fisher information matrix of an \code{mle} object.
#'
#' @param x the \code{mle} object to obtain the fisher information of.
#' @param ... pass additional arguments
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
