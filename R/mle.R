#' \code{make_mle} makes an \code{mle} object.
#'
#' @param theta.hat the MLE
#' @param loglike the log-likelihood of \code{theta.hat} given the data
#' @param score the score function evaluated at \code{theta.hat}
#' @param sigma the variance-covariance matrix of \code{theta.hat} given that data
#' @param info the information matrix of \code{theta.hat} given the data
#' @param obs observation (sample) data
#' @param sample_size number of observations in \code{obs}
#' @export
make_mle <- function(theta.hat,loglike=NULL,score=NULL,
                     sigma=NULL,info=NULL,obs=NULL,sample_size=NULL)
{
    structure(list(
        theta.hat=theta.hat,
        loglike=loglike,
        score=score,
        sigma=sigma,
        info=info,
        sample_size=sample_size),
        class=c("mle"))
}

#' Method for obtaining the number of observations in the sample used by
#' an \code{mle} object \code{x}.
#'
#' @param x the \code{mle} object to print
#' @param ... additional arguments to pass
#' @export
print.mle <- function(x,...)
{
    x$obs <- NULL
    print.default(x)
}

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

#' Method for obtaining the observations used by the \code{mle}.
#'
#' @param object the \code{mle} object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @export
obs.mle <- function(object,...) object$obs

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
        ci[i,] <- c(theta[j] - q*sqrt(sigma[j]),
                    theta[j] + q*sqrt(sigma[j]))
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
        mvtnorm::rmvnorm(n,theta,sigma,...)
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
#' However, assuming the regularity conditions, the MLE is asymptotically
#' unbiased, and thus for sufficiently large samples,
#' \code{mse(x.hat) = trace(vcov(theta.hat))}.
#'
#' @param x the \code{mle} object to compute the MSE of.
#' @param ... pass additional arguments
#' @export
mse.mle <- function(x,...)
{
    sum(diag(vcov(x)))
}

#' Computes the asymptotic bias of an \code{mle} object.
#'
#' The bias of an estimator is just \code{E(point(mle)-theta)} where \code{theta}
#' is the true parameter value. Assuming the regularity conditions, the bias
#' is 0.
#'
#' @param x the \code{mle} object to compute the bias of.
#' @param ... pass additional arguments
#' @export
bias.mle <- function(x,...)
{
    if (is.matrix(x))
        return(rep(0,nrow(x)))
    else
        return(rep(0,length(x)))
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

#' Function for obtaining a summary of \code{object}, which is a fitted
#' \code{mle} object.
#'
#' @param object the \code{mle} object
#' @param ... pass additional arguments
#'
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

#' Function for obtaining an estimate of the standard error of the MLE
#' \code{object}.
#'
#' @param object the MLE object
#' @export
se.mle <- function(object)
{
    sqrt(diag(vcov(object)))
}


