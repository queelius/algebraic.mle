
#' A function for estimating the empirical sampling distribution the
#' MLE given a MLE solver and a sample using the Bootstrap method.
#'
#' @param mle_solver given a data, find the MLE.
#' @param data data for resampling, where for each resample we generate an MLE
#' @param R bootstrap replicates
#' @param ... additional arguments to pass.
#' @importFrom boot boot
#' @export
mle_boot <- function(mle_solver,data,R=599,...)
{
    stopifnot(!is.function(mle_solver))

    boot(data=data,
         statistic=function(x,idx) point(mle_solver(x[idx,])),
         R=R,...)
}

#' A function for computing the sampling distribution of a statistic of the
#' MLE's sampling distribution using the Bootstrap method.
#'
#' @param mle a fitted \code{mle} object.
#' @param loglike.gen a generator for the log-likelihood function; it accepts
#'                observations and constructs the log-likelihood function
#' @param data data for generating MLEs for the bootstrap resampling.
#' @param R bootstrap replicates
#' @param method method for solving the MLE, defaults to numerically solving
#'               the root of the gradient of the log-likelihood using
#'               Newton-raphson.
#' @param ... additional arguments to pass.
#' @export
mle_boot_loglike <- function(mle,loglike.gen,data=NULL,R=NULL,method=mle_newton_raphson,...)
{
    solver <- function(xs) method(
        l=loglike.gen(xs),
        theta0=point(mle),...)

    boot(data=data,
         statistic=function(x,idx) point(solver(x[idx,])),
         R=R)
}

#' Method for obtaining the parameters of an \code{boot} object.
#'
#' @param x the \code{boot} object to obtain the parameters of.
#' @export
params.boot <- function(x) x$t0

#' Method for obtaining the parameter estimate of a \code{boot} object.
#'
#' @param x the \code{boot} object to obtain the parameters of.
#' @export
point.boot <- function(x) x$t0


#' Method for obtaining the number of parameters of an \code{boot} object.
#'
#' @param x the \code{boot} object to obtain the number of parameters of
#'
#' @export
nparams.boot <- function(x) length(x$t0)

#' Method for obtaining the number of observations in the sample used by
#' an \code{mle}.
#'
#' @param object the \code{mle} object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @importFrom stats nobs
#' @export
nobs.boot <- function(object,...) length(object$data)

#' Method for obtaining the observations used by the \code{mle}.
#'
#' @param object the \code{mle} object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @export
obs.boot <- function(object,...) object$data


#' Function to compute the confidence intervals of \code{mle} objects.
#'
#' @param object the \code{boot} object to compute the confidence intervals for
#' @param parm parameter indexes to compute the confidence intervals for,
#'             defaults to all. NOTE: not implemented so ignored for now.
#' @param level confidence level, defaults to 0.95 (alpha=.05)
#' @param type A vector of character strings representing the type of intervals
#'             required. The value should be any subset of the values
#'             \code{c("norm","basic", "stud", "perc", "bca")} or simply "all".
#' @param ... additional arguments to pass
#'
#' @importFrom stats confint
#' @importFrom boot boot.ci
#' @export
confint.boot <- function(object, parm=NULL, level=0.95, type="perc",...)
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

#' Method for sampling from an \code{boot} object.
#'
#' It creates a sampler for the \code{boot} object. It returns a function
#' that accepts a single parameter \code{n} denoting the number of samples
#' to draw from the \code{boot} object.
#'
#' @param x the \code{boot} object to create sampler for
#' @param ... additional arguments to pass
#' @export
sampler.boot <- function(x,...)
{
    function(n=1) sample(x$t,size=n,replace=T,...)
}

#' Computes the variance-covariance matrix of \code{boot} object.
#'
#' @param object the \code{boot} object to obtain the variance-covariance of
#' @param ... additional arguments to pass
#'
#' @importFrom stats cov
#' @export
vcov.boot <- function(object,...)
{
    cov(object$t)
}

#' Computes the estimate of the MSE of a \code{boot} object.
#'
#' @param x the \code{boot} object to compute the MSE of.
#' @param par if the true parameter value is known, you may provide it;
#'            otherwise we use the MLE of \code{par}.
#' @param ... pass additional arguments
#' @export
mse.boot <- function(x,par=NULL,...)
{
    if (is.null(par))
        par <- point(x)
    mean(rowSums(t(t(x$t)-as.vector(par))^2))
}

#' Computes the estimate of the bias of a \code{boot} object.
#'
#' @param x the \code{boot} object to compute the bias of.
#' @param par if the true parameter value is known, you may provide it;
#'            otherwise we use the MLE of \code{par}.
#' @param ... pass additional arguments
#' @export
bias.boot <- function(x,par=NULL,...)
{
    if (is.null(par))
        par <- point(x)

    m <- length(par)
    b <- numeric(m)
    for (j in 1:m)
        b[j] <- mean(x$t[,j]-par[j])
    b
}

#' Computes the point estimate of an \code{mle} object.
#'
#' @param x the \code{boot} object.
#' @param ... pass additional arguments
#' @export
point.boot <- function(x,...)
{
    x$t0
}

#' Function for obtaining a summary of \code{object}, which is a fitted
#' \code{boot} object.
#'
#' @param object the \code{boot} object
#' @param ... pass additional arguments
#' @export
summary.boot <- function(object,...)
{
    print(object)
    cat("The estimated mean is",point(object))
    cat("The estimated variance-covariance is\n")
    print(vcov(object))
    cat("The estimated mean squared error is",mse(object),"\n")
    cat("The estimated bias is",bias(object),"\n")
    print(confint(object))
}

