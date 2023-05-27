#' `mle` makes an `mle` object.
#'
#' @param theta.hat the MLE
#' @param loglike the log-likelihood of `theta.hat` given the data
#' @param score the score function evaluated at `theta.hat`
#' @param sigma the variance-covariance matrix of `theta.hat` given that data
#' @param info the information matrix of `theta.hat` given the data
#' @param obs observation (sample) data
#' @param nobs number of observations in `obs`
#' @param superclasses class (or classes) with `mle` as base
#' @export
mle <- function(theta.hat,
                loglike = NULL,
                score = NULL,
                sigma = NULL,
                info = NULL,
                obs = NULL,
                nobs = NULL,
                superclasses = NULL) {
    structure(
        list(
            theta.hat = theta.hat,
            loglike = loglike,
            score = score,
            sigma = sigma,
            info = info,
            obs = obs,
            nobs = nobs
        ),
        class = unique(c(superclasses, "mle"))
    )
}

#' Method for obtaining the number of observations in the sample used by
#' an `mle` object `x`.
#'
#' @param x the `mle` object to print
#' @param ... additional arguments to pass
#' @export
print.mle <- function(x, ...) {
    print(summary(x, ...))
}

#' Method for obtaining the parameters of an `mle` object.
#'
#' @param x the `mle` object to obtain the parameters of
#'
#' @export
params.mle <- function(x, ...) {
    x$theta
}

#' Method for obtaining the number of parameters of an `mle` object.
#'
#' @param x the `mle` object to obtain the number of parameters of
#'
#' @export
nparams.mle <- function(x) length(params(x))

#' Method for obtaining the AIC of an `mle` object.
#'
#' @param x the `mle` object to obtain the AIC of
#'
#' @export
aic.mle <- function(x) -2 * loglike(x) + 2 * nparams(x)

#' Method for obtaining the number of observations in the sample used by
#' an `mle`.
#'
#' @param object the `mle` object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @importFrom stats nobs
#' @export
nobs.mle <- function(object, ...) object$nobs

#' Method for obtaining the observations used by the `mle`.
#'
#' @param object the `mle` object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @export
obs.mle <- function(object, ...) object$obs

#' Method for obtaining the log-likelihood of an `mle` object.
#'
#' @param x the log-likelihood `l` evaluated at `x`, `l(x)`.
#' @param ... additional arguments to pass
#'
#' @export
loglike.mle <- function(x, ...) x$loglike

#' Function to compute the confidence intervals of `mle` objects.
#'
#' @param object the `mle` object to compute the confidence intervals for
#' @param parm parameter indexes to compute the confidence intervals for,
#'             defaults to all
#' @param level confidence level, defaults to 0.95 (alpha=.05)
#' @param ... additional arguments to pass
#'
#' @importFrom stats confint
#' @export
confint.mle <- function(object, parm = NULL, level = .95, ...) {
    stopifnot(is.numeric(level), level >= 0, level <= 1)
    V <- vcov(object, ...)
    if (is.null(V)) stop("No variance-covariance matrix available.")
    if (is.matrix(V)) {
        V <- diag(V)
    }
    sigma <- V

    theta <- point(object, ...)
    p <- length(theta)
    q <- stats::qnorm(level)
    if (is.null(parm)) {
        parm <- 1:p
    }

    parm <- parm[parm >= 1 & parm <= p]
    ci <- matrix(nrow = length(parm), ncol = 2)
    colnames(ci) <- c(
        paste0((1 - level) / 2 * 100, "%"),
        paste0((1 - (1 - level) / 2) * 100, "%")
    )

    i <- 1
    for (j in parm)
    {
        ci[i, ] <- c(
            theta[j] - q * sqrt(sigma[j]),
            theta[j] + q * sqrt(sigma[j])
        )
        i <- i + 1
    }

    if (is.null(names(theta))) {
        rownames(ci) <- paste0("param", parm)
    } else {
        rownames(ci) <- names(theta)[parm]
    }
    ci
}

#' Method for sampling from an `mle` object.
#'
#' It creates a sampler for the `mle` object. It returns a function
#' that accepts a single parameter `n` denoting the number of samples
#' to draw from the `mle` object.
#'
#' @param x the `mle` object to create sampler for
#' @param ... additional arguments to pass
#' @importFrom stats rnorm
#' @importFrom mvtnorm rmvnorm
#' @export
sampler.mle <- function(x, ...) {
    V <- vcov(x)
    mu <- point(x)

    if (nparams(x) == 1L) {
        sd <- sqrt(V)
        function(n = 1) stats::rnorm(n, mean = mu, sd = sd, ...)
    } else {
        function(n = 1) mvtnorm::rmvnorm(n, mu, V, ...)
    }
}


#' Computes the variance-covariance matrix of `mle` objects.
#'
#' @param object the `mle` object to obtain the variance-covariance of
#' @param ... additional arguments to pass
#'
#' @importFrom stats vcov
#' @export
vcov.mle <- function(object, ...) {
    object$sigma
}


#' Computes the MSE of an `mle` object.
#'
#' The MSE of an estimator is just the expected sum of squared differences,
#' e.g., if the true parameter value is `x` and we have an estimator `x.hat`,
#' then the MSE is
#' `mse(x.hat) = E[(x.hat-x)^2] = trace(vcov(x.hat)) + (bias(x.hat))^2`.
#'
#' Since `x` is not typically known, we normally must estimate the bias.
#' For sufficiently large samples, for the MLE assuming the regularity conditions,
#' the bias is given by
#' `mse(x.hat) = trace(vcov(x.hat)) + bias(x.hat)^2`, where
#' `bias(x.hat)` is an estimate of `bias(x.hat,x)`. Sometimes, we
#' can estimate the bias in closed form, but other times simulations must be
#' done, such as bootstrapping the bias. And, for really large samples, since
#' the MLE is asymptotically unbiased, `trace(vcov(x.hat))` may be a
#' reasonable estimate of the MSE.
#'
#' @param x the `mle` object to compute the MSE of.
#' @param theta true parameter value
#' @export
mse.mle <- function(x, theta = NULL) {
    sum(diag(vcov(x))) + sum(bias(x, theta)^2)
}

#' Computes the point estimate of an `mle` object.
#'
#' @param x the `mle` object.
#' @param ... pass additional arguments
#' @export
point.mle <- function(x, ...) {
    x$theta.hat
}

#' Function for obtaining the fisher information matrix of an `mle` object.
#'
#' @param x the `mle` object to obtain the fisher information of.
#' @param ... pass additional arguments
#' @export
fim.mle <- function(x, ...) {
    x$info
}

#' Function for obtaining a summary of `object`, which is a fitted
#' `mle` object.
#'
#' @param object the `mle` object
#' @param ... pass additional arguments
#'
#' @export
summary.mle <- function(object, ...) {
    structure(list(x = object),
        class = c("summary_mle", "summary")
    )
}

#' Function for printing a `summary` object for an `mle` object.
#'
#' @param object the `summary_mle` object
#' @param ... pass additional arguments
#'
#' @export
print.summary_mle <- function(object, ...) {
    cat("Maximum likelihood estimator of type", class(object$x)[1], "is normally distributed.\n")
    cat("The estimates of the parameters are given by:\n")
    print(point(object$x))

    SE <- se(object$x)
    if (!is.null(SE))
    {
        cat("The standard error is ", SE, ".\n")        
        cat("The asymptotic 95% confidence interval of the parameters are given by:\n")
        print(confint(object$x))
    }
    cat("The MSE of the estimator is ", mse(object$x), ".\n")
    if (!is.null(loglike(object$x))) {
        cat("The log-likelihood is ", loglike(object$x), ".\n")
        cat("The AIC is ", aic(object$x), ".\n")
    }
}

#' Function for obtaining an estimate of the standard error of the MLE
#' `object`.
#'
#' @param object the MLE object
#' @param ... pass additional arguments
#' @export
se.mle <- function(object, ...) {
    V <- vcov(object, ...)
    if (is.null(V)) return(NULL)
    sqrt(diag(as.matrix(V)))
}

#' Determine if an object `x` is an `mle` object.
#'
#' @param x the object to test
#' @export
is_mle <- function(x) {
    inherits(x, "mle")
}


#' Method for determining the orthogonal components of an `mle` object
#' `x`.
#'
#' @param x the `mle` object
#' @param tol the tolerance for determining if a number is close enough to zero
#' @param ... pass additional arguments
#'
#' @export
orthogonal.mle <- function(x, tol = sqrt(.Machine$double.eps), ...) {
    I <- fim(x, ...)
    if(is.null(I)) {
        return(NULL)
    }

    abs(fim(x, ...)) <= tol
}

#' score
#' 
#' Computes the score of an `mle` object.
#'
#' If reguarlity conditions are satisfied, it should be zero (or approximately,
#' if rounding errors occur).
#'
#' @inheritParams score
#' @method score mle
#' @export
score.mle <- function(x, ...) {
    x$score
}

#' Computes the bias of an `mle` object assuming the large sample
#' approximation is valid and the MLE regularity conditions are satisfied.
#' In this case, the bias is zero (or zero vector).
#' 
#' This is not a good estimate of the bias in general, but it's
#' arguably better than returning `NULL`.
#'
#' @param x the `mle` object to compute the bias of.
#' @param par true parameter value. normally, unknown (NULL), in which case
#'              we estimate the bias (say, using bootstrap)
#' @param ... additional arguments to pass
#'
#' @export
bias.mle <- function(x, par = NULL, ...) {
    rep(0, nparams(x, ...))
}

#' Estimate of predictive interval of `T|x` where `T` is a statistic and `x`
#' is an `mle` object that it is conditioned on. Generally, what this means
#' is that we want to find the distribution of `T` when we marginalize over
#' the distribution of `x`, since `x` is an uncertain estimate of the true
#' parameter value.
#'
#' Let
#'   `T|x ~ f(t|x)``
#' be the pdf of vector `T` given MLE `x` and
#'   `x ~ MVN(point(x),vcov(x))``
#' be the multivariate normal of the MLE `x`. Then,
#'   `(T,x) ~ f(t,x) = f(t|x) f(x)
#' is the joint distribution of `(T,x)`. To find `f(t)`,
#' we integrate `f(t,theta.hat)` over `x`.
#'
#' @param x `mle` object
#' @param samp `mle` sampler for parametric distribution that is a function
#'             of `x`
#' @param alpha (1-alpha)-predictive interval for `T|x`
#' @importFrom mvtnorm rmvnorm
#' @export
pred.mle <- function(x, samp, alpha = 0.05) {
    N <- 50000
    theta <- point(x)
    names(theta) <- NULL # maybe more efficient? benchmark this
    sigma <- vcov(x)
    thetas <- NULL
    if (length(theta) == 1) {
        thetas <- matrix(rnorm(N, theta, sd = sqrt(sigma)),nrow=N)
    } else {
        thetas <- rmvnorm(N, theta, sigma)
    }

    ob <- samp(1, theta)
    p <- length(ob)

    data <- matrix(nrow = N, ncol = p)
    data[1,] <- ob
    for (i in 2:N) {
        data[i,] <- samp(1, thetas[i,])
    }

    pi <- matrix(nrow=p, ncol=3)
    for (j in 1:p)
    {
        ## predictive interval for T_j|x
        dj <- sort(data[,j])
        pi[j,] <- c(mean(dj),quantile(dj, c(alpha / 2, 1 - alpha / 2)))
    }
    colnames(pi) <- c("mean", "lower", "upper")
    pi
}


