#' \code{mle} makes an \code{mle} object.
#'
#' @param theta.hat the MLE
#' @param loglike the log-likelihood of \code{theta.hat} given the data
#' @param score the score function evaluated at \code{theta.hat}
#' @param sigma the variance-covariance matrix of \code{theta.hat} given that data
#' @param info the information matrix of \code{theta.hat} given the data
#' @param obs observation (sample) data
#' @param nobs number of observations in \code{obs}
#' @param superclasses class (or classes) with \code{mle} as base
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
#' an \code{mle} object \code{x}.
#'
#' @param x the \code{mle} object to print
#' @param ... additional arguments to pass
#' @export
print.mle <- function(x, ...) {
    print(summary(x, ...))
}

#' Method for obtaining the parameters of an \code{mle} object.
#'
#' @param x the \code{mle} object to obtain the parameters of
#'
#' @export
params.mle <- function(x, ...) {
    x$theta
}

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
nobs.mle <- function(object, ...) object$nobs

#' Method for obtaining the observations used by the \code{mle}.
#'
#' @param object the \code{mle} object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @export
obs.mle <- function(object, ...) object$obs

#' Method for obtaining the log-likelihood of an \code{mle} object.
#'
#' @param x the log-likelihood \code{l} evaluated at \code{x}, \code{l(x)}.
#' @param ... additional arguments to pass
#'
#' @export
loglike.mle <- function(x, ...) x$loglike

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
confint.mle <- function(object, parm = NULL, level = .95, ...) {
    V <- vcov(object, ...)
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

#' Method for sampling from an \code{mle} object.
#'
#' It creates a sampler for the \code{mle} object. It returns a function
#' that accepts a single parameter \code{n} denoting the number of samples
#' to draw from the \code{mle} object.
#'
#' @param x the \code{mle} object to create sampler for
#' @param ... additional arguments to pass
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


#' Computes the variance-covariance matrix of \code{mle} objects.
#'
#' @param object the \code{mle} object to obtain the variance-covariance of
#' @param ... additional arguments to pass
#'
#' @importFrom stats vcov
#' @export
vcov.mle <- function(object, ...) {
    object$sigma
}


#' Computes the MSE of an \code{mle} object.
#'
#' The MSE of an estimator is just the expected sum of squared differences,
#' e.g., if the true parameter value is \code{x} and we have an estimator \code{x.hat},
#' then the MSE is
#' \code{mse(x.hat) = E[(x.hat-x)^2] = trace(vcov(x.hat)) + (bias(x.hat))^2}.
#'
#' Since \code{x} is not typically known, we normally must estimate the bias.
#' For sufficiently large samples, for the MLE assuming the regularity conditions,
#' the bias is given by
#' \code{mse(x.hat) = trace(vcov(x.hat)) + bias(x.hat)^2}, where
#' \code{bias(x.hat)} is an estimate of \code{bias(x.hat,x)}. Sometimes, we
#' can estimate the bias in closed form, but other times simulations must be
#' done, such as bootstrapping the bias. And, for really large samples, since
#' the MLE is asymptotically unbiased, \code{trace(vcov(x.hat))} may be a
#' reasonable estimate of the MSE.
#'
#' @param x the \code{mle} object to compute the MSE of.
#' @param theta true parameter value
#' @export
mse.mle <- function(x, theta = NULL) {
    sum(diag(vcov(x))) + sum(bias(x, theta)^2)
}

#' Computes the point estimate of an \code{mle} object.
#'
#' @param x the \code{mle} object.
#' @param ... pass additional arguments
#' @export
point.mle <- function(x, ...) {
    x$theta.hat
}

#' Function for obtaining the fisher information matrix of an \code{mle} object.
#'
#' @param x the \code{mle} object to obtain the fisher information of.
#' @param ... pass additional arguments
#' @export
fim.mle <- function(x, ...) {
    x$info
}

#' Function for obtaining a summary of \code{object}, which is a fitted
#' \code{mle} object.
#'
#' @param object the \code{mle} object
#' @param ... pass additional arguments
#'
#' @export
summary.mle <- function(object, ...) {
    structure(list(x = object),
        class = c("summary_mle", "summary")
    )
}

#' Function for printing a \code{summary} object for an \code{mle} object.
#'
#' @param object the \code{summary_mle} object
#' @param ... pass additional arguments
#'
#' @export
print.summary_mle <- function(object, ...) {
    cat("Maximum likelihood estimator of type", class(object$x)[1], "is normally distributed.\n")
    cat("The estimates of the parameters are given by:\n")
    print(point(object$x))
    cat("The fisher information matrix (FIM) is given by:\n")
    print(fim(object$x))
    cat("The variance-covariance matrix of the estimator is given by:\n")
    print(vcov(object$x))
    cat("The asymptotic 95% confidence interval of the parameters are given by:\n")
    print(confint(object$x))
    cat("The bias of the estimator is given by:\n")
    print(bias(object$x))
    cat("The MSE of the estimator is ", mse(object$x), ".\n")
    cat("The log-likelihood is ", loglike(object$x), ".\n")
    cat("The AIC is ", aic(object$x), ".\n")
    cat("The standard error is ", se(object$x), ".\n")
}

#' Function for obtaining an estimate of the standard error of the MLE
#' \code{object}.
#'
#' @param object the MLE object
#' @param ... pass additional arguments
#' @export
se.mle <- function(object, ...) {
    sqrt(diag(as.matrix(vcov(object, ...))))
}

#' Determine if an object \code{x} is an \code{mle} object.
#'
#' @param x the object to test
#' @export
is_mle <- function(x) {
    inherits(x, "mle")
}

#' Function for obtaining sample points for an \code{mle} object that is within
#' the \code{p}-probability region.the number of observations in the sample used by
#' an MLE object \code{x}.
#'
#' @param n the sample size
#' @param x the \code{mle} object
#' @param p the probability region
#'
#' @importFrom stats qchisq
#' @importFrom stats mahalanobis
#' @export
sample_mle_region <- function(n, x, p = .95) {
    stopifnot(p > 0.0 && p <= 1.0)
    stopifnot(n > 1L)
    stopifnot(is_mle(x))

    k <- nparams(x)
    crit <- stats::qchisq(p, k)
    nfo <- fim(x)
    mu <- point(x)

    i <- 0L
    samp <- sampler(x)
    data <- matrix(nrow = n, ncol = k)
    while (i < n) {
        x <- samp(1)
        d <- stats::mahalanobis(x, center = mu, cov = nfo, inverted = T)
        if (d <= crit) {
            i <- i + 1L
            data[i, ] <- x
        }
    }
    data
}

#' Method for determining the orthogonal components of an \code{mle} object
#' \code{x}.
#'
#' @param x the \code{mle} object
#' @param tol the tolerance for determining if a number is close enough to zero
#' @param ... pass additional arguments
#'
#' @export
orthogonal.mle <- function(x, tol = sqrt(.Machine$double.eps), ...) {
    abs(fim(x, ...)) <= tol
}

#' Computes the score of an \code{mle} object.
#'
#' If reguarlity conditions are satisfied, it should be zero (or approximately,
#' if rounding errors occur).
#'
#' @param x the \code{mle} object to compute the score of.
#' @param ... additional arguments to pass
#'
#' @export
score.mle <- function(x, ...) {
    x$score
}

#' Computes the bias of an \code{mle} object assuming the large sample
#' approximation is valid and the MLE regularity conditions are satisfied.
#' In this case, the bias is zero (or zero vector).
#'
#' @param x the \code{mle} object to compute the bias of.
#' @param par true parameter value. normally, unknown (NULL), in which case
#'              we estimate the bias (say, using bootstrap)
#' @param ... additional arguments to pass
#'
#' @export
bias.mle <- function(x, par = NULL, ...) {
    rep(0, nparams(x, ...))
}

#' Predictive pdf of T given mle.hat estimate for the parameters of T
#'
#' Let
#'   T|theta.hat ~ f(theta.hat)
#' be the pdf of T given parameter value theta.hat and
#'   theta.hat ~ normal(mean(theta.hat),vcov(theta.hat))
#' be the mvn of the MLE for theta.hat.
#'
#' Then,
#'   (T,theta.hat) ~ f(t,theta.hat) = f(t|theta.hat) f(theta.hat)
#' is the joint distribution of (T,theta.hat). To find f(t), we integrate
#' f(t,theta.hat) over theta.hat.
#'
#' @param x \code{mle} object
#' @param samp \code{mle} sampler for parametric distribution for the mle \code{x}
#' @param alpha (1-alpha)-predictive interval for T|x
#' @importFrom mvtnorm rmvnorm
#' @export
pred.mle <- function(x, samp, alpha = 0.05) {
    N <- 10000
    theta <- point(x)
    sigma <- vcov(x)

    thetas <- NULL
    if (length(theta) == 1) {
        thetas <- rnorm(N, theta, sd = sqrt(sigma))
    } else {
        thetas <- rmvnorm(N, theta, sigma)
    }

    ts <- sort(sapply(thetas, function(theta) samp(1, theta)))
    # predictive interval for T
    pi <- matrix(c(
        mean(ts),
        quantile(ts, c(alpha / 2, 1 - alpha / 2))
    ), nrow = 1)
    colnames(pi) <- c("mean", "lower", "upper")
    pi
}
