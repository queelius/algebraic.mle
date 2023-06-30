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
#' @param x the `mle` object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @importFrom stats nobs
#' @export
nobs.mle <- function(x, ...) x$nobs

#' Method for obtaining the observations used by the `mle` object `x`.
#'
#' @param x the `mle` object to obtain the number of observations for
#' @param ... additional arguments to pass
#' @export
obs.mle <- function(x, ...) x$obs

#' Method for obtaining the log-likelihood of an `mle` object.
#'
#' @param x the log-likelihood `l` evaluated at `x`, `l(x)`.
#' @param ... additional arguments to pass
#'
#' @export
loglike.mle <- function(x, ...) x$loglike

#' Function to compute the confidence intervals of `mle` objects.
#'
#' @param x the `mle` object to compute the confidence intervals for
#' @param level confidence level, defaults to 0.95 (alpha=.05)
#' @param ... additional arguments to pass
#' @export
confint.mle <- function(x, level = .95, use_t_dist = TRUE, ...) {
    stopifnot(is.numeric(level), level >= 0, level <= 1)
    V <- vcov(x, ...)
    if (is.null(V)) stop("No variance-covariance matrix available.")
    if (is.matrix(V)) {
        V <- diag(V)
    }
    sigma <- V

    theta <- point(x, ...)
    alpha <- (1 - level) / 2
    p <- length(theta)

    if (use_t_dist && is.null(nobs(x))) {
        warning("Unknown number of observations, using large sample approximation.")
        use_t_dist <- FALSE
    }

    if (use_t_dist) {
        q <- stats::qt(1 - alpha, df = nobs(x) - 1)
    } else {
        q <- stats::qnorm(1 - alpha)
    }
    ci <- matrix(nrow = p, ncol = 2)
    colnames(ci) <- c(paste0(alpha * 100, "%"),
                      paste0((1 - alpha) * 100, "%"))

    for (j in 1:p)
    {
        ci[j, ] <- c(
            theta[j] - q * sqrt(sigma[j]),
            theta[j] + q * sqrt(sigma[j]))
    }

    if (is.null(names(theta))) {
        rownames(ci) <- paste0("param", 1:p)
    } else {
        rownames(ci) <- names(theta)[1:p]
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
        stddev <- sqrt(V)
        function(n = 1) stats::rnorm(n, mean = mu, sd = stddev, ...)
    } else {
        function(n = 1) mvtnorm::rmvnorm(n, mu, V, ...)
    }
}


#' Computes the variance-covariance matrix of `mle` object `x`.
#'
#' @param x the `mle` object to obtain the variance-covariance of
#' @param ... additional arguments to pass
#'
#' @importFrom stats vcov
#' @export
vcov.mle <- function(x, ...) {
    x$sigma
}


#' Computes the MSE of an `mle` object.
#'
#' The MSE of an estimator is just the expected sum of squared differences,
#' e.g., if the true parameter value is `x` and we have an estimator `x.hat`,
#' then the MSE is
#' ```
#'     mse(x.hat) = E[(x.hat-x) %*% t(x.hat - x)] =
#'                       vcov(x.hat) + bias(x.hat, x) %*% t(bias(x.hat, x))
#' ```
#'
#' Since `x` is not typically known, we normally must estimate the bias.
#' Asymptotically, assuming the regularity conditions, the bias of an MLE
#' is zero, so we can estimate the MSE as `mse(x.hat) = vcov(x.hat)`, but for
#' small samples, this is not generally the case. If we can estimate the bias, then
#' we can replace the bias with an estimate of the bias.
#' 
#' Sometimes, we can estimate the bias analytically, but if not, we can use something
#' like the bootstrap. For example, if we have a sample of size `n`, we can bootstrap
#' the bias by sampling `n` observations with replacement, computing the MLE, and
#' then computing the difference between the bootstrapped MLE and the MLE. We can
#' repeat this process `B` times, and then average the differences to get an estimate
#' of the bias.
#'
#' @param x the `mle` object to compute the MSE of.
#' @param theta true parameter value, defaults to `NULL` for unknown. If `NULL`,
#'             then we let the bias method deal with it. Maybe it has a nice way
#'             of estimating the bias.
#' @export
mse.mle <- function(x, theta = NULL) {

    if (!is.null(theta)) {
        stopifnot(length(theta) == nparams(x))
    }

    b <- bias(x, theta)
    V <- vcov(x)
    if (nparams(x) == 1L) {
        res <- V + b^2
    } else {
        # multivariate MSE is a matrix whose diagonal elements are
        # the MSEs of the individual parameters. note that we are taking the
        # outer product of b with itself.
        res <- V + b %*% t(b)
    }
    res
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
    if (nparams(object$x) == 1L) {
        cat("The MSE of the estimator is ", mse(object$x), ".\n")
    } else {
        cat("The MSE of the individual componetns in a multivariate estimator is:\n")
        print(diag(mse(object$x)))
    }
    cat("The MSE of the estimator is ", mse(object$x), ".\n")
    if (!is.null(loglike(object$x))) {
        cat("The log-likelihood is ", loglike(object$x), ".\n")
        cat("The AIC is ", aic(object$x), ".\n")
    }
}

#' Function for obtaining an estimate of the standard error of the MLE
#' `x`.

#' @param x the bootstrapped MLE object
#' @param se.matrix if `TRUE`, return the square root of the variance-covariance
#' @param ... additional arguments to pass
#' @export
se.mle <- function(x, se.matrix = FALSE, ...) {

    V <- vcov(x, ...)
    if (is.null(V)) return(NULL)

    if (se.matrix && is.matrix(V)) {
        return(chol(V))
    } else {
        if (is.matrix(V)) {
            return(sqrt(diag(V)))
        } else {
            return(sqrt(V))
        }
    }
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
#' @param theta true parameter value. normally, unknown (NULL), in which case
#'              we estimate the bias (say, using bootstrap)
#' @param ... additional arguments to pass
#' @export
bias.mle <- function(x, theta = NULL, ...) {
    rep(0, nparams(x))
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


