#' Constructor for making `mle` objects, which provides a common interface
#' for maximum likelihood estimators.
#' 
#' This MLE makes the asymptotic assumption by default. Other MLEs,
#' like `mle_boot`, may not make this assumption.
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

#' Computes the variance-covariance matrix of `mle` object.
#' 
#' @param object the `mle` object to obtain the variance-covariance of
#' @param ... additional arguments to pass (not used)
#' @return the variance-covariance matrix
#' @export
vcov.mle <- function(object, ...) {
    object$sigma
}


#' Method for obtaining the parameters of an `mle` object.
#' 
#' @param x the `mle` object to obtain the parameters of
#' @export
params.mle <- function(x) {
    x$theta.hat
}

#' Method for obtaining the number of parameters of an `mle` object.
#' 
#' @param x the `mle` object to obtain the number of parameters of
#' @importFrom algebraic.dist params nparams
#' @export
nparams.mle <- function(x) { 
    length(params(x))
}

#' Method for obtaining the AIC of an `mle` object.
#'
#' @param x the `mle` object to obtain the AIC of
#' @importFrom algebraic.dist nparams
#' @export
aic.mle <- function(x) {
    -2 * loglik_val(x) + 2 * nparams(x)
}

#' Method for obtaining the number of observations in the sample used by
#' an `mle`.
#'
#' @param object the `mle` object to obtain the number of observations for
#' @param ... additional arguments to pass (not used)
#' @importFrom stats nobs
#' @export
nobs.mle <- function(object, ...) {
    object$nobs
}

#' Method for obtaining the observations used by the `mle` object `x`.
#'
#' @param x the `mle` object to obtain the number of observations for
#' @export
obs.mle <- function(x) {
    x$obs
}

#' Method for obtaining the log-likelihood of an `mle` object.
#'
#' @param x the log-likelihood `l` evaluated at `x`, `l(x)`.
#' @param ... additional arguments to pass
#' @return the log-likelihood of the fitted mle object `x`
#' @export
loglik_val.mle <- function(x, ...) {
    x$loglike
}

#' Function to compute the confidence intervals of `mle` objects.
#'
#' @param object the `mle` object to compute the confidence intervals for
#' @param parm the parameters to compute the confidence intervals for (not used)
#' @param level confidence level, defaults to 0.95 (alpha=.05)
#' @param use_t_dist logical, whether to use the t-distribution to compute
#'                   the confidence intervals.
#' @param ... additional arguments to pass
#' @importFrom stats qt qnorm vcov nobs
#' @importFrom algebraic.dist params
#' @export
confint.mle <- function(object, parm = NULL, level = .95,
    use_t_dist = FALSE, ...) {
    stopifnot(is.numeric(level), level >= 0, level <= 1)
    V <- vcov(object)
    if (is.null(V)) stop("No variance-covariance matrix available.")
    if (is.matrix(V)) {
        V <- diag(V)
    }
    sigma <- V

    theta <- params(object)
    alpha <- (1 - level) / 2
    p <- length(theta)

    if (use_t_dist && is.null(nobs(object))) {
        warning("Unknown number of observations, using large sample approximation.")
        use_t_dist <- FALSE
    }

    if (use_t_dist) {
        q <- qt(1 - alpha, df = nobs(object) - 1)
    } else {
        q <- qnorm(1 - alpha)
    }
    ci <- matrix(nrow = p, ncol = 2)
    colnames(ci) <- c(paste0(alpha * 100, "%"),
                      paste0((1 - alpha) * 100, "%"))

    for (j in 1:p) {
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
#' @importFrom stats rnorm vcov
#' @importFrom mvtnorm rmvnorm
#' @importFrom algebraic.dist params nparams
#' @export
sampler.mle <- function(x, ...) {
    V <- vcov(x)
    mu <- params(x)

    if (nparams(x) == 1L) {
        stddev <- sqrt(V)
        function(n = 1) rnorm(n, mean = mu, sd = stddev, ...)
    } else {
        function(n = 1) rmvnorm(n, mu, V, ...)
    }
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
#' @importFrom algebraic.dist nparams
#' @importFrom stats vcov
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

#' Function for obtaining the observed FIM of an `mle` object.
#'
#' @param x the `mle` object to obtain the FIM of.
#' @param ... pass additional arguments
#' @export
observed_fim.mle <- function(x, ...) {
    x$info
}

#' Function for obtaining a summary of `object`, which is a fitted
#' `mle` object.
#'
#' @param object the `mle` object
#' @param ... pass additional arguments
#' @export
summary.mle <- function(object, ...) {
    structure(list(x = object), class = c("summary_mle", "summary"))
}

#' Function for printing a `summary` object for an `mle` object.
#'
#' @param x the `summary_mle` object
#' @param ... pass additional arguments
#' @importFrom stats vcov
#' @importFrom algebraic.dist nparams params
#' @export
print.summary_mle <- function(x, ...) {
    cat("Maximum likelihood estimator of type", class(x$x)[1], "is normally distributed.\n")
    cat("The estimates of the parameters are given by:\n")
    print(params(x$x))

    SE <- se(x$x)
    if (!is.null(SE))
    {
        cat("The standard error is ", SE, ".\n")        
        cat("The asymptotic 95% confidence interval of the parameters are given by:\n")
        print(confint(x$x))
    }
    if (nparams(x$x) == 1L) {
        cat("The MSE of the estimator is ", mse(x$x), ".\n")
    } else {
        cat("The MSE of the individual components in a multivariate estimator is:\n")
        print(mse(x$x))
    }
    if (!is.null(loglik_val(x$x))) {
        cat("The log-likelihood is ", loglik_val(x$x), ".\n")
        cat("The AIC is ", aic(x$x), ".\n")
    }
}

#' Function for obtaining an estimate of the standard error of the MLE
#' object `x`.
#' 
#' @param x the MLE object
#' @param se.matrix if `TRUE`, return the square root of the variance-covariance
#' @param ... additional arguments to pass (not used)
#' @importFrom stats vcov
#' @export
se.mle <- function(x, se.matrix = FALSE, ...) {

    V <- vcov(x)
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
#' @export
orthogonal.mle <- function(x, tol = sqrt(.Machine$double.eps), ...) {
    I <- observed_fim(x, ...)
    if(is.null(I)) {
        return(NULL)
    }

    abs(I) <= tol
}

#' Computes the score of an `mle` object (score evaluated at the MLE).
#'
#' If reguarlity conditions are satisfied, it should be zero (or approximately,
#' if rounding errors occur).
#'
#' @param x the `mle` object to compute the score of.
#' @param ... additional arguments to pass (not used)
#' @export
score_val.mle <- function(x, ...) {
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
#' @importFrom algebraic.dist nparams
#' @export
bias.mle <- function(x, theta = NULL, ...) {
    rep(0, nparams(x))
}

#' Estimate of predictive interval of `T|data` using Monte Carlo integration.
#'
#' Let
#'   `T|x ~ f(t|x)``
#' be the pdf of vector `T` given MLE `x` and
#'   `x ~ MVN(params(x),vcov(x))``
#' be the estimate of the sampling distribution of the MLE for the parameters
#' of `T`. Then,
#'   `(T,x) ~ f(t,x) = f(t|x) f(x)
#' is the joint distribution of `(T,x)`. To find `f(t)` for a fixed `t`, we
#' integrate `f(t,x)` over `x` using Monte Carlo integration to find the
#' marginal distribution of `T`. That is, we:
#' 
#' 1. Sample from MVN `x`
#' 2. Compute `f(t,x)` for each sample
#' 3. Take the mean of the `f(t,x)` values asn an estimate of `f(t)`.
#' 
#' The `samp` function is used to sample from the distribution of `T|x`. It should
#' be designed to take 
#'
#' @param x an `mle` object.
#' @param samp The sampler for the distribution that is parameterized by the
#'             MLE `x`, i.e., `T|x`.
#' @param alpha (1-alpha)-predictive interval for `T|x`. Defaults to 0.05.
#' @param R number of samples to draw from the sampling distribution of `x`.
#'          Defaults to 50000.
#' @param ... additional arguments to pass into `samp`.
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats vcov quantile rnorm
#' @importFrom algebraic.dist params
#' @export
pred.mle <- function(x, samp, alpha = 0.05, R = 50000, ...) {

    theta <- params(x)
    sigma <- vcov(x)

    thetas <- NULL
    if (length(theta) == 1) {
        thetas <- matrix(rnorm(n = R, mean = theta, sd = sqrt(sigma)), nrow = R)
    } else {
        thetas <- rmvnorm(n = R, mean = theta, sigma = sigma)
    }

    ob <- samp(n = 1, theta, ...)
    p <- length(ob)

    data <- matrix(NA, nrow = R, ncol = p)
    data[1, ] <- ob
    for (i in 2:R) {
        data[i, ] <- samp(1, thetas[i, ], ...)
    }

    PI <- matrix(NA, nrow = p, ncol = 3)
    for (j in 1:p) {
        dj <- data[, j]
        PI[j, ] <- c(mean(dj), quantile(
            x = dj, probs = c(alpha / 2, 1 - alpha / 2)))
    }
    colnames(PI) <- c("mean", "lower", "upper")
    PI
}

#' Expectation operator applied to `x` of type `mle`
#' with respect to a function `g`. That is, `E(g(x))`.
#' 
#' Optionally, we use the CLT to construct a CI(`alpha`) for the
#' estimate of the expectation. That is, we estimate `E(g(x))` with
#' the sample mean and Var(g(x)) with the sigma^2/n, where sigma^2
#' is the sample variance of g(x) and n is the number of samples.
#' From these, we construct the CI.
#'
#' @param x `mle` object
#' @param g characteristic function of interest, defaults to identity
#' @param ... additional arguments to pass to `g`
#' @param control a list of control parameters:
#'  compute_stats - Whether to compute CIs for the expectations, defaults
#'                  to FALSE
#'  n             - The number of samples to use for the MC estimate,
#'                  defaults to 10000
#'  alpha         - The significance level for the confidence interval,
#'                  defaults to 0.05
#' @return If `compute_stats` is FALSE, then the estimate of the expectation,
#'        otherwise a list with the following components:
#'   value - The estimate of the expectation
#'   ci    - The confidence intervals for each component of the expectation
#'   n     - The number of samples
#' @importFrom algebraic.dist expectation_data expectation
#' @importFrom utils modifyList
#' @export
expectation.mle <- function(
    x,
    g = function(t) t,
    ...,
    control = list()) {

    defaults <- list(
      compute_stats = FALSE,
      n = 10000L,
      alpha = 0.05)
    control <- modifyList(defaults, control)

    stopifnot(is.numeric(control$n), control$n > 0,
      is.numeric(control$alpha), control$alpha > 0, control$alpha < 1)
  
    expectation_data(
      data = sampler(x)(control$n),
      g = g,
      ...,
      compute_stats = control$compute_stats,
      alpha = control$alpha)
  }

#' Method for obtaining the marginal distribution of an MLE
#' that is based on asymptotic assumptions:
#' 
#'  `x ~ MVN(params(x), inv(H)(x))`
#' 
#' where H is the (observed or expecation) Fisher information
#' matrix.
#' 
#' @param x The distribution object.
#' @param indices The indices of the marginal distribution to obtain.
#' @importFrom algebraic.dist marginal params obs
#' @importFrom stats vcov nobs
#' @export
marginal.mle <- function(x, indices) {
    if (length(indices) == 0) {
        stop("indices must be non-empty")
    }
    else if (any(indices < 0) || any(indices > nparams(x))) {
        stop("indices must be in [1, dim(x)]")
    }
    
    mle(theta.hat = params(x)[indices],
        sigma = vcov(x)[indices, indices],
        loglike = NULL,
        score = NULL,
        info = NULL,
        obs = obs(x),
        nobs = nobs(x))
}