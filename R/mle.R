#' Constructor for making `mle_fit` objects, which provides a common interface
#' for maximum likelihood estimators.
#'
#' This MLE makes the asymptotic assumption by default. Other MLEs,
#' like `mle_fit_boot`, may not make this assumption.
#'
#' @param theta.hat the MLE
#' @param loglike the log-likelihood of `theta.hat` given the data
#' @param score the score function evaluated at `theta.hat`
#' @param sigma the variance-covariance matrix of `theta.hat` given that data
#' @param info the information matrix of `theta.hat` given the data
#' @param obs observation (sample) data
#' @param nobs number of observations in `obs`
#' @param superclasses class (or classes) with `mle_fit` as base
#' @return An object of class \code{mle_fit}.
#' @examples
#' # MLE for normal distribution (mean and variance)
#' set.seed(123)
#' x <- rnorm(100, mean = 5, sd = 2)
#' n <- length(x)
#' mu_hat <- mean(x)
#' var_hat <- mean((x - mu_hat)^2)  # MLE of variance
#'
#' # Asymptotic variance-covariance of MLE
#' # For normal: Var(mu_hat) = sigma^2/n, Var(var_hat) = 2*sigma^4/n
#' sigma_matrix <- diag(c(var_hat/n, 2*var_hat^2/n))
#'
#' fit <- mle(
#'   theta.hat = c(mu = mu_hat, var = var_hat),
#'   sigma = sigma_matrix,
#'   loglike = sum(dnorm(x, mu_hat, sqrt(var_hat), log = TRUE)),
#'   nobs = n
#' )
#'
#' params(fit)
#' vcov(fit)
#' confint(fit)
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
        class = unique(c(superclasses, "mle_fit"))
    )
}

#' Print method for `mle_fit` objects.
#'
#' @param x the `mle_fit` object to print
#' @param ... additional arguments to pass
#' @return Invisibly returns \code{x}.
#' @export
print.mle_fit <- function(x, ...) {
    print(summary(x, ...))
}

#' Computes the variance-covariance matrix of `mle_fit` object.
#'
#' @param object the `mle_fit` object to obtain the variance-covariance of
#' @param ... additional arguments to pass (not used)
#' @return the variance-covariance matrix
#' @export
vcov.mle_fit <- function(object, ...) {
    object$sigma
}


#' Method for obtaining the parameters of an `mle_fit` object.
#'
#' @param x the `mle_fit` object to obtain the parameters of
#' @return Numeric vector of parameter estimates.
#' @export
params.mle_fit <- function(x) {
    x$theta.hat
}


#' Method for obtaining the number of parameters of an `mle_fit` object.
#'
#' @param x the `mle_fit` object to obtain the number of parameters of
#' @return Integer number of parameters.
#' @importFrom algebraic.dist params nparams
#' @export
nparams.mle_fit <- function(x) {
    length(params(x))
}

#' Log-likelihood of an \code{mle_fit} object.
#'
#' Returns a \code{"logLik"} object with \code{df} (number of parameters) and
#' \code{nobs} attributes, enabling \code{AIC()} and \code{BIC()} via
#' \code{stats::AIC.default}.
#'
#' @param object the \code{mle_fit} object
#' @param ... additional arguments (not used)
#' @return A \code{"logLik"} object.
#' @importFrom stats logLik
#' @export
logLik.mle_fit <- function(object, ...) {
    val <- object$loglike
    if (is.null(val)) val <- NA_real_
    structure(val,
              df = nparams(object),
              nobs = nobs(object),
              class = "logLik")
}

#' Extract coefficients from an \code{mle_fit} object.
#'
#' Delegates to \code{\link[algebraic.dist]{params}()}.
#'
#' @param object the \code{mle_fit} object
#' @param ... additional arguments (not used)
#' @return Named numeric vector of parameter estimates.
#' @importFrom stats coef
#' @export
coef.mle_fit <- function(object, ...) {
    params(object)
}

#' Method for obtaining the number of observations in the sample used by
#' an `mle_fit`.
#'
#' @param object the `mle_fit` object to obtain the number of observations for
#' @param ... additional arguments to pass (not used)
#' @return Integer number of observations, or NULL if not available.
#' @importFrom stats nobs
#' @export
nobs.mle_fit <- function(object, ...) {
    object$nobs
}

#' Method for obtaining the observations used by the `mle_fit` object `x`.
#'
#' @param x the `mle_fit` object to obtain the number of observations for
#' @return The observation data used to fit the MLE, or NULL if not stored.
#' @export
obs.mle_fit <- function(x) {
    x$obs
}

#' Function to compute the confidence intervals of `mle_fit` objects.
#'
#' @param object the `mle_fit` object to compute the confidence intervals for
#' @param parm the parameters to compute the confidence intervals for (not used)
#' @param level confidence level, defaults to 0.95 (alpha=.05)
#' @param use_t_dist logical, whether to use the t-distribution to compute
#'                   the confidence intervals.
#' @param ... additional arguments to pass
#' @return Matrix of confidence intervals with columns for lower and upper bounds.
#' @importFrom stats qt qnorm vcov nobs
#' @importFrom algebraic.dist params
#' @export
confint.mle_fit <- function(object, parm = NULL, level = .95,
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
    se <- sqrt(sigma)
    ci <- cbind(theta - q * se, theta + q * se)
    colnames(ci) <- c(paste0(alpha * 100, "%"),
                      paste0((1 - alpha) * 100, "%"))

    if (is.null(names(theta))) {
        rownames(ci) <- paste0("param", seq_len(p))
    } else {
        rownames(ci) <- names(theta)[seq_len(p)]
    }
    ci
}

#' Method for sampling from an `mle_fit` object.
#'
#' It creates a sampler for the `mle_fit` object. It returns a function
#' that accepts a single parameter `n` denoting the number of samples
#' to draw from the `mle_fit` object.
#'
#' @param x the `mle_fit` object to create sampler for
#' @param ... additional arguments to pass
#' @return A function that takes parameter \code{n} and returns \code{n} samples
#'   from the asymptotic distribution of the MLE.
#' @importFrom stats rnorm vcov
#' @importFrom mvtnorm rmvnorm
#' @importFrom algebraic.dist params nparams
#' @export
sampler.mle_fit <- function(x, ...) {
    sampler(mle_to_dist(x), ...)
}

#' Computes the MSE of an `mle_fit` object.
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
#' @param x the `mle_fit` object to compute the MSE of.
#' @param theta true parameter value, defaults to `NULL` for unknown. If `NULL`,
#'             then we let the bias method deal with it. Maybe it has a nice way
#'             of estimating the bias.
#' @return The MSE as a scalar (univariate) or matrix (multivariate).
#' @importFrom algebraic.dist nparams
#' @importFrom stats vcov
#' @export
mse.mle_fit <- function(x, theta = NULL) {

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

#' Function for obtaining the observed FIM of an `mle_fit` object.
#'
#' @param x the `mle_fit` object to obtain the FIM of.
#' @param ... pass additional arguments
#' @return The observed Fisher Information Matrix, or NULL if not available.
#' @export
observed_fim.mle_fit <- function(x, ...) {
    x$info
}

#' Function for obtaining a summary of `object`, which is a fitted
#' `mle_fit` object.
#'
#' @param object the `mle_fit` object
#' @param ... pass additional arguments
#' @return An object of class \code{summary_mle_fit}.
#' @export
summary.mle_fit <- function(object, ...) {
    structure(list(x = object), class = c("summary_mle_fit", "summary"))
}

#' Function for printing a `summary` object for an `mle_fit` object.
#'
#' @param x the `summary_mle_fit` object
#' @param ... pass additional arguments
#' @return Invisibly returns \code{x}.
#' @importFrom stats vcov
#' @importFrom algebraic.dist nparams params
#' @export
print.summary_mle_fit <- function(x, ...) {
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
    ll <- logLik(x$x)
    if (!is.na(as.numeric(ll))) {
        cat("The log-likelihood is ", as.numeric(ll), ".\n")
        cat("The AIC is ", AIC(x$x), ".\n")
    }
}

#' Function for obtaining an estimate of the standard error of the MLE
#' object `x`.
#'
#' @param x the MLE object
#' @param se.matrix if `TRUE`, return the square root of the variance-covariance
#' @param ... additional arguments to pass (not used)
#' @return Vector of standard errors, or matrix if \code{se.matrix = TRUE}, or
#'   NULL if variance-covariance is not available.
#' @importFrom stats vcov
#' @export
se.mle_fit <- function(x, se.matrix = FALSE, ...) {

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

#' Determine if an object `x` is an `mle_fit` object.
#'
#' @param x the object to test
#' @return Logical TRUE if \code{x} is an \code{mle_fit} object, FALSE otherwise.
#' @examples
#' fit <- mle(theta.hat = c(mu = 5), sigma = matrix(0.1))
#' is_mle(fit)       # TRUE
#' is_mle(list(a=1)) # FALSE
#' @export
is_mle <- function(x) {
    inherits(x, "mle_fit")
}

#' Method for determining the orthogonal components of an `mle_fit` object
#' `x`.
#'
#' @param x the `mle_fit` object
#' @param tol the tolerance for determining if a number is close enough to zero
#' @param ... pass additional arguments
#' @return Logical matrix indicating which off-diagonal FIM elements are
#'   approximately zero (orthogonal parameters), or NULL if FIM unavailable.
#' @export
orthogonal.mle_fit <- function(x, tol = sqrt(.Machine$double.eps), ...) {
    I <- observed_fim(x, ...)
    if(is.null(I)) {
        return(NULL)
    }

    abs(I) <= tol
}

#' Computes the score of an `mle_fit` object (score evaluated at the MLE).
#'
#' If reguarlity conditions are satisfied, it should be zero (or approximately,
#' if rounding errors occur).
#'
#' @param x the `mle_fit` object to compute the score of.
#' @param ... additional arguments to pass (not used)
#' @return The score vector evaluated at the MLE, or NULL if not available.
#' @export
score_val.mle_fit <- function(x, ...) {
    x$score
}

#' Computes the bias of an `mle_fit` object assuming the large sample
#' approximation is valid and the MLE regularity conditions are satisfied.
#' In this case, the bias is zero (or zero vector).
#'
#' This is not a good estimate of the bias in general, but it's
#' arguably better than returning `NULL`.
#'
#' @param x the `mle_fit` object to compute the bias of.
#' @param theta true parameter value. normally, unknown (NULL), in which case
#'              we estimate the bias (say, using bootstrap)
#' @param ... additional arguments to pass
#' @return Numeric vector of zeros (asymptotic bias is zero under regularity
#'   conditions).
#' @importFrom algebraic.dist nparams
#' @export
bias.mle_fit <- function(x, theta = NULL, ...) {
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
#' @param x an `mle_fit` object.
#' @param samp The sampler for the distribution that is parameterized by the
#'             MLE `x`, i.e., `T|x`.
#' @param alpha (1-alpha)-predictive interval for `T|x`. Defaults to 0.05.
#' @param R number of samples to draw from the sampling distribution of `x`.
#'          Defaults to 50000.
#' @param ... additional arguments to pass into `samp`.
#' @return Matrix with columns for mean, lower, and upper bounds of the
#'   predictive interval.
#' @importFrom mvtnorm rmvnorm
#' @importFrom stats vcov quantile rnorm
#' @importFrom algebraic.dist params
#' @export
pred.mle_fit <- function(x, samp, alpha = 0.05, R = 50000, ...) {

    theta <- params(x)
    thetas <- as.matrix(sampler(x)(R), nrow = R)

    ob <- samp(n = 1, theta, ...)
    p <- length(ob)

    data <- matrix(NA, nrow = R, ncol = p)
    data[1, ] <- ob
    for (i in seq_len(R)[-1]) {
        data[i, ] <- samp(1, thetas[i, ], ...)
    }

    PI <- t(apply(data, 2, function(dj) {
        c(mean(dj), quantile(x = dj, probs = c(alpha / 2, 1 - alpha / 2)))
    }))
    colnames(PI) <- c("mean", "lower", "upper")
    PI
}

#' Expectation operator applied to `x` of type `mle_fit`
#' with respect to a function `g`. That is, `E(g(x))`.
#' 
#' Optionally, we use the CLT to construct a CI(`alpha`) for the
#' estimate of the expectation. That is, we estimate `E(g(x))` with
#' the sample mean and Var(g(x)) with the sigma^2/n, where sigma^2
#' is the sample variance of g(x) and n is the number of samples.
#' From these, we construct the CI.
#'
#' @param x `mle_fit` object
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
expectation.mle_fit <- function(
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
#' @return An \code{mle_fit} object representing the marginal distribution for the
#'   selected parameter indices.
#' @importFrom algebraic.dist marginal params obs
#' @importFrom stats vcov nobs
#' @export
marginal.mle_fit <- function(x, indices) {
    if (length(indices) == 0) {
        stop("indices must be non-empty")
    }
    else if (any(indices < 1) || any(indices > nparams(x))) {
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

## ── Distribution methods (delegate to asymptotic normal/mvn) ──────────────

#' PDF of the asymptotic distribution of an MLE.
#'
#' Returns a closure computing the probability density function of the
#' asymptotic normal (univariate) or multivariate normal (multivariate)
#' distribution of the MLE.
#'
#' @param x An \code{mle_fit} object.
#' @param ... Additional arguments (not used).
#' @return A function computing the PDF at given points.
#' @importFrom stats density
#' @importFrom algebraic.dist normal mvn
#' @export
density.mle_fit <- function(x, ...) {
    density(mle_to_dist(x))
}

#' CDF of the asymptotic distribution of an MLE.
#'
#' @param x An \code{mle_fit} object.
#' @param ... Additional arguments (not used).
#' @return A function computing the CDF at given points.
#' @importFrom algebraic.dist cdf
#' @export
cdf.mle_fit <- function(x, ...) {
    cdf(mle_to_dist(x))
}

#' Quantile function of the asymptotic distribution of an MLE.
#'
#' @param x An \code{mle_fit} object.
#' @param ... Additional arguments (not used).
#' @return A function computing quantiles for given probabilities.
#' @importFrom algebraic.dist inv_cdf
#' @export
inv_cdf.mle_fit <- function(x, ...) {
    inv_cdf(mle_to_dist(x))
}

#' Support of the asymptotic distribution of an MLE.
#'
#' @param x An \code{mle_fit} object.
#' @return A support object (interval for normal distributions).
#' @importFrom algebraic.dist sup
#' @export
sup.mle_fit <- function(x) {
    sup(mle_to_dist(x))
}

#' Dimension (number of parameters) of an MLE.
#'
#' @param x An \code{mle_fit} object.
#' @return Integer; the number of parameters.
#' @export
dim.mle_fit <- function(x) {
    nparams(x)
}

#' Mean of the asymptotic distribution of an MLE.
#'
#' Returns the parameter estimates, which are the mean of the asymptotic
#' distribution.
#'
#' @param x An \code{mle_fit} object.
#' @param ... Additional arguments (not used).
#' @return Numeric vector of parameter estimates.
#' @export
mean.mle_fit <- function(x, ...) {
    params(x)
}

#' Conditional distribution from an MLE.
#'
#' Computes the conditional distribution of a subset of parameters given
#' observed values for other parameters. Uses the closed-form Schur complement
#' for the multivariate normal. Returns a \code{dist} object (not an \code{mle_fit}),
#' since conditioning loses the MLE context.
#'
#' @param x An \code{mle_fit} object with at least 2 parameters.
#' @param P Optional predicate function for Monte Carlo fallback.
#' @param ... Additional arguments forwarded to \code{P}.
#' @param given_indices Integer vector of conditioned parameter indices.
#' @param given_values Numeric vector of observed values.
#' @return A \code{normal}, \code{mvn}, or \code{empirical_dist} object.
#' @importFrom algebraic.dist conditional
#' @export
conditional.mle_fit <- function(x, P = NULL, ...,
                            given_indices = NULL, given_values = NULL) {
    conditional(mle_to_dist(x), P = P, ...,
                given_indices = given_indices, given_values = given_values)
}