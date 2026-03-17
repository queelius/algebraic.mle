#' Computes the distribution of `g(x)` where `x` is an `mle_fit` object.
#'
#' By the invariance property of the MLE, if `x` is an `mle_fit` object,
#' then under the right conditions, asymptotically, `g(x)` is normally
#' distributed,
#'     g(x) ~ normal(g(params(x)),sigma)
#' where `sigma` is the variance-covariance of `f(x)`
#'
#' We provide two different methods for estimating the
#' variance-covariance of `f(x)`:
#'     method = "delta" -> delta method
#'     method = "mc" -> monte carlo method
#'
#' @param x an `mle_fit` object
#' @param g a deterministic, differentiable function mapping a numeric parameter
#'   vector to a numeric vector. The Jacobian is computed numerically via
#'   \code{numDeriv::jacobian} when \code{method = "delta"}.
#' @param ... additional arguments to pass to the `g` function
#' @param n number of samples to take to estimate distribution of `g(x)` if
#'         `method == "mc"`.
#' @param method method to use to estimate distribution of `g(x)`,
#'               "delta" or "mc".
#' @return An \code{mle_fit} object of class \code{mle_fit_rmap} representing the
#'   transformed MLE with variance estimated by the specified method.
#' @examples
#' # MLE for normal distribution
#' set.seed(123)
#' x <- rnorm(100, mean = 5, sd = 2)
#' n <- length(x)
#' fit <- mle(
#'   theta.hat = c(mu = mean(x), var = var(x)),
#'   sigma = diag(c(var(x)/n, 2*var(x)^2/n)),
#'   nobs = n
#' )
#'
#' # Transform: compute MLE of standard deviation (sqrt of variance)
#' # Using delta method
#' g <- function(theta) sqrt(theta[2])
#' sd_mle <- rmap(fit, g, method = "delta")
#' params(sd_mle)
#' se(sd_mle)
#' @importFrom stats cov vcov nobs
#' @importFrom numDeriv jacobian
#' @importFrom algebraic.dist rmap sampler params
#' @importFrom MASS ginv
#' @export
rmap.mle_fit <- function(x, g, ...,
    n = 1000, method = c("mc", "delta")) {
    stopifnot(is.numeric(n), n > 0, is_mle(x), is.function(g))
    n <- as.integer(n)

    method <- match.arg(method)
    theta <- unname(params(x))

    if (method == "mc") {
        samp <- as.matrix(sampler(x)(n), nrow = n)
        p <- length(g(samp[1, ], ...))
        g_samp <- matrix(nrow = n, ncol = p)
        for (i in seq_len(n)) g_samp[i, ] <- g(samp[i, ], ...)
        g_sigma <- cov(g_samp)
    } else {
        J <- jacobian(func = g, x = theta, ...)
        g_sigma <- J %*% vcov(x) %*% t(J)
    }

    mle(theta.hat = g(theta),
        sigma = g_sigma,
        info = ginv(g_sigma),
        nobs = nobs(x),
        superclasses = "mle_fit_rmap")
}
