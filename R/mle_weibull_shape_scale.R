#' MLE of the rate parameter when we assume the sample is i.i.d. and drawn
#' from the weibull distribution.
#'
#' @param x a sample of observations
#' @param k initial estimate of shape parameter k
#' @param eps we numerically solve the MLE equation, `|old-new| <= eps` is stopping condition
#' @param keep_obs Boolean, specifies whether to keep observations
#' @return an `mle` object.
#' @importFrom MASS ginv
#' @export
mle_weibull_shape_scale <- function(x, k0 = 1, eps = 1e-7, keep_obs = F) {
    n <- length(x)
    stopifnot(n > 0)
    stopifnot(k0 > 0)
    stopifnot(eps > 0)

    s <- mean(log(x))
    repeat    {
        k1 <- k0
        k0 <- 1 / (sum(x^k1 * log(x)) / sum(x^k1) - s)
        if (abs(k0 - k1) < eps) {
            break
        }
    }
    theta.hat <- c(k0, mean(x^k0)^(1 / k0))
    l <- weibull_shape_scale_loglike(x)(theta.hat)
    names(theta.hat) <- c("shape", "scale")
    fim <- weibull_shape_scale_fim(x)(theta.hat)
    rownames(fim) <- names(theta.hat)
    colnames(fim) <- names(theta.hat)
    sigma <- ginv(fim)
    rownames(sigma) <- names(theta.hat)
    colnames(sigma) <- names(theta.hat)

    mle(
        theta.hat = theta.hat,
        loglike = l,
        score = weibull_shape_scale_score(x)(theta.hat),
        sigma = sigma,
        info = fim,
        obs = if (keep_obs) x else NULL,
        nobs = n,
        superclasses = c("mle_weibull_shape_scale")
    )
}

#' log-likelihood function generator given data `x` for the weibull
#' distribution.
#'
#' The returned log-likelihood function takes a single vector `theta` of
#' size `2` (at least) where the first component is the shape parameter
#' `k` and the second component is the scale parameter `lambda`.
#'
#' It can be used in statistical models or optimization algorithms to estimate
#' the parameters of the Weibull distribution.
#'
#' @param x data
#' @export
weibull_shape_scale_loglike <- function(x) {
    n <- length(x)
    stopifnot(n > 0) # we need at least one observation

    stopifnot(all(x >= 0)) # if any values in `x` negative, the two
    # parameter weibull distribution is not a good fit to
    # the data. we stop if so to avoid taking its log,
    # which is undefined for negative values.
    function(theta) {
        n * log(theta[1] / theta[2]) +
            (theta[1] - 1) * sum(log(x / theta[2])) -
            sum((x / theta[2])^theta[1])
    }
}

#' score function generator given data `x` for the weibull
#' distribution given a simple random sample.
#'
#' @param x data
#' @export
weibull_shape_scale_score <- function(x) {
    n <- length(x)
    function(theta) {
        c(
            n / theta[1] + sum(log(x / theta[2])) - sum((x / theta[2])^theta[1] * log(x / theta[2])),
            theta[1] / theta[2] * (sum((x / theta[2])^theta[1]) - n)
        )
    }
}

#' log-likelihood function generator given data `x` for the weibull
#' distribution
#'
#' @param x data
#' @export
weibull_shape_scale_fim <- function(x) {
    n <- length(x)
    function(theta) {
        k <- theta[1]
        a <- theta[2]
        d <- sum((x / a)^k * log(x / a)^2)
        matrix(c(
            n / k^2 + d,
            n / a - d / a,
            n / a - d / a,
            n * k / a^2 + k * (k + 1) / a^2 * sum((x / a)^k)
        ), nrow = 2)
    }
}
