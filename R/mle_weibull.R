#' MLE of the rate parameter when we assume the sample is i.i.d. and drawn
#' from the weibull distribution.
#'
#' @param x a sample of observations
#' @param k initial estimate of shape parameter k
#' @param rel_tol when numerically solving the MLE equation,
#'                `|new-old|/|new| <= rel_tol` is stopping condition
#' @param keep_obs Boolean, specifies whether to keep observations
#' @return an `mle` object.
#' @importFrom MASS ginv
#' @export
mle_weibull <- function(x, k0 = 1, rel_tol = 1e-7, keep_obs = FALSE) {
    n <- length(x)
    stopifnot(n > 0)
    stopifnot(k0 > 0)
    stopifnot(rel_tol > 0)

    log_x <- log(x)
    s <- mean(log_x)
    repeat    {
        k1 <- k0
        k0 <- 1 / (sum(x^k1 * log_x) / sum(x^k1) - s)
        if (abs(k1 - k0) / abs(k1) < rel_tol) {
            break
        }
    }
    par.hat <- c(k0, mean(x^k0)^(1 / k0))
    ll <- weibull_loglike(x)(par.hat)
    fim <- weibull_fim(x)(par.hat)
    sigma <- ginv(fim)
    cnames <- c("shape", "scale")
    names(par.hat) <- cnames
    rownames(fim) <- cnames
    colnames(fim) <- cnames
    rownames(sigma) <- cnames
    colnames(sigma) <- cnames

    mle(theta.hat = par.hat,
        loglike = ll,
        score = weibull_score(x)(par.hat),
        sigma = sigma,
        info = fim,
        obs = if (keep_obs) x else NULL,
        nobs = n,
        superclasses = c("mle_weibull")
    )
}

#' A log-likelihood function generator given data `x` for the weibull
#' distribution.
#'
#' The returned log-likelihood function takes a single vector `par` of
#' size `2` (at least) where the first component is the shape parameter
#' `k` and the second component is the scale parameter `lambda` such
#' that E[X] = lambda * GAMMA(1 + 1/k).
#'
#' @param x data
#' @return log-likelihood function for the weibull distribution
#'         with shape parameter `k` and scale parameter `lambda`
#'         given data `x`. (Accepts a parameter vector `par` of size 2.)
#' @export
weibull_loglike <- function(x) {
    n <- length(x)
    stopifnot(n > 0, all(x > 0)) # we need at least one observation
                                # and all observations must be positive
    sum_log_x <- sum(log(x))
    function (par) {
        k <- par[1] # shape
        l <- par[2] # scale
        n * log(k / l^k) +
            (k - 1) * sum_log_x - sum((x / l)^k)
    }
}

#' score function generator given data `x` for the weibull
#' distributison given a simple random sample.
#'
#' @param x data
#' @export
weibull_score <- function(x) {
    n <- length(x)
    stopifnot(n > 0, all(x > 0)) # we need at least one observation
                                 # and all observations must be positive
    s <- sum(log(x))
    function(par) {
        k <- par[1] # shape
        l <- par[2] # scale
        x_l <- x / l
        x_l_k <- x_l^k
        c(s + n / k - n * log(l) - sum(x_l_k * log(x_l)),
          k * sum(x_l_k) / l - n * k / l)
    }
}

#' weibull_fim
#' 
#' Observed Fisher information generator given data `x` for the weibull
#' distribution
#'
#' @param x observed data
#' @export
weibull_fim <- function(x) {
    n <- length(x)
    function(par) {
        k <- par[1]
        l <- par[2]
        x_l_k <- (x / l)^k
        log_x_l <- log(x / l)
        p_k_l <- n / l - sum(x_l_k * (1 + k * log_x_l)) / l
        matrix(c(
            # d^2/dk^2
            n / k^2 + sum(x_l_k * log_x_l^2),
            # d^2/dkdl
            p_k_l,
            # d^2/dkdl
            p_k_l,
            # d^2/dl^2
            -n * k / l^2 + k * (k + 1) / l^2 * sum(x_l_k)
        ), nrow = 2)
    }
}
