#' MLE of the (mu,var) parameter vector when we assume the sample is i.i.d. and
#' drawn from the normal distribution.
#'
#' Of course, the draws are unlikely to be normal, but the normal distribution
#' is often a good model. Hypothesis testing, such as relative
#' likelihoods, can be used to assess the appropriateness of the normal
#' model to the data.
#'
#' @param x a sample of observations
#' @param keep_obs whether to store observations in `mle_normal` object
#' @return an `mle` object.
#' @export
mle_normal <- function(x, keep_obs = TRUE) {
    n <- length(x)
    stopifnot(n > 0) # we need at least one observation

    theta.hat <- c(mean(x), (1 - 1 / n) * var(x))
    names(theta.hat) <- c("mu", "var")
    info.hat <- normal_fim(x)(theta.hat)

    mle(
        theta.hat = theta.hat,
        loglike = normal_loglike(x)(theta.hat),
        score = normal_score(x)(theta.hat),
        sigma = ginv(info.hat),
        info = info.hat,
        obs = if (keep_obs) x else NULL,
        nobs = n,
        superclasses = c("mle_normal")
    )
}

#' log-likelihood function generator given data `x` for the normal
#' distribution
#'
#' @param x data
#' @export
normal_loglike <- function(x) {
    # more efficiently computes:
    #    sum(dnorm(x,mean=theta[1],sd=sqrt(theta[2]),log=T))
    n <- length(x)
    S <- sum(x)
    S2 <- sum(x^2)
    function(theta) {
        -n / 2 * log(2 * pi * theta[2]) -
            1 / (2 * theta[2]) * (S2 - 2 * theta[1] * S + n * theta[1]^2)
    }
}

#' Score function generator given data `x` for the normal
#' distribution
#'
#' @param x data (simple random sample)
#' @export
normal_score <- function(x) {
    n <- length(x)
    S <- sum(x)
    S2 <- sum(x^2)
    function(theta) {
        matrix(
            c(0.5 / theta[2] * (S - n * theta[1]),
              -0.5 * n / theta[2] + 0.5 / theta[2]^2 * (S2 - 2 * theta[1] *
                    S + n * theta[1]^2)),
            nrow = 2
        )
    }
}

#' Fisher information matrix generator given data `x` for the normal
#' distribution
#'
#' @param x data (simple random sample)
#' @param observed whether to return the observed fisher information, default is
#'                 `TRUE`
#' @export
normal_fim <- function(x, observed = TRUE) {
    n <- length(x)
    if (observed) {
        S <- sum(x)
        S2 <- sum(x^2)
        return(function(theta) {
            covar <- 0.5 / theta[2]^2 * (S - n * theta[1])
            matrix(
                c(
                    n / theta[2],
                    covar,
                    covar,
                    -0.5 * n / theta[2]^2 + 1 / theta[2]^3 * (S2 - 2 * theta[1] * S + n * theta[1]^2)
                ),
                nrow = 2
            )
        })
    } else {
        return(function(theta) {
            matrix(
                c(n / theta[2],    0,
                0,                 0.5 * n / theta[2]^2),
                nrow = 2
            )
        })
    }
}

#' Computes the bias of an `mle_normal` object.
#' 
#' This function computes the bias of an `mle_normal` object.
#' It is a method for the generic function `bias`.
#' 
#' @param x An `mle_normal` object (subclass of `mle`) to compute the bias of.
#' @param par The true parameter value. If this is unknown (NULL), the bias is
#'            estimated.
#' @param ... Additional arguments (currently unused).
#' @return A numeric vector of length 2, the bias of the `mle_normal` estimator.
#' @seealso \code{\link{bias}} for the generic function.
#' @export
bias.mle_normal <- function(x, par = NULL, ...) {
    if (is.null(par)) {
        c(0, -1 / nobs(x) * point(x)[2])
    } else {
        c(0, -1 / nobs(x) * par[2])
    }
}