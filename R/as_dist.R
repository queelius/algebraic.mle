#' Convert an MLE to its asymptotic distribution.
#'
#' Internal helper that creates a \code{\link[algebraic.dist]{normal}} (univariate)
#' or \code{\link[algebraic.dist]{mvn}} (multivariate) distribution from an
#' \code{mle} object using the asymptotic normality of the MLE.
#'
#' @param x An \code{mle} object with a non-NULL variance-covariance matrix.
#' @return A \code{normal} or \code{mvn} distribution object.
#' @keywords internal
mle_to_dist <- function(x) {
    V <- vcov(x)
    if (is.null(V)) stop("Cannot convert to distribution: no variance-covariance matrix available.")
    mu <- params(x)
    if (nparams(x) == 1L) {
        normal(mu = as.numeric(mu), var = as.numeric(V))
    } else {
        if (!is.matrix(V)) V <- as.matrix(V)
        mvn(mu = as.numeric(mu), sigma = V)
    }
}

#' Convert an MLE to a distribution object.
#'
#' Converts an \code{mle} object to its asymptotic
#' \code{\link[algebraic.dist]{normal}} or \code{\link[algebraic.dist]{mvn}}
#' distribution. The MLE must have a variance-covariance matrix available.
#'
#' @param x An \code{mle} object.
#' @param ... Additional arguments (not used).
#' @return A \code{normal} (univariate) or \code{mvn} (multivariate) distribution.
#' @examples
#' fit <- mle(theta.hat = c(mu = 5), sigma = matrix(0.25))
#' d <- as_dist(fit)
#' mean(d)   # 5
#' vcov(d)   # 0.25
#'
#' fit2 <- mle(theta.hat = c(a = 1, b = 2), sigma = diag(c(0.1, 0.2)))
#' d2 <- as_dist(fit2)
#' mean(d2)  # c(1, 2)
#' @importFrom algebraic.dist as_dist normal mvn
#' @export
as_dist.mle <- function(x, ...) {
    mle_to_dist(x)
}

#' Convert a bootstrap MLE to an empirical distribution.
#'
#' Converts an \code{mle_boot} object to an
#' \code{\link[algebraic.dist]{empirical_dist}} built from the bootstrap
#' replicates.
#'
#' @param x An \code{mle_boot} object.
#' @param ... Additional arguments (not used).
#' @return An \code{empirical_dist} object.
#' @examples
#' set.seed(123)
#' x <- rexp(50, rate = 2)
#' rate_mle <- function(data, indices) 1 / mean(data[indices])
#' boot_result <- boot::boot(data = x, statistic = rate_mle, R = 200)
#' fit <- mle_boot(boot_result)
#' d <- as_dist(fit)
#' mean(d)
#' @importFrom algebraic.dist as_dist empirical_dist
#' @export
as_dist.mle_boot <- function(x, ...) {
    empirical_dist(x$t)
}
