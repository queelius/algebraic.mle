#' Label a confidence interval matrix with column and row names.
#'
#' @param ci A two-column matrix of lower and upper bounds.
#' @param theta Named numeric vector of parameter estimates.
#' @param alpha The lower-tail probability (e.g., 0.025 for a 95 percent CI).
#' @return The \code{ci} matrix with column names set to percentile labels
#'   and row names set to parameter names (or \code{"param1"}, \code{"param2"}, ...).
#' @keywords internal
label_ci <- function(ci, theta, alpha) {
    colnames(ci) <- c(paste0(alpha * 100, "%"),
                      paste0((1 - alpha) * 100, "%"))
    rownames(ci) <- if (is.null(names(theta))) {
        paste0("param", seq_len(length(theta)))
    } else {
        names(theta)
    }
    ci
}

#' Function to compute the confidence intervals from a variance-covariance matrix
#'
#' @param sigma either the variance-covariance matrix or the vector of variances
#'              of the parameter estimator
#' @param theta the point estimate
#' @param level confidence level, defaults to 0.95 (alpha=.05)
#' @return Matrix of confidence intervals with rows for each parameter and
#'   columns for lower and upper bounds.
#' @examples
#' # Compute CI for a bivariate parameter
#' theta <- c(mu = 5.2, sigma2 = 4.1)
#' vcov_matrix <- diag(c(0.1, 0.5))  # Variance of estimators
#'
#' confint_from_sigma(vcov_matrix, theta)
#' confint_from_sigma(vcov_matrix, theta, level = 0.99)
#' @importFrom stats qnorm
#' @export
confint_from_sigma <- function(sigma, theta, level = .95) {
    stopifnot(is.numeric(level), level >= 0, level <= 1)
    if (is.matrix(sigma)) sigma <- diag(sigma)
    alpha <- (1 - level) / 2
    q <- stats::qnorm(1 - alpha)
    se <- sqrt(sigma)
    ci <- cbind(theta - q * se, theta + q * se)
    label_ci(ci, theta, alpha)
}
