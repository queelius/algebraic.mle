#' Function to compute the confidence intervals from a variance-covariance matrix
#'
#' @param sigma either the variance-covariance matrix or the vector of variances
#'              of the parameter estimator
#' @param theta the point estimate
#' @param level confidence level, defaults to 0.95 (alpha=.05)
#' @importFrom stats qnorm
#' @export
confint_from_sigma <- function(sigma, theta, level = .95) {
    stopifnot(is.numeric(level), level >= 0, level <= 1)
    if (is.matrix(sigma)) {
        sigma <- diag(sigma)
    }
    alpha <- (1 - level) / 2
    p <- length(theta)
    q <- stats::qnorm(1 - alpha)
    ci <- matrix(nrow = p, ncol = 2)
    colnames(ci) <- c(paste0(alpha * 100, "%"),
                      paste0((1 - alpha) * 100, "%"))

    for (j in 1:p) {
        ci[j, ] <- c(
            theta[j] - q * sqrt(sigma[j]),
            theta[j] + q * sqrt(sigma[j])
        )
    }

    if (is.null(names(theta))) {
        rownames(ci) <- paste0("param", 1:p)
    } else {
        rownames(ci) <- names(theta)[1:p]
    }
    ci
}


get_first_attr <- function(xs, g, props) {
    for (x in xs)
    {
        y <- g(x)
        for (prop in props)
        {
            if (!is.null(prop(y))) {
                return(prop(y))
            }
        }
    }
    NULL
}
