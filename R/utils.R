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


welford_stats <- function() {
    n <- 0
    mean <- 0
    M2 <- 0

    list(
        update = function(x) {
            n <<- n + 1
            delta <- x - mean
            mean <<- mean + delta / n
            delta2 <- x - mean
            M2 <<- M2 + delta * delta2
        },
        get_mean = function() {
            mean
        },
        get_variance = function() {
            if (n < 2) NA else M2 / (n - 1)
        }
    )
}

confidence_intervals <- function(V, theta, parm = NULL, level = 0.95, ...) {
    sigma <- diag(V)
    p <- length(theta)
    q <- stats::qnorm(level)
    if (is.null(parm)) {
        parm <- 1:p
    }

    parm <- parm[parm >= 1 & parm <= p]
    ci <- matrix(nrow = length(parm), ncol = 2)
    colnames(ci) <- c(
        paste((1 - level) / 2 * 100, "%"),
        paste((1 - (1 - level) / 2) * 100, "%")
    )

    i <- 1
    for (j in parm)
    {
        ci[i, ] <- c(
            theta[j] - q * sqrt(sigma[j]),
            theta[j] + q * sqrt(sigma[j])
        )
        i <- i + 1
    }
    rownames(ci) <- parm
    ci
}
