#' @export
entropy <- function(p) {
    -sum(p * log2(p))
}

#' @export
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

#' @export
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

#' @export
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

#' @export
clip_step <- function(step, max_norm=1) {
    step_norm <- sqrt(sum(step^2))
    if (step_norm > max_norm)
        step <- step * max_norm / step_norm
    step
}

#' @export
newton_raphson_l_bfgs_b <- function(
        f, gr, x0, eps = 1e-5,
        lower = -Inf, upper = Inf,
        max_iter = 1000, ...)
{
    optim(
        par = x0, fn = f, gr = gr,
        method = "L-BFGS-B",
        lower = lower, upper = upper,
        control = list(maxit = max_iter, pgtol = eps),
        ...)
}

#' @export
backtracking_line_search <- function(
    f,
    dir,
    x0,
    max_step = 1,
    fix = NULL,
    sup = function(x) TRUE,
    debug = FALSE,
    r = 0.5)
{
    norm <- function(x) sqrt(sum(x^2))
    d.norm <- norm(dir)
    stopifnot(d.norm > 0)

    x.orig <- x0
    sol.max <- f(x0)
    sol.argmax <- x0
    eta <- max_step / d.norm

    while (eta > 1e-8) {
        x <- x0 + eta * dir
        if (debug) {
            cat("x = ", x, ", eta = ", eta, "\n")
        }

        if (!sup(x)) {
            if (debug) {
                cat("x = ", x, " not in domain of f\n")
            }

            if (!is.null(fix)) {
                if (debug) {
                    cat("fixing x to ", fix(x), "\n")
                }

                # adjust eta for the fixed x
                x.fixed <- fix(x)
                eta <- sqrt(sum((x.fixed - x.orig)^2))
                stopifnot(norm(x.fixed-x.orig) <= max_step)
                x <- x.fixed
            }
        }

        if (sup(x)) {
            fx <- f(x)
            if (debug && is.nan(fx)) {
                cat("x = ", x, " is NaN\n")
            }
            if (!is.nan(fx) && fx > sol.max) {
                sol.argmax <- x
                sol.max <- fx
                break
            }
        }

        eta <- r * eta
    }

    list(argmax = sol.argmax, max = sol.max, found_better = eta > 1e-10)
}

