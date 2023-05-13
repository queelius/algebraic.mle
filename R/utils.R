

entropy <- function(p) {
    -sum(p * log2(p))
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

clip_step <- function(step, max_norm = 1) {
    step_norm <- sqrt(sum(step^2))
    if (step_norm > max_norm) {
        step <- step * max_norm / step_norm
    }
    step
}

newton_raphson_l_bfgs_b <- function(f, gr, x0, eps = 1e-5,
                                    lower = -Inf, upper = Inf,
                                    max_iter = 1000, ...) {
    optim(
        par = x0, fn = f, gr = gr,
        method = "L-BFGS-B",
        lower = lower, upper = upper,
        control = list(maxit = max_iter, pgtol = eps),
        ...
    )
}

backtracking_line_search <- function(
    f,
    dir,
    x0,
    max_step,
    fix,
    sup,
    debug,
    max_iter,
    min_eta = 1e-8,
    r) {
    norm <- function(x) sqrt(sum(x^2))
    d.norm <- norm(dir)
    stopifnot(d.norm > 0)

    x.orig <- x0
    sol.max <- f(x0)
    sol.argmax <- x0
    eta <- max_step / d.norm
    iter <- 0L
    found_better <- FALSE

    while (eta >= min_eta && iter != max_iter) {
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
                if (norm(x.fixed - x.orig) > max_step) {
                    warning("fixed x yields a step size > max_step")
                }
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
                found_better <- TRUE
                break
            }
        }

        iter <- iter + 1L
        eta <- r * eta
    }

    list(argmax = sol.argmax, max = sol.max,
         found_better = found_better)
}


local_minimize_ls <- function(
    f, step, x0, sup = function(theta) TRUE,
    max_iter = 1000L, eta = 1, eps = 1e-5) {
    x <- x0
    f1 <- f(x)
    iter <- 0L
    convergence <- FALSE
    repeat {
        r <- eta
        d <- step(x0)
        f0 <- f(x0)
        repeat        {
            if (iter > max_iter) {
                break
            }
            iter <- iter + 1L
            x <- x0 - r * d
            if (sup(x)) {
                f1 <- f(theta)
                if (f1 <= f0) break
            }
            r <- r / 2
        }

        if (max(abs(x0 - x)) < eps) {
            convergence <- TRUE
            break
        }
        if (iter > max_iter) {
            break
        }

        x0 <- x
    }

    list(
        iter = iter, value = f1, par = x,
        convergence = convergence
    )
}

maximize_loglike_ls <- function(loglik, theta0, info, scr, eta = 1,
    clip = 1, max_iter = 10000L, eps = 1e-4, sup = function(theta) TRUE) {
    theta <- theta0
    l1 <- NULL
    iter <- 0L
    repeat {
        r <- eta
        I <- info(theta0)

        # we start off with a max step size of ||d|| = eta
        d <- clip_step(ginv(I) %*% scr(theta0), clip / eta)
        l0 <- loglik(theta0)
        repeat        {
            if (iter > max_iter) {
                break
            }
            iter <- iter + 1L
            theta <- theta0 + r * d
            if (sup(theta)) {
                l1 <- loglik(theta)
                if (l0 <= l1) {
                    break
                }
            }
            r <- r / 2
        }

        if (max(abs(theta0 - theta)) < eps) {
            convergence <- TRUE
            break
        }
        if (iter > max_iter) {
            break
        }

        theta0 <- theta
    }

    list(par = theta, value = l1, iter = iter, convergence = convergence)
}

grad_descent <- function(f, x0, df,
                         sup = function(theta) TRUE,
                         eps = 1e-10, lr = 1, debug = FALSE,
                         max_iter = 10000L) {
    iter <- 1L
    x1 <- NULL
    converged <- FALSE
    while (iter != max_iter) {
        f0 <- f(x0)
        g0 <- df(x0)
        if (debug) {
            cat("iter = ", iter, ", x0 = ", x0, ", f(x0) = ", f0, "\n")
        }
        eta <- lr
        good <- FALSE

        x1 <- x0
        while (iter != max_iter) {
            iter <- iter + 1L
            x <- x0 - eta * g0
            if (sup(x) && f(x) <= f0) {
                x1 <- x
                good <- TRUE
                break
            }
            eta <- eta / 2
        }

        if (!good) {
            break
        }

        if (max(abs(x1 - x0)) < eps) {
            converged <- TRUE
            break
        }
        x0 <- x1
    }
    return(list(param = x1, iter = iter, converged = converged))
}
