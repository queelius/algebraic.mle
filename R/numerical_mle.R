#' sim_anneal
#' 
#' This function implements the simulated annealing algorithm,
#' which is a global optimization algorithm that is useful for
#' finding a good starting point for a local optimization algorithm.

#' We do not return this as an MLE object because, to be a good
#' estimate of the MLE, the gradient of `f` evaluated
#' at its solution should be close to zero, assuming the MLE
#' is interior to the domain of `f`. However, since this algorithm
#' is not guided by gradient information, it is not sensitive to
#' the gradient of `f` and instead only seeks to maximize `f`.
#' 
#' This also works for discrete optimization problems.
#' 
#' @param f Objective function to maximize, `f : R^d -> R`
#' @param x0 Initial guess
#' @param t_init Initial temperature
#' @param t_end Final temperature
#' @param sup Support function, returns TRUE if x is in the domain of f
#' @param neigh Neighborhood function, returns a random neighbor of x
#' @param max_iter Maximum number of iterations, used instead of t_end
#'                 if not NULL, defaults to NULL
#' @param iter_per_temp Number of iterations per temperature
#' @param ... Additional arguments to neigh
#' @param alpha Cooling factor
#' @param debug If TRUE, print debugging information to the console
#' @param trace If TRUE, track the history of positions and values
#' @return list with best solution (argmax) and its corresponding
#'         objective function value (max), and optionally trace_x and trace_fx
#' @export
sim_anneal <- function(
    f,
    x0,
    t_init = 100,
    t_end = 1e-3,
    alpha = 0.99,
    iter_per_temp = 100,
    sup = function(x) TRUE,
    neigh = function(x) x + rnorm(length(x)),
    max_iter = NULL,
    debug = FALSE,
    trace = FALSE, ...) {
    stopifnot(is.function(f), is.function(sup), is.function(neigh),
                t_init > 0, alpha > 0, alpha < 1, iter_per_temp > 0,
                is.logical(debug), is.logical(trace), sup(x0))

    if (!is.null(max_iter)) {
        stopifnot(max_iter > 0)
        if (!is.null(t_end)) {
            warning("t_end is ignored when max_iter is not NULL")
        }
        t_end <- t_init * alpha^max_iter
    }
    stopifnot(t_end > 0, t_end < t_init)

    argmax <- x0
    fmax <- f(argmax)
    fx0 <- fmax
    t <- t_init
    iter <- 0L
    trace_x <- matrix(nrow=0,ncol=length(x0))
    trace_fx <- c()

    while (t > t_end) {
        iter <- iter + 1L
        x <- neigh(x0,...)
        if (!sup(x)) {
            if (debug) {
                cat("x = (", x ,") not in support, skipping\n")
            }
            next
        }
        fx <- f(x)
        if (is.nan(fx)) {
            if (debug) {
                cat("f(x=",x,") is NaN, skipping\n")
            }
            next
        }

        if (fx > fmax) {
            if (debug) {
                cat("temp:", t, ", max:", fmax, ", argmax:", argmax, "\n")
            }

            argmax <- x
            fmax <- fx
        }

        if (exp((fx - fx0) / t) > runif(1)) {
            if (debug) {
                cat("temp:", t, ", val:", fx, ", x:", x, "\n")
            }
            x0 <- x
            fx0 <- fx
            trace_x <- rbind(trace_x, x)
            trace_fx <- c(trace_fx, fx)
        }
        if (iter %% iter_per_temp == 0) {
            t <- t * alpha
        }
    }

    list(argmax = argmax, max = fmax, trace_x = trace_x, trace_fx = trace_fx)
}

#' mle_optim
#'
#' This function takes the output of `optim` and turns it into an `mle` object.
#' @param res the output of `optim`
#' @importFrom MASS ginv
#' @export
mle_optim <- function(res) {
    sigma <- NULL
    info <- NULL
    if (!is.null(res$hessian))
    {
        info <- -res$hessian
        sigma <- ginv(info)
    }
    theta.hat <- mle(theta.hat=res$par,
        loglike=res$value,
        score=NULL,
        sigma=sigma,
        info=info,
        obs=NULL,
        nobs=NULL,
        superclasses=c("mle_optim","numerical_mle"))
    theta.hat$counts <- res$counts
    theta.hat$convergence <- res$convergence
    theta.hat$converged <- res$convergence == 0
    theta.hat$message <- res$message
    theta.hat$method <- res$method
    theta.hat
}

#' mle_local_search
#' 
#' This assumes the MLE is an interior point and that you have provided
#' an initial guess `theta0` that is near it. Use a global search method
#' like `sim_anneal` to find a good initial guess.
#' 
#' General local search method. It doesn't do anything fancy, but it is
#' more easily tweakable than `optim` in the `stats` package and can be
#' made to generate results when `optim` fails.
#' I use it for debugging and testing, but I also use it to implement
#' in a more transparent way the local search methods for newton-raphson
#' and gradient ascent.
#'
#' @param l function, log-likelihood function of type `R^p -> R`
#' @param theta0 numeric, initial guess, a column vector of size `p`
#' @param dir function, promising direction function of type `R^p -> R^p`
#' @param sup predicate function, domain of support function of type
#'            `R^p -> {TRUE,FALSE}` for the log-likelihood function `l`
#' @param eta numeric, learning rate, defaults to 1
#' @param max_iter integer, maximum number of iterations
#' @param max_iter_ls integer, maximum number of iterations for the line search
#' @param eps numeric, tolerance for convergence
#' @param r numeric, backtracking line search parameter, defaults to 0.8
#' @param proj function, projection function of type `R^p -> R^p`
#' @param tol function, tolerance function of type `R^p -> R`
#' @param debug logical, output debugging information if `TRUE`;
#'              defaults to `FALSE`
#' @param trace logical, if `TRUE` store the path of the search in the
#'             `trace` attribute of the output
#' @export
mle_local_search <- function(
    ll,
    theta0,
    dir,
    eps = 1e-5,
    proj = NULL,
    sup = NULL,
    tol = NULL,
    eta = 1,
    r = 0.8,
    max_iter = 1000L,
    max_iter_ls = 1000L,
    debug = FALSE,
    trace = FALSE) {
    stopifnot(eps > 0, eps < 1, eta > 0, r > 0, r < 1,
        max_iter > 0, max_iter_ls > 0,
        is.function(ll), is.function(dir))
    if (is.null(proj)) proj <- function(theta) theta
    stopifnot(is.function(proj))
    if (is.null(tol)) tol <- function(dtheta) abs(max(dtheta))
    stopifnot(is.function(tol))
    if (is.null(sup)) sup <- function(theta) TRUE
    stopifnot(is.function(sup))

    theta <- theta0
    max <- ll(theta0)
    trace_out <- list()

    for (iter in 1:max_iter)
    {
        d <- dir(theta0)
        if (debug) {
            cat("dir =", d, "theta0 =", theta0, ", loglike =", max, "\n")
        }

        # we use backtracking for an approximate line search
        res <- backtracking_line_search(
            f = ll,
            dir = d,
            x0 = theta0,
            max_step = eta,
            sup = sup,
            fix = proj,
            debug = debug,
            r = 0.8,
            max_iter = max_iter_ls,
        )

        theta <- res$argmax
        max <- res$max

        if (!res$found_better) {
            warning("backtracking line search failed to find a better solution")
            break
        }
        if (trace) {
            trace_out <- c(trace_out, list(theta))
        }
        if (tol(theta - theta0) < eps) {
            break
        }
        theta0 <- theta
    }

    theta.hat <- mle(
        theta.hat = theta,
        loglike = max,
        superclasses = c("mle_local_search","numerical_mle"))
    theta.hat$iter <- iter
    theta.hat$converged <- iter < max_iter
    theta.hat$learning_rate <- eta
    if (trace) {
        theta.hat$trace <- trace_out
    }

    theta.hat
}


#' MLE method using the Newton-Raphson method
#'
#' @param ll log-likelihood function of type `R^p -> R`
#' @param theta0 initial guess of theta with `p` components
#' @param info FIM function of type `R^p -> R^{p-by-q}`.
#' @param score score Function of type `R^p -> R^p`, the gradient of `ll``
#' @param inverted if `TRUE` `info` is covariance instead of FIM
#' @param eps numeric, stopping criterion
#' @param sup domain of support function of type `R^p -> {T,F}` for `ll`
#' @param eta numeric, learning rate, defaults to 1
#' @param r backing tracking parameter
#' @param max_iter maximum number of iterations
#' @param debug logical, output debugging info if `TRUE`; defaults to `FALSE`
#'
#' @importFrom MASS ginv
#' @export
mle_newton_raphson <- function(ll,
                               theta0,
                               info,
                               score,
                               eps = 1e-5,
                               sup = NULL,
                               eta = 1,
                               max_iter = 0L,
                               r = 0.5,
                               proj = NULL,
                               tol = NULL,
                               debug = FALSE,
                               inverted = FALSE,
                               trace = FALSE) {
    V <- NULL
    if (inverted) {
        V <- info
        info <- function(theta) ginv(V(theta))
    }
    else {
        V <- function(theta) ginv(info(theta))
    }
    dir <- function(theta) V(theta) %*% score(theta)
    res <- mle_local_search(
        ll = ll,
        theta0 = theta0,
        dir = dir,
        eps = eps,
        sup = sup,
        eta = eta,
        proj = proj,
        tol = tol,
        max_iter = max_iter,
        debug = debug,
        trace = trace)
    class(res) <- c("mle_newton_raphson", class(res))
    res$score <- score(res$theta.hat)
    res$sigma <- V(res$theta.hat)
    res$info <- info(res$theta)
    res
}

#' mle_gradient_ascent
#' 
#' MLE method using gradient ascent
#'
#' @param theta0 initial guess of theta with `p` components
#' @param ll log-likelihood function
#' @param score score function of type `R^p -> R^p`, the gradient of `ll`
#' @param eps numeric, stopping criterion
#' @param sup domain of support function of type `R^p -> {TRUE,FALSE}` for `ll`
#' @param eta numeric, learning rate, defaults to 1
#' @param r numeric, backing tracking parameter
#' @param max_iter integer, maximum iterations
#' @param tol tolerance function of type `R^p -> R`
#' @param proj projection function of type `R^p -> R^p`
#' @param debug logical, output debugging info if `TRUE`; defaults to `FALSE'
#' @importFrom MASS ginv
#' @importFrom numDeriv hessian
#' @export
mle_gradient_ascent <- function(ll,
                                theta0,
                                score,
                                eps = 1e-5,
                                sup = NULL,
                                eta = 1,
                                r = 0.5,
                                tol = NULL,
                                proj = NULL,
                                max_iter = 0L,
                                debug = FALSE,
                                trace = FALSE) {
    res <- mle_local_search(
        ll = ll,
        theta0 = theta0,
        dir = score,
        eps = eps,
        sup = sup,
        tol = tol,
        eta = eta,
        r = r,
        proj = proj,
        max_iter = max_iter,
        debug = debug,
        trace = trace)
    class(res) <- c("mle_gradient_ascent", class(res))
    res$info <- -hessian(ll, res$theta.hat)
    res$score <- score(res$theta.hat)
    res$sigma <- ginv(res$info)
    res
}

#' Function to determine whether an `mle` object has converged.
#'
#' @param x the `mle` object
#' @param ... additional arguments to pass
#' @export
is_converged <- function(x) {
    stopifnot(is(x, "numerical_mle"))
    theta.hat$converged
}

