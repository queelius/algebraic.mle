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
#' @param options List of options
#' * `t_init` Initial temperature
#' * `t_end` Final temperature
#' * `alpha` Cooling factor
#' * `iter_per_temp` Number of iterations per temperature
#' * `max_iter` Maximum number of iterations, used instead of t_end
#'             if not NULL, defaults to NULL
#' * `debug` If TRUE, print debugging information to the console
#' * `trace` If TRUE, track the history of positions and values
#' * `sup` Support function, returns TRUE if x is in the domain of f
#' * `neigh` Neighborhood function, returns a random neighbor of x
#' * `...` Additional arguments to neigh
#' @return list with best solution (argmax) and its corresponding
#'         objective function value (max), and optionally path
#' @importFrom stats runif
#' @export
sim_anneal <- function(f, x0, options=list(), ...) {
    defaults <- list(
        t_init = 100,
        t_end = 1e-6,
        alpha = 0.95,
        iter_per_temp = 100,
        max_iter = 1e20,
        debug = FALSE,
        trace = FALSE,
        sup = function(x) TRUE,
        neigh = function(x) x + rnorm(length(x))
    )
    options <- modifyList(defaults, options)
    stopifnot(options$t_init > 0,
              options$alpha > 0, options$alpha < 1,
              options$iter_per_temp > 0,
              (is.null(options$max_iter) || options$max_iter > 0),
              is.function(f),
              is.function(options$neigh),
              is.function(options$sup),
              is.logical(options$debug),
              is.logical(options$trace))

    if (!options$sup(x0)) {
        stop("initial guess x0 not in support")
    }

    argmax <- x0
    fmax <- f(argmax)
    fx0 <- fmax
    t <- options$t_init
    iter <- 0L
    path <- matrix(nrow=0,ncol=length(x0))

    while (t > options$t_end) {
        iter <- iter + 1L
        x <- options$neigh(x0,...)
        if (!options$sup(x)) {
            if (options$debug) {
                cat("x = (", x ,") not in support, skipping\n")
            }
            next
        }
        fx <- f(x)
        if (is.nan(fx)) {
            if (options$debug) {
                cat("f(x=",x,") is NaN, skipping\n")
            }
            next
        }

        if (fx > fmax) {
            if (options$debug) {
                cat("temp:", t, ", max:", fmax, ", argmax:", argmax, "\n")
            }

            argmax <- x
            fmax <- fx
        }

        if (exp((fx - fx0) / t) > runif(1)) {
            if (options$debug) {
                cat("temp:", t, ", val:", fx, ", x:", x, "\n")
            }
            x0 <- x
            fx0 <- fx
            if (options$trace) {
                path <- rbind(path, x)
            }
        }
        if (iter %% options$iter_per_temp == 0) {
            t <- t * options$alpha
        }
    }

    sol <- list(argmax = argmax, max = fmax, options = options)
    if (options$trace) {
        sol$path <- path
    }
    sol
}

#' mle_optim
#'
#' This function takes the output of `optim` and turns it into an `mle` object.
#' 
#' @param sol the output of `optim`
#' @return a `numerical_mle` object, specialized for `optim` (stats package)
#' solutions.
#' @importFrom MASS ginv
#' @export
mle_optim <- function(sol) {
    sigma <- NULL
    info <- NULL
    if (!is.null(sol$hessian))
    {
        info <- -sol$hessian
        sigma <- ginv(info)
    }
    mle.sol <- mle(
        theta.hat=sol$par,
        loglike=sol$value,
        score=NULL,
        sigma=sigma,
        info=info,
        obs=NULL,
        nobs=NULL,
        superclasses=c("mle_optim","numerical_mle"))
    mle.sol$converged <- sol$convergence == 0
    mle.sol$optim_data <- sol
    mle.sol
}

#' mle_local_search
#' 
#' This assumes the MLE is an interior point and that you have provided
#' an initial guess `theta0` that is near it. Use a global search method
#' like `sim_anneal` to find a good initial guess.
#' 
#' @param ll function, log-likelihood function
#' @param theta0 numeric, initial guess
#' @param dir function, promising direction function
#' @param options list, options for the local search:
#' * `sup` predicate function, domain of support for log-likelihood
#' * `eta` numeric, learning rate, defaults to 1
#' * `max_iter` integer, maximum number of iterations, defaults to 1000
#' * `max_iter_ls` integer, maximum number of iterations for the line search,
#' defaults to 1000
#' * `abs_tol` numeric, tolerance for convergence, defaults to NULL
#' (use `rel_tol`)
#' * `rel_tol` numeric, relative tolerance for convergence, defaults to 1e-5
#' * `r` numeric, backtracking line search parameter, defaults to 0.8
#' * `proj` function, projection function to enforce domain of support
#' * `norm` function, we pass the difference of successive theta updates
#' to the norm to check for convergence, defaults to the infinity norm.
#' * `debug` logical, output debugging information if TRUE; default FALSE
#' * `trace` logical, if TRUE store the path of the search in the `path`
#' attribute of the output
#' @return an `mle` object with the following additional attributes:
#' * `iter` integer, the number of iterations
#' * `converged` logical, whether the algorithm converged
#' * `options` list, the options used for the search
#' * `path` matrix, the path of the search if `trace` is TRUE
#' @export
mle_local_search <- function(
    ll,
    dir,
    theta0,
    options = list())
{
    defaults <- list(
        abs_tol = NULL,
        rel_tol = 1e-5,
        proj = function(theta) theta,
        sup = function(theta) TRUE,
        norm = function(dtheta) max(abs(dtheta)),
        eta = 1,
        r = 0.8,
        max_iter = 1000L,
        max_iter_ls = 1000L,
        debug = FALSE,
        trace = FALSE
    )

    options <- modifyList(defaults, options)
    stopifnot(options$eta > 0, options$r > 0, options$r < 1,
              options$max_iter > 0, options$max_iter_ls > 0,
              is.logical(options$debug),
              is.logical(options$trace),
              is.function(ll), is.function(dir),
              is.function(options$proj),
              is.function(options$norm),
              is.function(options$sup))

    if (is.null(options$rel_tol)) {
        stopifnot(!is.null(options$abs_tol), options$abs_tol > 0)
    }
    else {
        stopifnot(!is.null(options$rel_tol), options$rel_tol > 0)
    }

    theta <- theta0
    max <- ll(theta0)
    converged <- FALSE
    path <- NULL
    if (options$trace) {
        path <- matrix(nrow=options$max_iter, ncol=length(theta0))
    }
    iter <- 1L
    repeat
    {
        d <- dir(theta0)
        if (options$debug) {
            cat("dir =", d, "theta0 =", theta0, ", loglike =", max, "\n")
        }

        # we use backtracking for an approximate line search
        res <- backtracking_line_search(
            f = ll,
            dir = d,
            x0 = theta0,
            max_step = options$eta,
            sup = options$sup,
            fix = options$proj,
            debug = options$debug,
            r = options$r,
            max_iter = options$max_iter_ls)

        theta <- res$argmax
        max <- res$max

        if (!res$found_better) {
            if (options$debug) {
                cat("backtracking line search failed to find a better solution\n")
            }
            warning("backtracking line search failed to find a better solution")
            break
        }
        if (options$trace) {
            path[iter,] <- as.vector(theta)
        }

        if (is.null(options$abs_tol)) {
            if (options$norm(theta - theta0) < options$rel_tol * norm(theta)) {
                converged <- TRUE
                break
            }            
        }
        else {
            if (options$norm(theta - theta0) < options$abs_tol) {
                converged <- TRUE
                break
            }
        }
        theta0 <- theta
        iter <- iter + 1L
        if (iter > options$max_iter) {
            break
        }
    }

    if (iter == options$max_iter) {
        if (options$debug) {
            cat("maximum number of iterations reached\n")
        }
        warning("maximum number of iterations reached")
    }

    sol <- mle(theta.hat = theta, loglike = max,
        superclasses = c("mle_local_search", "numerical_mle"))
    sol$iter <- iter
    sol$converged <- converged
    sol$options <- options

    if (options$trace) {
        sol$path <- path[1:iter, , drop = FALSE]
    }
    sol
}

#' MLE method using the Newton-Raphson method
#'
#' @param ll log-likelihood function
#' @param theta0 initial guess
#' @param info FIM function
#' @param score score Function, the gradient of log-likelihood
#' @param inverted logical, if TRUE `info` is covariance instead of FIM
#' @param options list, options for the local search, see `mle_local_search`
#' @return an object of class `mle_newton_raphson`, which is an `mle` object
#' @importFrom MASS ginv
#' @export
mle_newton_raphson <- function(
    ll,
    score,
    info,
    theta0,
    inverted = FALSE,
    options = list())
{
    stopifnot(is.function(ll),
              is.function(info),
              is.function(score))

    V <- NULL
    if (inverted) {
        V <- info
        info <- function(theta) ginv(V(theta))
    }
    else {
        V <- function(theta) ginv(info(theta))
    }
    dir <- function(theta) V(theta) %*% score(theta)
    sol <- mle_local_search(
        ll = ll,
        theta0 = theta0,
        dir = dir,
        options = options)

    class(sol) <- c("mle_newton_raphson", class(sol))
    sol$score <- score(sol$theta.hat)
    sol$sigma <- V(sol$theta.hat)
    sol$info <- info(sol$theta)
    sol
}

#' mle_gradient_ascent
#' 
#' MLE method using gradient ascent
#'
#' @param theta0 initial guess of theta with `p` components
#' @param ll log-likelihood function
#' @param score score function of type `R^p -> R^p`, the gradient of `ll`
#' @importFrom MASS ginv
#' @importFrom numDeriv hessian
#' @export
mle_gradient_ascent <- function(ll, theta0, score, options) {
    sol <- mle_local_search(
        ll = ll,
        theta0 = theta0,
        dir = score,
        options = options)
    class(sol) <- c("mle_gradient_ascent", class(sol))
    sol$info <- -hessian(ll, sol$theta.hat)
    sol$score <- score(sol$theta.hat)
    sol$sigma <- ginv(sol$info)
    sol
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

