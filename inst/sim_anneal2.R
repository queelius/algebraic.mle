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
#' @param par Initial guess
#' @param fn Objective function to maximize
#' @param control List of optional arguments
#' @param ... Additional arguments that may be passed; loads into options,
#'            and is also passed into `neigh`
#' @describeIn sim_anneal control
#' @field t_init Initial temperature, defaults to 100
#' @field t_end Final temperature, defaults to 0
#' @field alpha Numberic, cooling factor, defaults to .95
#' @field REPORT The frequency of reports if control$debug > 0.
#'        Defaults to every 100 iterations.
#' @field it_per_temp Ineger, iterations per temperature, defaults to 100
#' @field maxit Integer, maximum number of iterations, defaults to 100000
#' @field accept_p Acceptance probability function, defaults to
#'        `runif(1) < exp((val0 - val1) / temp)`, where `val0` is the
#'        function value at the current position and `val1` is the
#'        function value at the proposed position. We pass `temp`
#'        (temperature), `val1` (new value at candidate position `par1`),
#'        `val0` (old value at current position), `it` (iteration), and
#'        `...` as arguments to `accept_p`.
#' @field fnscale Scaling factor for `fn`, defaults to 1. If negative,
#'        then turns the problem into a maximization problem.
#' @field sup Support function, returns TRUE if `par` is in the domain of `fn`
#' @field proj Projection function, returns a vector in the domain of `fn`
#' @field neigh Neighborhood function, returns a random neighbor of `par`,
#'              defaults to `par + rnorm(length(par))`. We pass `par`
#'              `temp`, `it`, `value`, and `...` as arguments to `neigh`.
#' @field trace logical, whether to store current changes in position
#'              and associated other values in a `trace_info` matrix.
#' @return list, members include:
#'            `par`: best solution
#'            `value`: function value `fn(par)`
#'            `fn_count`: a count of `fn` invocations
#'            `accepted`: a count of accepted moves
#'            `trace_info`: a matrix of trace information (optional)
#' @importFrom stats runif
#' @export
sim_anneal <- function(par, fn, control = list(),
    logger_config = default_log_config, ...)
{
    m <- length(par)
    defaults <- list(
        t_init = 100,
        t_end = 1e-3,
        fnscale = 1,
        alpha = 0.95,
        it_per_temp = 100L,
        maxit = 10000L,
        accept_p = function(temp, val1, val0, ...) {
            val1 < val0 || (runif(1) < exp((val0 - val1) / temp))
        },
        trace = FALSE,
        debug = 0L,
        REPORT = 100L,
        trace_info_size_inc = 10000L,
        sup = function(par, ...) TRUE,
        proj = function(par, ...) par,
        neigh = function(par, ...) par + rnorm(m))

    control <- modifyList(defaults, control)
    control <- modifyList(control, list(...))
    stopifnot(is.numeric(control$REPORT), is.numeric(control$debug),
              is.numeric(control$it_per_temp), is.numeric(control$maxit),
              is.numeric(control$t_init), is.numeric(control$t_end),
              is.numeric(control$alpha), is.function(control$neigh),
              is.function(control$sup), is.function(control$accept_p),
              is.function(control$proj), is.logical(control$trace),
              control$t_init > 0, control$t_end > 0,
              control$t_end <= control$t_init,
              control$alpha > 0, control$alpha < 1,
              control$trace_info_size_inc > 0,
              control$it_per_temp > 0, control$REPORT > 0,
              control$maxit > 0, control$maxit != Inf)

    if (!is.function(fn)) {
        stop("`fn` must be a function (to maximize)")
    }
    if (!control$sup(par)) {
        stop("initial guess `par` not in support")
    }

    par0 <- par
    val0 <- fn(par0) / control$fnscale
    best_par <- par0
    best_val <- val0
    accepted <- 0L
    t <- control$t_init
    it <- 0L

    trace_info <- NULL
    if (control$trace) {
        cnames <- c(paste0("par", 1:m), "value", "it", "temp", "best")
        trace_info <- matrix(NA,
            nrow = control$trace_info_size_inc,
            ncol = length(cnames))
        colnames(trace_info) <- cnames
        trace_app <- trace_info
        append_trace <- function(data) {
            if (accepted > nrow(trace_info)) {
                trace_info <<- rbind(trace_info, trace_app)
            }
            trace_info[accepted, ] <<- data
        }
    }

    debug_out <- function() {
        par_str <- paste(sprintf("%.5f", par0), collapse=", ")
        message(paste(
            "Iteration", sprintf("%4d", it), "|",
            "Temp:", sprintf("%.4f", t), "|",
            "Solution: (", par_str, ")", "|",
            "Value:", sprintf("%.4f", control$fnscale) * val0, "|",
            "Best: ", sprintf("%1d", best_val == val0), sep = " "))
    }

    repeat {

        if (it > control$maxit) {
            break
        }
        it <- it + 1L

        if (it %% control$it_per_temp == 0) {
            t <- t * control$alpha
        }

        if (!is.null(control$t_end) && t < control$t_end) {
            break
        }

        par1 <- control$proj(control$neigh(
            par = par0,
            temp = t,
            value = val0,
            it = it,
            ...))

        if (!control$sup(par1)) {
            if (control$debug > 2L) {
                # pretty print the vector `par1`
                message("projected neighbor not in the support: ",
                        paste(sprintf("%.5f", par1), collapse = ", "))
            }
            next
        }

        val1 <- fn(par1) / control$fnscale
        if (is.nan(val1)) {
            if (control$debug > 2L) {
                message("`fn` was not-a-number at ", par1)
            }
            next
        }

        if (val1 < best_val) {

            if (control$debug > 0L) {
                par_str <- paste(sprintf("%.5f", par1), collapse = ", ")
                message("Iteration ", it, " | ",
                        "Found better ", val1, " at (", par_str, ")")
            }
            best_par <- par1
            best_val <- val1
        }

        if (control$accept_p(t, val1, val0, it, ...)) {
            par0 <- par1
            val0 <- val1
            accepted <- accepted + 1L

            if (control$trace) {
                append_trace(c(par0, control$fnscale * val0, it,
                    t, best_val == val0))
            }
        }

        if (control$debug > 1L && it %% control$REPORT == 0) {
            debug_out()
        }
    }

    res <- list(par = best_par,
                value = control$fnscale * best_val,
                fn_count = it,
                accepted = accepted)

    if (control$trace) {
        res$trace_info <- trace_info[1:accepted, ]
    }
    res
}
