#' newton_raphson
#'
#' Performs Newton-Raphson to find a solution that finds a local
#' maxima of the objective function `fn`. We assume the maxima
#' is an interior point of the support and that an initial guess par
#' is near the local maxima.
#' 
#' Use a global search method like `sim_anneal` to find a good initial
#' guess that is near a global maximum of `fn`.
#' 
#' @param fn function, objective function to maximized
#' @param par numeric, initial guess
#' @param gr function, gradient function
#' @param hess function, hessian function
#' @param control list, control for the local search, see function description.
#' @describeIn newton_raphson control
#' @field eta numeric, learning rate
#' @field fnscale numeric, scaling factor for `fn`. If negative, then
#'       turns the problem into a maximization problem.
#' @field inverted logical, whether the hessian function is already inverted
#' @field maxit integer, maximum number of iterations, defaults to 1000
#' @field convergence function, check for convergence, defaults to
#'        `||gr(par)||` and `||(par1 - par0)||` being approximately zero
#'        (both tests use absolute tolerance,  see `abs_tol`).
#' @field proj function, projection function to enforce domain of support
#' @field trace logical, keep track of trace information
#' @field abs_tol numeric, absolute tolerance for convergence, which is used
#'        by the default `convergence` function
#' @field REPORT integer, frequency of tracer reports, defautls to every
#'        10 iterations
#' @return list, designed to be mostly consistent with `optim` in stats.
#'         list elements:
#'            `par`: best solution
#'            `value`: function value `fn(par)`
#'            `gr`: gradient at `par`
#'            `hessian`: hessian at `par`
#'            `convergence`: convergence status, 0 if converged, 1 if not
#'                           (consistent with `optim`)
#'            `count`: a count of `fn` and `gr` evaluations
#'                     (consistent with `optim`)
#'            `trace_info`: a matrix of trace information (optional)
#' @importFrom MASS ginv
#' 
#' @export
newton_raphson <- function(
    par,
    fn,
    gr,
    hess,
    control = list(),
    ...) {
    defaults <- list(
        convergence = NULL,
        proj = function(x) x,
        tol = 1e-6,
        eta = 1,
        debug = 0L,
        r = 0.5,
        inverted = FALSE,
        fnscale = 1,
        maxit = 100L,
        trace = FALSE,
        trace_info_size_inc = 1000L,
        REPORT = 10L,
        sup = function(par) TRUE,
        test = function(par0, par1, gr0, tol, ...) {
            return(all(abs(par1 - par0) < tol) &&
                   all(abs(gr0) < tol))
        }
    )
    control <- modifyList(defaults, control)
    control <- modifyList(control, list(...))

    stopifnot(
        is.numeric(control$tol), is.numeric(control$eta),
        is.numeric(control$fnscale), is.numeric(control$maxit),
        is.numeric(control$REPORT), is.numeric(par),
        is.function(fn), is.function(gr), is.function(hess),
        is.function(control$test), is.function(control$proj),
        is.logical(control$trace), is.logical(control$inverted),
        is.function(control$sup), is.numeric(control$trace_info_size_inc),
        control$maxit > 0, control$eta > 0, control$maxit != Inf,
        control$r > 0, control$abs_tol > 0, control$REPORT > 0,
        control$fnscale != 0, control$trace_info_size_inc >= 1)

    par0 <- par
    par1 <- NULL
    val0 <- fn(par0) / control$fnscale
    hess0 <- NULL
    inv_hess0 <- NULL
    it <- 0L
    count <- 0L
    m <- length(par)
    trace_info <- NULL
    append_trace <- NULL

    if (control$trace) {
        cnames <- c("value", paste0("par", 1:m), paste0("gr", 1:m))
        trace_info <- matrix(NA,
            nrow = control$trace_info_size_inc,
            ncol = length(cnames))
        colnames(trace_info) <- cnames
        trace_app <- trace_info

        append_trace <- function(data) {
            if (count > nrow(trace_info)) {
                trace_info <<- rbind(trace_info, trace_app)
            }
            trace_info[count, ] <<- data
        }
    }

    debug_out <- function(par0, gr0, val0, it, count) {
        par_str <- paste(sprintf("%.5f", par0), collapse = ", ")
        gr_str <- paste(sprintf("%.5f", gr0), collapse = ", ")
        message(paste(
            "Iteration", sprintf("%4d", it), "|",
            "Count", sprintf("%4d", count), "|",
            "Solution: (", par_str, ")", "|",
            "Value:", sprintf("%.5f", val0), "|",
            "Gradient: (", gr_str, ")"))
    }

    par1 <- NULL

    repeat {
        gr0 <- gr(par0)
        #if (control$inverted) {
        #    inv_hess0 <- hess(par0)
        #} else {
            hess0 <- hess(par0)
        #    inv_hess0 <- ginv(hess0)
        #}
        #d0 <- inv_hess0 %*% gr0
        d0 <- ginv(hess0) %*% gr0
        #d0 <- gr0
        alpha <- control$eta
        while (it <= control$maxit) {
            it <- it + 1L
            par1 <- control$proj(par0 - alpha * d0)
            cat("it:", it, "par1: ", par1,
                ", alpha: ", alpha, ", d0: (", d0, ")'\n")
            if (!control$sup(par1)) {
                alpha <- control$r * alpha
                next
            }
            val1 <- fn(par1) / control$fnscale
            if (is.nan(val1)) {
                warning("NaN value encountered")
                next
            }
            cat("val1: ", val1, ", val0: ", val0, "\n")
            if (val1 <= val0) {
                break
            }
            alpha <- control$r * alpha
        }

        converged <- control$test(
            par0 = par0, par1 = par1, gr0 = gr0, tol = control$tol, ...)

        par0 <- par1
        val0 <- val1

        count <- count + 1L
        if (control$trace) {
            append_trace(c(control$fnscale * val0, par0, gr0))
        }

        #if (control$debug > 0L && count %% control$REPORT == 0) {
        debug_out(par0, gr0, val0, it, count)
        #}

        if (it > control$maxit || converged) {
            break
        }
    }

    res <- list(par = par0, gr = gr0, value = control$fnscale * val0,
                it = it, count = c(it, count),
                hessian = hess0, inv_hess = inv_hess0,
                convergence = as.integer(converged))

    if (control$trace) {
        res$trace_info <- trace_info[1:count, ]
    }
    res

}



















newton_raphson2 <- function(
    par,
    fn,
    gr,
    hess,
    control = list(),
    ...) {

    defaults <- list(
        proj = function(x) x,
        sup = function(par) TRUE,
        tol = 1e-6,
        eta = 1,
        r = 0.5,
        maxit = 100L
    )
    control <- modifyList(defaults, control)
    control <- modifyList(control, list(...))

    stopifnot(
        is.numeric(control$tol), is.numeric(control$eta),
        is.numeric(control$maxit),
        is.numeric(par),
        is.function(fn), is.function(gr), is.function(hess),
        is.function(control$proj),
        is.function(control$sup),
        control$maxit > 0, control$eta > 0,
        control$r > 0, control$r < 1)

    par0 <- par
    par1 <- NULL
    val0 <- fn(par0)
    it <- 0L
    repeat {
        gr0 <- gr(par0)
        d0 <- ginv(hess(par0)) %*% gr0
        alpha <- control$eta
        while (it <= control$maxit) {
            it <- it + 1L
            par1 <- control$proj(par0 - alpha * d0)
            if (!control$sup(par1)) {
                alpha <- control$r * alpha
                next
            }
            val1 <- fn(par1)
            if (is.nan(val1)) {
                warning("NaN value encountered")
                next
            }
            if (val1 >= val0) {
                break
            }
            alpha <- control$r * alpha
        }

        converged <- all(abs(par1 - par0) < control$tol)
        par0 <- par1
        val0 <- val1

        if (it > control$maxit || converged) {
            break
        }
    }

    list(par = par0, value = val1, it = it)
}
