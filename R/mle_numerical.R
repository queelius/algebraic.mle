#' mle_optim
#'
#' This function takes the output of `optim`, `newton_raphson`, or `sim_anneal`
#' and turns it into an `mle_numerical` (subclass of `mle`) object.
#' 
#' @param sol the output of `optim` or `newton_raphson`
#' @param options list, options for things like sigma and FIM
#' @return a `numerical_mle` object.
#' @importFrom MASS ginv
#' @export
mle_numerical <- function(sol, options = list(), superclasses = NULL) {
    if (is.null(options$hessian)) {
        if (!is.null(sol$hessian)) {
            options$hessian <- sol$hessian
        } else if (!is.null(options$sigma)) {
            options$hessian <- -ginv(options$sigma)
        }
    }
    if (is.null(options$sigma) && !is.null(options$hessian)) {
        options$sigma <- -ginv(options$hessian)
    }
    sol_mle <- mle(
        theta.hat = sol$par,
        loglike = sol$value,
        score = NULL,
        sigma = options$sigma,
        info = if (is.null(options$hessian)) NULL else -options$hessian,
        obs = options$obs,
        nobs = options$nobs,
        superclasses = c(superclasses, "mle_numerical"))
    sol_mle$converged <- sol$convergence == 0
    sol_mle$sol <- sol
    sol_mle$sol$call <- match.call()
    sol_mle
}