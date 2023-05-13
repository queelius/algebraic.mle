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
mle_numerical <- function(sol, options = list()) {
    if (is.null(options$hessian)) {
        if (!is.null(options$sigma)) {
            options$hessian <- -ginv(options$sigma)
        } else if (!is.null(sol$sigma)) {
            options$hessian <- sol$hessian
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
        info = -options$hessian,
        obs = options$obs,
        nobs = options$nobs,
        superclasses = c("mle_numerical"))
    sol_mle$converged <- sol$convergence == 0
    sol_mle
}