#' This function takes the output of `optim`, `newton_raphson`, or `sim_anneal`
#' and turns it into an `mle_numerical` (subclass of `mle`) object.
#'
#' @param sol the output of `optim` or `newton_raphson`
#' @param options list, options for things like sigma and FIM
#' @param superclasses list, superclasses to add to the `mle_numerical` object
#' @return An object of class \code{mle_numerical} (subclass of \code{mle}).
#' @examples
#' # Fit exponential distribution using optim
#' set.seed(123)
#' x <- rexp(100, rate = 2)
#'
#' # Log-likelihood for exponential distribution
#' loglik <- function(rate) {
#'   if (rate <= 0) return(-Inf)
#'   sum(dexp(x, rate = rate, log = TRUE))
#' }
#'
#' # Optimize (maximize by setting fnscale = -1)
#' result <- optim(
#'   par = 1,
#'   fn = loglik,
#'   method = "Brent",
#'   lower = 0.01, upper = 10,
#'   hessian = TRUE,
#'   control = list(fnscale = -1)
#' )
#'
#' # Wrap in mle_numerical
#' fit <- mle_numerical(result, options = list(nobs = length(x)))
#' params(fit)
#' se(fit)
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