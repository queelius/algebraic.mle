#' Combine independent MLEs for the same parameter.
#'
#' Given multiple independent MLEs that estimate the same parameter \eqn{\theta},
#' produces an optimally weighted combination using inverse-variance (Fisher
#' information) weighting.
#'
#' The combined estimator has:
#' \itemize{
#'   \item \code{theta.hat}: \eqn{(\sum I_i)^{-1} \sum I_i \hat\theta_i}
#'   \item \code{sigma}: \eqn{(\sum I_i)^{-1}}
#'   \item \code{info}: \eqn{\sum I_i}
#'   \item \code{nobs}: sum of individual sample sizes
#' }
#'
#' When the Fisher information matrix is not directly available but the
#' variance-covariance matrix is, the FIM is computed as \code{ginv(vcov)}.
#'
#' For the legacy interface that accepts a list, see \code{\link{mle_weighted}}.
#'
#' @param x An \code{mle} object, or a list of \code{mle} objects.
#' @param ... Additional \code{mle} objects to combine.
#' @return An \code{mle} object representing the optimally weighted combination.
#' @seealso \code{\link{mle_weighted}}, \code{\link{joint}}
#' @examples
#' # Three independent estimates of the same rate
#' fit1 <- mle(theta.hat = c(lambda = 2.1), sigma = matrix(0.04), nobs = 50L)
#' fit2 <- mle(theta.hat = c(lambda = 1.9), sigma = matrix(0.02), nobs = 100L)
#' fit3 <- mle(theta.hat = c(lambda = 2.0), sigma = matrix(0.03), nobs = 70L)
#'
#' comb <- combine(fit1, fit2, fit3)
#' params(comb)
#' se(comb)
#' @export
combine <- function(x, ...) {
    UseMethod("combine", x)
}

#' @rdname combine
#' @export
combine.list <- function(x, ...) {
    if (length(list(...)) > 0L) {
        stop("When passing a list, do not pass additional arguments.")
    }
    combine_mles(x)
}

#' @rdname combine
#' @importFrom MASS ginv
#' @importFrom algebraic.dist params nparams
#' @importFrom stats vcov nobs
#' @export
combine.mle <- function(x, ...) {
    dots <- list(...)
    if (length(dots) == 0L) return(x)
    combine_mles(c(list(x), dots))
}

#' Internal workhorse for combine.
#' @param mles A list of mle objects.
#' @return An mle object.
#' @keywords internal
combine_mles <- function(mles) {
    if (!all(vapply(mles, is_mle, logical(1)))) {
        stop("All arguments to combine() must be mle objects.")
    }
    if (length(mles) == 1L) return(mles[[1L]])

    # Get FIMs, falling back to ginv(vcov) when needed
    fims <- lapply(mles, function(m) {
        I <- observed_fim(m)
        if (!is.null(I)) return(if (is.matrix(I)) I else matrix(I, 1, 1))
        V <- vcov(m)
        if (is.null(V)) {
            stop("combine() requires either a Fisher information matrix or ",
                 "variance-covariance matrix for each mle object.")
        }
        if (!is.matrix(V)) V <- matrix(V, 1, 1)
        ginv(V)
    })

    info_combined <- Reduce(`+`, fims)
    sigma_combined <- ginv(info_combined)
    theta_combined <- as.vector(
        sigma_combined %*%
        Reduce(`+`, Map(`%*%`, fims, lapply(mles, function(m) {
            p <- params(m)
            matrix(p, ncol = 1)
        })))
    )
    names(theta_combined) <- names(params(mles[[1L]]))

    # Sum nobs (use 0L for NULL, then convert back)
    nobs_vals <- vapply(mles, function(m) {
        n <- nobs(m)
        if (is.null(n)) NA_integer_ else as.integer(n)
    }, integer(1))
    nobs_combined <- if (anyNA(nobs_vals)) NULL else sum(nobs_vals)

    mle(theta.hat = theta_combined,
        loglike = NULL,
        score = NULL,
        sigma = sigma_combined,
        info = info_combined,
        obs = NULL,
        nobs = nobs_combined)
}
