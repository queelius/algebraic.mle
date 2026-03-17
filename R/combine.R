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
#' @param x An \code{mle_fit} object, or a list of \code{mle_fit} objects.
#' @param ... Additional \code{mle_fit} objects to combine.
#' @return An \code{mle_fit} object representing the optimally weighted combination.
#' @seealso \code{\link{joint}}
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
combine.mle_fit <- function(x, ...) {
    dots <- list(...)
    if (length(dots) == 0L) return(x)
    combine_mles(c(list(x), dots))
}

#' Internal workhorse for combine.
#' @param mles A list of mle_fit objects.
#' @return An mle_fit object.
#' @keywords internal
combine_mles <- function(mles) {
    if (!all(vapply(mles, is_mle, logical(1)))) {
        stop("All arguments to combine() must be mle_fit objects.")
    }
    if (length(mles) == 1L) return(mles[[1L]])

    dims <- vapply(mles, nparams, integer(1))
    if (length(unique(dims)) > 1L) {
        stop("All MLEs must have the same number of parameters for combine(). ",
             "Got: ", paste(dims, collapse = ", "))
    }

    fims <- lapply(mles, fim_or_ginv_vcov)

    info_combined <- Reduce(`+`, fims)
    sigma_combined <- ginv(info_combined)

    # Inverse-variance weighted combination: (sum I_i)^{-1} sum(I_i theta_i)
    weighted_sum <- Reduce(`+`, Map(
        function(I, m) I %*% matrix(params(m), ncol = 1),
        fims, mles))
    theta_combined <- as.vector(sigma_combined %*% weighted_sum)
    names(theta_combined) <- names(params(mles[[1L]]))

    nobs_vals <- vapply(mles, function(m) {
        n <- nobs(m)
        if (is.null(n)) NA_integer_ else as.integer(n)
    }, integer(1))
    nobs_combined <- if (anyNA(nobs_vals)) NULL else sum(nobs_vals)

    mle(theta.hat = theta_combined,
        sigma = sigma_combined,
        info = info_combined,
        nobs = nobs_combined)
}

#' Extract FIM from an MLE, falling back to ginv(vcov).
#'
#' @param m An mle_fit object.
#' @return A matrix (FIM).
#' @keywords internal
fim_or_ginv_vcov <- function(m) {
    I <- observed_fim(m)
    if (!is.null(I)) {
        if (!is.matrix(I)) I <- matrix(I, 1, 1)
        return(I)
    }
    V <- vcov(m)
    if (is.null(V)) {
        stop("combine() requires either a Fisher information matrix or ",
             "variance-covariance matrix for each mle_fit object.")
    }
    if (!is.matrix(V)) V <- matrix(V, 1, 1)
    ginv(V)
}
