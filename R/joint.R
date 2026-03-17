#' Compose independent MLEs into a joint MLE.
#'
#' Given two or more independent MLEs with disjoint parameter sets, produces
#' a joint MLE with block-diagonal variance-covariance structure.
#'
#' The joint MLE has:
#' \itemize{
#'   \item \code{theta.hat}: concatenation of all parameter vectors
#'   \item \code{sigma}: block-diagonal from individual vcov matrices
#'   \item \code{loglike}: sum of log-likelihoods (when all available)
#'   \item \code{info}: block-diagonal from individual FIMs (when all available)
#'   \item \code{score}: concatenation of score vectors (when all available)
#'   \item \code{nobs}: NULL (different experiments have no shared sample size)
#' }
#'
#' @param x An \code{mle_fit} object.
#' @param ... Additional \code{mle_fit} objects to join.
#' @return An \code{mle_fit} object representing the joint MLE.
#' @examples
#' # Two independent experiments
#' fit_rate <- mle(theta.hat = c(lambda = 2.1), sigma = matrix(0.04), nobs = 50L)
#' fit_shape <- mle(theta.hat = c(k = 1.5, s = 3.2),
#'                  sigma = matrix(c(0.1, 0.02, 0.02, 0.3), 2, 2), nobs = 100L)
#'
#' # Joint MLE: 3 params, block-diagonal covariance
#' j <- joint(fit_rate, fit_shape)
#' params(j)   # c(lambda = 2.1, k = 1.5, s = 3.2)
#' vcov(j)     # 3x3 block-diagonal
#'
#' # Existing algebra works on the joint:
#' marginal(j, 2:3)   # recover shape params
#' as_dist(j)          # MVN for distribution algebra
#' @export
joint <- function(x, ...) {
    UseMethod("joint", x)
}

#' Build a block-diagonal matrix from a list of square matrices.
#'
#' @param blocks List of matrices (scalars are promoted to 1x1 matrices).
#' @param dims Integer vector of block dimensions.
#' @return A square matrix with blocks placed along the diagonal.
#' @keywords internal
block_diag <- function(blocks, dims) {
    p <- sum(dims)
    result <- matrix(0, p, p)
    offset <- 0L
    for (i in seq_along(blocks)) {
        idx <- offset + seq_len(dims[i])
        B <- blocks[[i]]
        if (!is.matrix(B)) B <- matrix(B, 1, 1)
        result[idx, idx] <- B
        offset <- offset + dims[i]
    }
    result
}

#' Collect a field from a list of MLEs, returning NULL if any are NULL.
#'
#' @param mles List of mle_fit objects.
#' @param extractor Function that extracts the field from a single MLE.
#' @return A list of field values, or NULL if any value is NULL.
#' @keywords internal
collect_or_null <- function(mles, extractor) {
    vals <- lapply(mles, extractor)
    if (any(vapply(vals, is.null, logical(1)))) NULL else vals
}

#' @rdname joint
#' @importFrom algebraic.dist params nparams
#' @importFrom stats vcov nobs
#' @importFrom MASS ginv
#' @export
joint.mle_fit <- function(x, ...) {
    mles <- c(list(x), list(...))

    if (length(mles) < 2L) {
        stop("joint() requires at least 2 mle_fit objects.")
    }
    if (!all(vapply(mles, is_mle, logical(1)))) {
        stop("All arguments to joint() must be mle_fit objects.")
    }

    all_names <- unlist(lapply(mles, function(m) names(params(m))))
    if (anyDuplicated(all_names)) {
        dups <- all_names[duplicated(all_names)]
        stop("Parameter names must be disjoint. Overlapping: ",
             paste(unique(dups), collapse = ", "))
    }

    vcovs <- lapply(mles, vcov)
    if (any(vapply(vcovs, is.null, logical(1)))) {
        stop("All mle_fit objects must have a variance-covariance matrix for joint().")
    }

    dims <- vapply(mles, nparams, integer(1))

    loglikes <- collect_or_null(mles, function(m) m$loglike)
    scores <- collect_or_null(mles, score_val)
    fims <- collect_or_null(mles, observed_fim)

    mle(theta.hat = unlist(lapply(mles, params)),
        loglike = if (!is.null(loglikes)) sum(unlist(loglikes)),
        score = if (!is.null(scores)) unlist(scores),
        sigma = block_diag(vcovs, dims),
        info = if (!is.null(fims)) block_diag(fims, dims),
        obs = NULL,
        nobs = NULL)
}
