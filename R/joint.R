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
#' @param x An \code{mle} object.
#' @param ... Additional \code{mle} objects to join.
#' @return An \code{mle} object representing the joint MLE.
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

#' @rdname joint
#' @importFrom algebraic.dist params nparams
#' @importFrom stats vcov nobs
#' @importFrom MASS ginv
#' @export
joint.mle <- function(x, ...) {
    mles <- c(list(x), list(...))

    if (length(mles) < 2L) {
        stop("joint() requires at least 2 mle objects.")
    }
    if (!all(vapply(mles, is_mle, logical(1)))) {
        stop("All arguments to joint() must be mle objects.")
    }

    # Validate disjoint parameter names
    all_names <- unlist(lapply(mles, function(m) names(params(m))))
    if (anyDuplicated(all_names)) {
        dups <- all_names[duplicated(all_names)]
        stop("Parameter names must be disjoint. Overlapping: ",
             paste(unique(dups), collapse = ", "))
    }

    # All must have vcov
    vcovs <- lapply(mles, vcov)
    if (any(vapply(vcovs, is.null, logical(1)))) {
        stop("All mle objects must have a variance-covariance matrix for joint().")
    }

    # Concatenate parameters
    theta_joint <- unlist(lapply(mles, params))

    # Build block-diagonal vcov
    dims <- vapply(mles, nparams, integer(1))
    p <- sum(dims)
    sigma_joint <- matrix(0, p, p)
    offset <- 0L
    for (i in seq_along(mles)) {
        idx <- offset + seq_len(dims[i])
        V <- vcovs[[i]]
        if (!is.matrix(V)) V <- matrix(V, 1, 1)
        sigma_joint[idx, idx] <- V
        offset <- offset + dims[i]
    }

    # Sum log-likelihoods (NULL if any missing)
    loglikes <- lapply(mles, loglik_val)
    loglike_joint <- if (any(vapply(loglikes, is.null, logical(1)))) {
        NULL
    } else {
        sum(unlist(loglikes))
    }

    # Block-diagonal FIM (NULL if any missing)
    fims <- lapply(mles, observed_fim)
    info_joint <- if (any(vapply(fims, is.null, logical(1)))) {
        NULL
    } else {
        info <- matrix(0, p, p)
        offset <- 0L
        for (i in seq_along(mles)) {
            idx <- offset + seq_len(dims[i])
            I <- fims[[i]]
            if (!is.matrix(I)) I <- matrix(I, 1, 1)
            info[idx, idx] <- I
            offset <- offset + dims[i]
        }
        info
    }

    # Concatenate scores (NULL if any missing)
    scores <- lapply(mles, score_val)
    score_joint <- if (any(vapply(scores, is.null, logical(1)))) {
        NULL
    } else {
        unlist(scores)
    }

    mle(theta.hat = theta_joint,
        loglike = loglike_joint,
        score = score_joint,
        sigma = sigma_joint,
        info = info_joint,
        obs = NULL,
        nobs = NULL)
}
