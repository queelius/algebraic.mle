#' Linearly transform an \code{mle} object \code{x} by multiplying it
#' on the LHS by a matrix \code{A}.
#'
#' @param A a (non-random) matrix
#' @param x an \code{mle} object to linearly transform
#' @return a \code{mle_linear_transform} object
#' @export
mle_linear_transform <- function(A,x)
{
    t.theta.hat <- A %*% point(x)
    t.sigma <- t(A) %*% vcov(x) %*% A

    # the origional observations are `nobs(x)`, which is how `x` was
    # fitted as an mle object.  The transformed observations are
    # `t.obs` are the observations that were needed to yield an
    # mle `A %*% x`. NOTE: TEST THIS OUT!
    t.obs <- lapply(obs(x), function(o) A %*% o)

    mle(
        theta.hat = t.theta.hat,
        loglike = NULL,
        score = NULL,
        sigma = t.sigma,
        info = ginv(t.sigma),
        obs = t.obs,
        nobs = nobs(x),
        superclasses = c("mle_linear_transform"))
}