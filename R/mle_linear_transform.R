#' Linearly transform an \code{mle} object \code{x} by multiplying it
#' on the LHS by a matrix \code{A}.
#'
#' @param A a (non-random) matrix
#' @param x an \code{mle} object to linearly transform
#' @export
mle_linear_transform <- function(A,x)
{
    mu <- A %*% point(x)
    sigma <- t(A) %*% point(x) %*% A

    res <- structure(list(
        theta.hat=mu,
        info=MASS::ginv(sigma),
        sigma=sigma),
        class=unique(c("linear_transformed_mle",class(x))))
    attributes(res) <- attributes(x)
    res
}
