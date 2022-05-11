#' Linearly transform an \code{mle} object \code{x} by multiplying it
#' on the LHS by a matrix \code{A}.
#'
#' @param A a matrix
#' @param x an \code{mle} object to linearly transform
#' @importFrom matlib %*%
#' @export
mle_linear_transform(A,x)
{
    mu <- A %*% point(x)
    sigma <- t(A) %*% point(x) %*% A

    res <- structure(list(
        theta.hat=mu,
        info=MASS::ginv(sigma),
        sigma=sigma),
        class=unique(c("transformed_mle",class(x))))
    attributes(res) <- attributes(x)
    res
}
