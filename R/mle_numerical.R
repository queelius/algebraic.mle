#' General iterative MLE method
#'
#' @param l log-likelihood function of type \eqn{R^p \mapsto R}
#' @param theta0 initial guess of \eqn{\theta} with \eqn{p} components
#' @param dir promising direction function of type \eqn{R^p \mapsto R^p}
#' @param stop_cond stopping condition function of type \eqn{R^p \times R^p \mapsto \{T,F\}}
#' @param sup domain of support function of type \eqn{R^p \mapsto \{T,F\}} for the log-likelihood function \eqn{l}
#' @param eta learning rate, defaults to 1
#' @param r backing tracking parameter
#' @param max_iter maximum number of iterations
#' @param debug Boolean, output debugging information if TRUE; defaults to FALSE
#'
#' @importFrom numDeriv grad
#' @importFrom numDeriv hessian
#' @export
mle_iterative <- function(
        l,
        theta0,
        dir=NULL,
        stop_cond=NULL,
        sup=NULL,
        eta=1,
        r=0.5,
        max_iter=0L,
        debug=F)
{
    if (is.null(dir)) dir <- function(theta) grad(l,theta)
    if (is.null(stop_cond)) stop_cond <- function(theta1,theta0) abs(max(theta1-theta0)) < 1e-3
    if (is.null(sup)) sup <- function(theta) all(theta > 0)

    theta1 <- theta0
    iter <- 1L
    repeat
    {
        # we use backtracking for an approximate line search
        eta0 <- eta
        d <- dir(theta0)
        l0 <- l(theta0)
        repeat
        {
            theta1 <- theta0 + eta0*d
            if (sup(theta1) && l(theta1) > l0)
                break
            eta0 <- r*eta0
        }

        if (debug)
            print(theta1)

        if (iter == max_iter || stop_cond(theta1,theta0))
            break

        iter <- iter + 1L
        theta0 <- theta1
    }

    theta.hat <- make_mle(theta1,l(theta1))
    theta.hat$iter <- iter
    theta.hat$max_iter <- max_iter
    theta.hat$learning_rate <- eta
    class(theta.hat) <- unique(c("mle_numerical",class(theta.hat)))
    theta.hat
}

#' Iterative MLE method using the Newton-Raphson method
#'
#' @param l log-likelihood function of type \code{R^p -> R}
#' @param theta0 initial guess of theta with \code{p} components
#' @param info information matrix function of type \code{R^p -> R^{p-by-q}}
#' @param score score function of type \code{R^p -> R^p}
#' @param stop_cond stopping condition function of type \code{(R^p,R^p) -> \{T,F\}}
#' @param sup domain of support function of type \code{R^p -> \{T,F\}} for the log-likelihood function \code{l}
#' @param eta learning rate, defaults to 1
#' @param r backing tracking parameter
#' @param max_iter maximum number of iterations
#' @param debug Boolean, output debugging information if TRUE; defaults to FALSE
#'
#' @importFrom MASS ginv
#' @export
mle_newton_raphson <- function(
        l,
        theta0,
        info,
        score,
        stop_cond=NULL,
        sup=NULL,
        eta=1,
        r=0.5,
        max_iter=0L,
        debug=F)
{
    res <- mle_iterative(
        l,
        theta0,
        dir=function(theta) ginv(info(theta)) %*% score(theta),
        stop_cond=stop_cond,
        sup=sup,
        eta=eta,
        r=r,
        max_iter=max_iter,
        debug=debug)

    res$score <- score(res$theta.hat)
    res$info <- info(res$theta.hat)
    res$sigma <- ginv(res$info)
    res
}

#' Iterative MLE method using gradient ascent
#'
#' @param theta0 initial guess of theta with \eqn{p} components
#' @param l log-likelihood function
#' @param score score function of type \eqn{R^p -> R^p}
#' @param stop_cond stopping condition function of type \eqn{R^p \times R^p \mapsto \{T,F\}}
#' @param sup domain of support function of type \eqn{R^p \mapsto \{true,false\}} for the log-likelihood function \eqn{l}
#' @param eta learning rate, defaults to 1
#' @param r backing tracking parameter
#' @param max_iter maximum iterations
#' @param debug Boolean, output debugging information if TRUE; defaults to FALSE
#'
#' @importFrom MASS ginv
#' @importFrom numDeriv hessian
#' @importFrom numDeriv grad
#' @export
mle_gradient_ascent <- function(
        l,
        theta0,
        score=function(theta) grad(l,theta),
        stop_cond=NULL,
        sup=NULL,
        eta=1,
        r=0.5,
        max_iter=0L,
        debug=F)
{
    res <- mle_iterative(
        l,
        theta0,
        dir=score,
        stop_cond=stop_cond,
        sup=sup,
        eta=eta,
        r=r,
        max_iter=max_iter,
        debug=debug)

    res$score <- score(res$theta.hat)
    res$info <- -hessian(l,res$theta.hat)
    res$sigma <- ginv(res$info)
    res
}
