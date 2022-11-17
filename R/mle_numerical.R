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
            if (sup(theta1) && l(theta1) >= l0)
                break
            eta0 <- r*eta0
        }

        if (debug)
            cat("theta =",theta1, ", loglike =",l(theta1),"\n")

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
#' @param info information matrix function of type \code{R^p -> R^{p-by-q}}, defaults to taking the negative of the hessian of \code{l}.
#' @param score score function of type \code{R^p -> R^p}, defaults to taking the gradient of \code{l}.
#' @param inverted if T, then \code{info} is inverted (covariance matrix instead of info)
#' @param stop_cond stopping condition function of type \code{(R^p,R^p) -> \{T,F\}}
#' @param sup domain of support function of type \code{R^p -> \{T,F\}} for the log-likelihood function \code{l}
#' @param eta learning rate, defaults to 1
#' @param r backing tracking parameter
#' @param max_iter maximum number of iterations
#' @param debug Boolean, output debugging information if TRUE; defaults to FALSE
#'
#' @importFrom MASS ginv
#' @importFrom numDeriv grad
#' @importFrom numDeriv hessian
#' @importFrom numDeriv jacobian
#' @export
mle_newton_raphson <- function(
        l,
        theta0,
        info=NULL,
        score=NULL,
        stop_cond=NULL,
        sup=NULL,
        eta=1,
        r=0.5,
        max_iter=0L,
        debug=F,
        inverted=F)
{
    stopifnot(!inverted || !is.null(info))

    if (is.null(info))
        info <- ifelse(is.null(score),
                       function(theta) -hessian(l,theta),
                       function(theta) -jacobian(score,theta))
    if (is.null(score))
        score <- function(theta) grad(l,theta)

    res <- mle_iterative(
        l,
        theta0,
        dir=ifelse(inverted,
                   function(theta) info(theta) %*% score(theta),
                   function(theta) ginv(info(theta)) %*% score(theta)),
        stop_cond=stop_cond,
        sup=sup,
        eta=eta,
        r=r,
        max_iter=max_iter,
        debug=debug)

    res$score <- score(res$theta.hat)
    if (inverted)
    {
        res$sigma <- info(res$theta.hat)
        res$info <- ginv(res$sigma)
    }
    else
    {
        res$info <- info(res$theta.hat)
        res$sigma <- ginv(res$info)
    }
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
        score=NULL,
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
        dir=ifelse(is.null(score),function(theta) grad(l,theta),score),
        stop_cond=stop_cond,
        sup=sup,
        eta=eta,
        r=r,
        max_iter=max_iter,
        debug=debug)

    res$info <- ifelse(is.null(score),
                       -hessian(l,res$theta.hat),
                       -jacobian(score,res$theta.hat))
    res$score <- ifelse(is.null(score),
                        grad(l,res$theta.hat),
                        score(res$theta.hat))
    res$sigma <- ginv(res$info)
    res
}
