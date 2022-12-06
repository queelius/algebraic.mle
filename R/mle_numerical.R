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
#' @export
mle_iterative <- function(
        l,
        theta0,
        dir=NULL,
        eps=1e-5,
        sup=NULL,
        eta=1,
        r=0.5,
        max_iter=0L,
        debug=F)
{
    stopifnot(eps > 0)
    if (is.null(dir)) dir <- function(theta) grad(l,theta)
    if (is.null(sup)) sup <- function(theta) all(theta > 0)

    theta1 <- theta0
    iter <- 1L
    repeat
    {
        # we use backtracking for an approximate line search
        eta0 <- eta
        d <- dir(theta0)
        l0 <- l(theta0)

        if (debug)
            cat("dir =", d, "theta0 =",theta0, ", loglike =",l0,"\n")

        repeat
        {
            theta1 <- theta0 + eta0*d

            if (debug)
                cat("theta1 =",theta1,"\n")

            if (sup(theta1))
            {
                l1 <- l(theta1)

                if (is.na(l1))
                    cat("NA:", l1, "theta1 =", theta1,"\n")

                if (!is.na(l1) && l1 >= l0)
                    break
            }
            eta0 <- r*eta0
            if (iter == max_iter)
                break

            iter <- iter + 1L
        }

        if (iter == max_iter || abs(max(theta1-theta0)) < eps)
            break

        theta0 <- theta1
    }

    theta.hat <- make_mle(
        theta.hat=theta1,
        loglike=l(theta1),
        superclasses=c("mle_numerical"))
    theta.hat$iter <- iter
    theta.hat$max_iter <- max_iter
    theta.hat$learning_rate <- eta
    theta.hat
}

#' Iterative MLE method using the Newton-Raphson method
#'
#' @param l log-likelihood function of type \code{R^p -> R}
#' @param theta0 initial guess of theta with \code{p} components
#' @param info information matrix function of type \code{R^p -> R^{p-by-q}}, defaults to taking the negative of the hessian of \code{l}.
#' @param score score function of type \code{R^p -> R^p}, defaults to taking the gradient of \code{l}.
#' @param inverted if T, then \code{info} is inverted (covariance matrix instead of info)
#' @param eps stopping condition
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
        eps=1e-5,
        sup=NULL,
        eta=1,
        r=0.5,
        max_iter=0L,
        debug=F,
        inverted=F)
{
    res <- mle_iterative(
        l=l,
        theta0=theta0,
        dir=ifelse(inverted,
                   function(theta) info(theta) %*% score(theta),
                   function(theta) ginv(info(theta)) %*% score(theta)),
        eps=eps,
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
#' @param eps stopping condition
#' @param sup domain of support function of type \eqn{R^p \mapsto \{true,false\}} for the log-likelihood function \eqn{l}
#' @param eta learning rate, defaults to 1
#' @param r backing tracking parameter
#' @param max_iter maximum iterations
#' @param debug Boolean, output debugging information if TRUE; defaults to FALSE
#'
#' @importFrom MASS ginv
#' @export
mle_gradient_ascent <- function(
        l,
        theta0,
        score,
        eps=1e-5,
        sup=NULL,
        eta=1,
        r=0.5,
        max_iter=0L,
        debug=F)
{
    res <- mle_iterative(
        l=l,
        theta0=theta0,
        dir=score,
        eps=eps,
        sup=sup,
        eta=eta,
        r=r,
        max_iter=max_iter,
        debug=debug)

    res$info <- -hessian(l,res$theta.hat)
    #res$info <- -jacobian(score,res$theta.hat)
    res$score <- score(res$theta.hat)
    res$sigma <- ginv(res$info)
    res
}





#' Iterative MLE method using gradient ascent
#'
#' @param mle_solver MLE solver, takes a starting point
#' @param start_gen starting point generator
#' @param trials number of trials to run
#' @export
mle_solver_random_restarts <- function(
        mle_solver,
        start_gen,
        trials=10)
{
    mle <- NULL
    loglik <- -Inf

    for (i in 1:trials)
    {
        mle.candidate <- mle_solver(start_gen())
        l <- loglike(mle.candidate)
        if (l > loglik)
        {
            loglik <- l
            mle <- mle.candidate
        }
    }

    mle
}
