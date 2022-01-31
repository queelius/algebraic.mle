#' Fisher scoring algorithm.
#'
#' @param theta0 initial guess of theta with \eqn{p} components
#' @param info information matrix function of type \eqn{R^p -> R^{p \times q}}
#' @param score score function of type \eqn{R^p -> R^p}
#' @param eps stopping condition
#' @param max_iter maximum number of iterations
#' @export
fisher_scoring <- function(theta0,info,score,eps=1e-5,max_iter=250L)
{
    theta1 <- theta0
    nfo <- NULL
    sigma <- NULL

    for (iter in 1L:max_iter)
    {
        nfo <- info(theta0)
        sigma <- MASS::ginv(nfo)
        s <- score(theta0)
        theta1 <- theta0 + sigma %*% s
        if (max(abs(s)) < eps)
            break
        theta0 <- theta1
    }

    structure(list(
        theta.hat=theta1,
        sigma=sigma,
        score=s,
        info=nfo,
        eps=eps,
        iter=iter,
        max_iter=max_iter),
        class=c("mle_numerical","mle","estimate"))
}


#' Gradient ascent
#'
#' @param theta0 initial guess of theta with \eqn{p} components
#' @param loglike log-likelihood function
#' @param score score function of type \eqn{R^p -> R^p}
#' @param eps stopping condition
#' @param max_iter maximum iterations
#' @export
gradient_ascent <- function(theta0,
                            loglike,
                            score=NULL,
                            eps=1e-3,
                            max_iter=250L)
{
    if (is.null(score))
        score <- function(theta) numDeriv::grad(loglike,theta)

    i <- 1L
    repeat
    {
        s <- score(theta0)
        theta1 <- theta0 + s
        if (i == max_iter || max(abs(s) < eps))
        {
            structure(list(
                theta.hat=theta0,
                score=s,
                eps=eps,
                iter=i),
                class=c("mle_numerical","mle","estimate"))
        }


        #alpha <- stats::optimise(ls,c(0,max_alpha),maximum=T,p=p)$maximum
        theta0 <- theta1
        i <- i + 1L
    }
}
