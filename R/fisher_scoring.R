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

    for (iter in 1:max_iter)
    {
        nfo <- info(theta0)
        sigma <- MASS::ginv(nfo)
        s <- score(theta0)
        theta1 <- theta0 + sigma %*% s
        if (max(abs(theta1-theta0)) < eps)
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
#' @param normalize whether to normalize the jumps to a magnitude of max_alpha
#' @param max_alpha maximum alpha step
#' @param randomize randomize steps according to some energy function that
#'                  decreases over time
#' @export
gradient_ascent <- function(theta0,
                            loglike,
                            score=NULL,
                            eps=1e-3,
                            normalize=F,
                            max_alpha=3,
                            randomize=F)
{
    if (is.null(score))
        score <- function(theta) numDeriv::grad(loglike,theta)

    energy <- function(iter) { exp(-iter) }
    norm <- function(x) { x / sqrt(sum(x^2)) }
    ls <- function(t,p) { loglike(theta0 + t * p) }
    #stop_cond <- function(s) { sqrt(sum(s^2)) < eps }
    stop_cond <- function(s) { max(abs(s)) < eps }
    proj <- function(theta)
    {
        for (j in 1:length(theta))
        {
            if (theta[j] < 0)
                theta[j] <- stats::runif(1,1,5)
        }
        theta
    }

    i <- 1
    repeat
    {
        p <- score(theta0)
        #if (stop_cond(p))

        #if (randomize) p <- p + energy(i)*stats::rnorm(n=length(theta0))
        #if (normalize) p <- norm(p)

        theta1 <- proj(theta0 + p)
        if (stop_cond(theta0-theta1))
        {
            structure(list(
                theta.hat=theta0,
                score=p,
                eps=eps,
                iter=i),
                class=c("mle_numerical","mle","estimate"))
        }


        #alpha <- stats::optimise(ls,c(0,max_alpha),maximum=T,p=p)$maximum
        #theta0 <- proj(theta0 + alpha * p)
        i <- i + 1
    }
}
