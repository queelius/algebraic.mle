#' Fisher scoring algorithm.
#'
#' @param theta0 initial guess of theta with \eqn{p} components
#' @param info information matrix function of type \eqn{R^p -> R^{p \times q}}
#' @param score score function of type \eqn{R^p -> R^p}
#' @param eps stopping condition
#' @param max_iter maximum number of iterations
#'
#' @importFrom dplyr %>%
#'
#' @export
fisher_scoring <- function(theta0,info,score,eps=1e-5,max_iter=250L)
{
    theta1 <- theta0
    nfo <- NULL
    sigma <- NULL

    for (iter in 1:max_iter)
    {
        print("---")
        nfo <- info(theta0)
        print(nfo)
        sigma <- MASS::ginv(nfo)
        print(sigma)
        s <- score(theta0)
        print(s)
        theta1 <- theta0 + sigma %*% s
        print(theta0)
        if (max(abs(theta1-theta0)) < eps)
            break
        theta0 <- theta1
        print("###")
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

