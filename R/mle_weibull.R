#' MLE of the rate parameter when we assume the sample is i.i.d. and drawn
#' from the weibull distribution.
#'
#' @param x a sample of observations
#' @return an \code{mle} object.
#' @export
mle_weibull <- function(x,eps=1e-5)
{
    n <- length(x)
    stopifnot(n > 0)

    k.hat <- 1
    s <- mean(log(x))
    repeat
    {
        k.old <- k.hat
        k.hat <- (sum(x^k.old*log(x))/sum(x^k.old) - s)^(-1)
        if (abs(k.hat-k.old) < eps)
            break
    }
    theta.hat <- c(k.hat,mean(x^k.hat)^(1/k.hat))
    make_mle(
        theta.hat=theta.hat,
        loglike=weibull_loglike(x)(theta.hat),
        score=weibull_score(x)(theta.hat),
        sigma=ginv(weibull_fisher_info(x)(theta.hat)),
        info=weibull_fisher_info(x)(theta.hat),
        obs=x,
        sample_size=n)
}

#' log-likelihood function generator given data \code{x} for the weibull
#' distribution
#'
#' @param x data
#' @export
weibull_loglike <- function(x)
{
    n <- length(x)
    s <- sum(log(x))
    function(theta) n*log(theta[1]) - n*theta[1]*log(theta[2]) +
        (theta[1]-1)*s - sum(x^theta[1])/theta[2]^theta[1]
}

#' score function generator given data \code{x} for the weibull
#' distribution given a simple random sample.
#'
#' @param x data
#' @export
weibull_score <- function(x)
{
    n <- length(x)
    s <- sum(log(x))

    function(theta) c(
        n/theta[1]-n*log(theta[2]) + s -
            sum((x/theta[2])^theta[1]*log(x/theta[2])),
        -n*theta[1]/theta[2] + theta[1]/theta[2]^(theta[1]+1)*sum(x^theta[1])
    )
}

#' log-likelihood function generator given data \code{x} for the weibull
#' distribution
#'
#' @param x data
#' @export
weibull_fisher_info <- function(x)
{
}
