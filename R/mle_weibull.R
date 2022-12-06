#' MLE of the rate parameter when we assume the sample is i.i.d. and drawn
#' from the weibull distribution.
#'
#' @param x a sample of observations
#' @param k initial estimate of shape parameter k
#' @param eps we numerically solve the MLE equation, \code{|old-new| <= eps} is stopping condition
#' @param keep_obs Boolean, specifies whether to keep observations
#' @return an \code{mle} object.
#' @export
mle_weibull <- function(x,k=1,eps=1e-5,keep_obs=F)
{
    n <- length(x)
    stopifnot(n > 0)
    stopifnot(k > 0)
    stopifnot(eps > 0)

    s <- mean(log(x))
    repeat
    {
        k1 <- k
        k0 <- 1/(sum(x^k1*log(x))/sum(x^k1) - s)
        if (abs(k-k1) < eps)
            break
    }
    theta.hat <- c(k,mean(x^k)^(1/k))
    make_mle(
        theta.hat=theta.hat,
        loglike=weibull_loglike(x)(theta.hat),
        score=weibull_score(x)(theta.hat),
        sigma=ginv(weibull_fisher_info(x)(theta.hat)),
        info=weibull_fisher_info(x)(theta.hat),
        obs=ifelse(keep_obs,x,NULL),
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

    function(theta) matrix(c(
        n/theta[1]-n*log(theta[2]) + s -
            sum((x/theta[2])^theta[1]*log(x/theta[2])),
        theta[1]/theta[2]^(theta[1]+1)*sum(x^theta[1])-n*theta[1]/theta[2]
    ),nrow=2)
}

#' log-likelihood function generator given data \code{x} for the weibull
#' distribution
#'
#' @param x data
#' @export
weibull_fisher_info <- function(x)
{
    n <- length(x)
    function(theta)
    {
        k <- theta[1]
        lam <- theta[2]
        d <- n/lam-1/lam^(k+1)*(1+k*log(lam))*sum(x^k)-k/lam^(k+1)*sum(x^k*log(x))

        matrix(c(
            n/k^2 + sum((x/lam)^k*(log(x/lam))^2),
            d,d,
            -n*k/lam^2+k*(k+1)/lam^(k+2)*sum(x^k)
        ),nrow=2)
    }
}
