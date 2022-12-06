#' MLE of the (mu,sigma) parameter vector when we assume the sample is i.i.d. and
#' drawn from the normal distribution.
#'
#' Of course, the draws are unlikely to be normal, but the normal distribution
#' is often a good model. Hypothesis testing, such as relative
#' likelihoods, can be used to assess the appropriateness of the normal
#' model to the data.
#'
#' @param x a sample of observations
#' @param keep_obs whether to store observations in \code{mle_normal} object
#' @return an \code{mle} object.
#' @export
mle_normal <- function(x,keep_obs=F)
{
    n <- length(x)
    stopifnot(n > 0)
    theta.hat <- mle_normal_mu_var(x)
    info.hat <- normal_fisher_info(x)(theta.hat)

    make_mle(
        theta.hat=theta.hat,
        loglike=normal_loglike(x)(theta.hat),
        score=normal_score(x)(theta.hat),
        sigma=ginv(info.hat),
        info=info.hat,
        obs=ifelse(keep_obs,x,NULL),
        sample_size=n,
        superclasses=c("mle_normal"))
}

#' MLE of the (mu,sigma) parameter vector when we assume the sample is i.i.d. and
#' drawn from the normal distribution.
#'
#' @param x a sample of observations
#' @export
mle_normal_mu_var <- function(x)
{
    mu <- mean(x)
    sigma <- mean((x-mu)^2)
    matrix(c(mu,sigma),nrow=2)
}

#' log-likelihood function generator given data \code{x} for the normal
#' distribution
#'
#' @param x data
#' @export
normal_loglike <- function(x)
{
    # more efficiently computes:
    #    sum(dnorm(x,mean=theta[1],sd=sqrt(theta[2]),log=T))
    n <- length(x)
    S <- sum(x)
    S2 <- sum(x^2)
    function(theta) -n/2*log(2*pi*theta[2]) -
        1/(2*theta[2])*(S2-2*theta[1]*S+n*theta[1]^2)
}

#' Score function generator given data \code{x} for the normal
#' distribution
#'
#' @param x data (simple random sample)
#' @export
normal_score <- function(x)
{
    n <- length(x)
    S <- sum(x)
    S2 <- sum(x^2)
    function(theta)
        matrix(c(0.5/theta[2]*(S-n*theta[1]),
                 -0.5*n/theta[2]+0.5/theta[2]^2*(S2-2*theta[1]*S+n*theta[1]^2)),
               nrow=2)
}


#' Fisher information matrix generator given data \code{x} for the normal
#' distribution
#'
#' @param x data (simple random sample)
#' @param observed whether to return the observed fisher information, default is TRUE
#' @export
normal_fisher_info <- function(x,observed=T)
{
    n <- length(x)
    if (observed)
    {
           S <- sum(x)
           S2 <- sum(x^2)
           return(function(theta)
           {
               covar <- 0.5/theta[2]^2*(S-n*theta[1])
               matrix(c(n/theta[2],
                        covar,
                        covar,
                        -0.5*n/theta[2]^2+1/theta[2]^3*(S2-2*theta[1]*S+n*theta[1]^2)),
                      nrow=2)
           })
    }
    else
    {
        return(function(theta)
            matrix(c(n/theta[2],
                     0,0,
                     0.5*n/theta[2]^2),
                   nrow=2))
    }
}

#' Computes the bias of an \code{mle_normal} object.
#'
#' @param x the \code{mle} object to compute the bias of.
#' @param theta the true parameter value
#' @export
bias.mle_normal <- function(x,theta)
{
    matrix(c(0,-1/x$sample_size*theta[2]),nrow=2)
}


