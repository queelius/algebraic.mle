#' Function for obtaining a prediction for a function \code{f} with
#' respect to an \code{mle} object \code{x}.
#'
#' We use Monte-Carlo simulation to get \code{N} sample points within the
#' \code{p}-confidence region and taking the maximum value with respect to
#' some distance measure.
#'
#' @param N the sample size
#' @param x the \code{mle} object
#' @param p the value that defines the confidence region
#' @param f the function to predict the value of with respect to input \code{x}
#' @param D the distance measure, defaults to identity function
#'
#' @export
pred <- function(x,f,N,p=.95,D=function(x) x)
{
    stopifnot(N >= 1)
    stopifnot(inherits(x,"mle"))
    stopifnot(p > 0 && p <= 1)

    mu <- point(x)
    xs <- as.matrix(sample_mle_region(N-1,x,p))
    R <- nrow(xs)

    y.max <- f(mu)
    y.min <- y.max
    d.max <- D(y.max)
    d.min <- d.max

    for (i in seq_len(R))
    {
        y <- f(xs[i,])
        d.trial <- D(y)
        if (d.max < d.trial)
        {
            d.max <- d.trial
            y.max <- y
        }
        if (d.min > d.trial)
        {
            d.min <- d.trial
            y.min <- y
        }
    }

    c(y.min,y.max)
}

#' Function for obtaining an empirical sampling distribution from the asymptotic
#' distribution of the \code{mle} object \code{x} when applied to a function
#' \code{g}.
#'
#' @param n the sample size
#' @param x the \code{mle} object
#' @param g the function to predict the value of with respect to input \code{x}
#' @param alpha confidence-level, (1-alpha)
#' @param ... pass additional arguments
#'
#' @export
pred.interval <- function(x,g,n=1000,alpha=.05,...)
{
    stopifnot(is.integer(n) && n >= 1L)
    stopifnot(is_mle(x))
    stopifnot(is.function(g))

    g.samp <- apply(as.matrix(sampler(x,...)(n),nrow=n),1,g)
    g.hat <- g(point(x))

    k <- ncol(g.samp)
    intervals <- matrix(nrow=k,ncol=3)
    for (j in 1:k)
        intervals[j,] <- c(g.hat[j],quantile(sort(g.samp[,j]),c(alpha/2,1-alpha/2)))
    colnames(intervals) <- c("mle","lower","upper")
    intervals
}

