#' chi-squar goodness-of-fit test
#'
#' @param obs a sample of observations
#' @param cdf the cdf of the reference distribution
#' @param nbreaks number of bins
#' @param ... additional arguments to pass
#' @return a hypothesis test object
#'
#' @importFrom stats chisq.test
#' @importFrom zoo rollapply
#' @export
chisqr.test <- function(obs,
                        cdf,
                        nbreaks=ceiling(sqrt(length(obs))),
                        ...)
{
    h <- hist(obs,right=F,breaks=nbreaks)
    t <- chisq.test(
        h$counts,
        p=rollapply(cdf(h$breaks,...),
                         2,
                         \(x) x[2]-x[1]),
        rescale.p=T,
        simulate.p.value=T)

    t$data.name=deparse(substitute(obs))
    t$size=length(obs)
    t$nbreaks=nbreaks

    t
}

#' Cramer Von Mises hypothesis test
#'
#' @param obs a sample of observations
#' @param ref a sample from the reference distribution
#' @return hypothesis test object
#' @importFrom CDFt CramerVonMisesTwoSamples
#' @export
cramer.test <- function(obs,ref)
{
    stat <- CramerVonMisesTwoSamples(obs,ref)
    structure(list(
        p.value=exp(-stat)/6.0,
        stat=stat,
        obs.size=length(obs),
        ref.size=length(ref)),
        class=c("cramer.test","hypothesis.test"))
}

#' General two sample hypothesis test.
#'
#' @param sampler random sampler for reference distribution
#' @param obs a sample of observations
#' @param cdf a cdf for the reference distribution, defaults to NULL
#' @param n sample size (for sampler)
#' @param method the method for the two sample hypothesis test
#' @return hypothesis test object
#' @importFrom CDFt KolmogorovSmirnov
#' @export
two.sample.test <- function(sampler,obs,cdf=NULL,n=10000,method="cramer")
{
    if (method=="cramer")
    {
        stopifnot(!is.null(sampler))
        ref <- sampler(n,...)
        return(cramer.test(obs,ref))
    }
    else if (method=="ks")
    {
        stopifnot(!is.null(sampler))
        ref <- sampler(n,...)
        return(ks.test(obs,ref))
    }
    else if (method=="ks2")
    {
        stopifnot(!is.null(sampler))
        ref <- sampler(n,...)
        return(KolmogorovSmirnov(obs,ref))
    }
    else if (method=="chisqr")
    {
        stopifnot(!is.null(cdf))
        return(chisqr.test(obs,cdf,...))
    }
}
