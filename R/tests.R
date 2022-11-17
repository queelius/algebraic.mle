#' MLE of the (mu,sigma) parameter vector when we assume the sample is i.i.d. and
#' drawn from the normal distribution.
#'
#' Of course, the draws are unlikely to be normal, but the normal distribution
#' is often a good model. Hypothesis testing, such as relative
#' likelihoods, can be used to assess the appropriateness of the normal
#' model to the data.
#'
#' @param obs a sample of observations
#' @param obs a sample of observations from reference distribution
#' @param nbreaks number of bins
#' @return a hypothesis test object
#' @export
chisqr.test <- function(obs,
                        cdf,
                        nbreaks=ceiling(sqrt(length(obs))),
                        ...)
{
    h <- hist(obs,right=F,breaks=nbreaks)
    t <- stats::chisq.test(
        h$counts,
        p=zoo::rollapply(cdf(h$breaks,...),
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
#' @export
cramer.test <- function(obs,ref)
{
    stat <- CDFt::CramerVonMisesTwoSamples(obs,ref)
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
        return(CDFt::KolmogorovSmirnov(obs,ref))
    }
    else if (method=="chisqr")
    {
        stopifnot(!is.null(cdf))
        return(chisqr.test(obs,cdf,...))
    }
}
