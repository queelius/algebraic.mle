#' chi-squar goodness-of-fit test
#'
#' @param obs a sample of observations
#' @param cdf the cdf of the reference distribution
#' @param options  a list of options to pass to the test function
#' @param breaks the number of breaks to use in the histogram
#' @param rescale.p logical, whether to rescale the p-value
#' @param simulate.p.value logical, whether to simulate the p-value
#' @param nbreaks the number of breaks to use in the histogram
#' @param ... additional arguments to pass
#' @return a hypothesis test object
#'
#' @importFrom stats chisq.test
#' @importFrom zoo rollapply
#' @export
chisqr.test <- function(obs,
                        cdf,
                        options = list(),
                        ...)
{
    defaults <- list(
        breaks = NULL,
        rescale.p = TRUE,
        simulate.p.value = TRUE,
        nbreaks = ceiling(sqrt(length(obs))))

    options <- modifyList(defaults,options)
    options <- modifyList(options,list(...))

    if (is.null(options$breaks)) {
        breaks = options$nbreaks
    }
    else {
        breaks = "Sturges"
    }

    h <- hist(obs, right = FALSE, breaks = breaks)
    t <- chisq.test(
        h$counts,
        p=rollapply(cdf(h$breaks,...), 2,
                         \(x) x[2]-x[1]),
        rescale.p = options$rescale.p,
        simulate.p.value = options$simulate.p.value)

    t$data.name <- deparse(substitute(obs))
    t$size <- length(obs)
    t$nbreaks <- nbreaks
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
two.sample.test <- function(sampler,obs,cdf=NULL,options=list(),...)
{
    defaults <- list(
        n=1000,
        method="cramer")

    options <- modifyList(defaults,options)
    options <- modifyList(options,list(...))

    if (options$method=="cramer")
    {
        stopifnot(!is.null(sampler))
        ref <- sampler(options$n,...)
        return(cramer.test(obs,ref))
    }
    else if (options$method=="ks")
    {
        stopifnot(!is.null(sampler))
        ref <- sampler(options$n,...)
        return(ks.test(obs,ref))
    }
    else if (options$method=="ks2")
    {
        stopifnot(!is.null(sampler))
        ref <- sampler(n,...)
        return(KolmogorovSmirnov(obs,ref))
    }
    else if (options$method=="chisqr")
    {
        stopifnot(!is.null(cdf))
        return(chisqr.test(obs=obs,cdf=cdf,...))
    }
}
