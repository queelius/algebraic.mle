% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/algebraic.mle.R
\name{observed_fim}
\alias{observed_fim}
\title{Generic method for computing the observed FIM
of an `mle` object.}
\usage{
observed_fim(x, ...)
}
\arguments{
\item{x}{the object to obtain the fisher information of}

\item{...}{additional arguments to pass}
}
\description{
Fisher information is a way of measuring the amount of
information that an observable random variable `X`
carries about an unknown parameter `theta`
upon which the probability of `X` depends.
}
\details{
The inverse of the Fisher information matrix
is the variance-covariance of the MLE for
`theta`.

Some MLE objects do not have an observed FIM,
e.g., if the MLE's sampling distribution was
bootstrapped.
}
