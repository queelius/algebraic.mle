#' Log-Likelihood Function Generator for the Exponential Distribution
#' with a rate parameter that is a function of predictors in `df`,
#'
#'    `resp(df) ~ exp(rate(df, beta))`
#'
#' where `beta` is a vector of parameters.
#' 
#' The function returned by this function is suitable for use with
#' `optim` or `nlm` to find the maximum likelihood estimate of `beta`.
#' 
#' Note that the expected value is `1/rate(...)`.
#' 
#' @param df a numeric vector representing a sample of observations.
#' @param resp a function that returns the response variable given a row
#'             from a data frame.
#' @param rate a function that returns the rate given a row from a data frame
#'             and a parameter vector `beta`
#' @return Returns a function that computes the conditional log-likelihood
#' @export
exp_conditional_rate_loglike <- function(df, resp, rate) {

    function(beta) {
        sum(dexp(x = resp(df), rate = rate(df, beta), log = TRUE))
    }
}



#' Log-Likelihood Function Generator for the Weibull Distribution
#' with a shape and scale parameters that are a function of predictors in `df`,
#'
#'    `resp(df) ~ weibull(shape(df, beta), scale(df, beta))`
#'
#' where `beta` is a vector of parameters.
#' 
#' The function returned by this function is suitable for use with
#' `optim` or `nlm` to find the maximum likelihood estimate of `beta`.
#' 
#' Note that the expected value given the data is
#' `scale(...) * gamma(1 + 1/shape(...))`.
#' 
#' @param df a numeric vector representing a sample of observations.
#' @param resp a function that returns the response variable given a data frame.
#' @param shape a function that returns the shape given a data frame
#'             and a parameter vector `beta`
#' @param scale a function that returns the scale given a data frame and a
#'              parameter vector `beta`
#' @return Returns a function that computes the conditional log-likelihood
#' @export
weibull_conditional_shape_scale_loglike <-
    function(df, resp, shape, scale) {
    function(beta) {
        sum(dweibull(x = resp(df),
            shape = shape(df, beta),
            scale = scale(df, beta), log = TRUE))
    }
}