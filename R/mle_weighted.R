#' mle_weighted takes a list of MLEs for some parameter gamma and then
#' combines them to form the MLE for gamma over all the information in the
#' samples used to generate the list of MLEs.
#'
#' @param mles A list of MLEs, all for the same parameter.
#' @return an object of type \code{mle_weighted}, which inherits from \code{mle}.
#' @export
mle_weighted <- function(mles)
{
    infos <- purrr::map(mles,fisher_info)
    info.wt <- purrr::reduce(infos,`+`)
    theta.wt <- purrr::reduce(infos %*% mles,`+`)
    cov.wt <- MASS::ginv(info.wt)
    structure(list(
        theta.hat=theta.wt,
        info=info.wt,
        sigma=cov.wt),
        class=c("mle_weighted","mle","estimate"))
}
