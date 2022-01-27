#' \code{mle_weighted} takes a list of MLEs for some parameter theta and then
#' combines them to form the MLE for theta over all the information in the
#' samples used to generate the MLEs in the list.
#'
#' @param mles A list of MLEs, all for the same parameter.
#' @return an object of type \code{mle_weighted}, which inherits from \code{mle}.
#' @export
mle_weighted <- function(mles)
{
    A <- fisher_info(mles[[1]])
    B <- A %*% point(mles[[1]])
    info.wt <- A

    for (i in 2:length(mles))
    {
        A <- fisher_info(mles[[i]])
        info.wt <- info.wt + A
        B <- B + A %*% point(mles[[i]])
    }

    cov.wt <- solve(info.wt,diag(1,nrow(info.wt)))
    #cov.wt <- MASS::ginv(info.wt)
    theta.wt <- cov.wt %*% B

    structure(list(
        theta.hat=theta.wt,
        info=info.wt,
        sigma=cov.wt),
        class=c("mle_weighted","mle","estimate"))
}
