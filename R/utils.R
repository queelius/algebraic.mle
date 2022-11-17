confidence_intervals <- function(V, theta, parm=NULL, level=0.95, ...)
{
    sigma <- diag(V)
    p <- length(theta)
    q <- stats::qnorm(level)
    if (is.null(parm))
        parm <- 1:p

    parm <- parm[parm >= 1 & parm <= p]
    ci <- matrix(nrow=length(parm),ncol=2)
    colnames(ci) <- c(paste((1-level)/2*100,"%"),
                      paste((1-(1-level)/2)*100,"%"))

    i <- 1
    for (j in parm)
    {
        ci[i,] <- c(theta[j] - q*sqrt(sigma[j]),
                    theta[j] + q*sqrt(sigma[j]))
        i <- i + 1
    }
    rownames(ci) <- parm
    ci
}

