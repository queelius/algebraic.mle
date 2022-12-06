library(plyr)

do.boots <- function(n)
{
    xs <- rnorm(n,0,1)

    var <- sum((xs - mean(xs))^2)/n

    var.b <- raply(10000,
    {
        xs.b <- sample(xs, n, replace=T)
        sum((xs.b - mean(xs.b))^2)/n
    })

    # estimated bias
    list(bias.hat=mean(var.b)-var, bias=((n-1)/n)-1)
}


perc_error <- function(x,theta) abs((x$theta.hat-theta)/theta)*100



xs <- seq(100, 1000, 100)
ys <- numeric(length(xs))
i <- 1L
for (x in xs) { ys[i] <- do.boots(x)$bias.hat; i <- i + 1L }

bias.true <- (xs-1)/xs-1

