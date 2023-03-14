library(algebraic.mle)
library(plyr)
library(MASS)
library(tibble)
library(ggplot2)

pe <- function(obs,expected) { abs((obs-expected)/expected)*100 }

mu <- 5
sigma <- 3
theta <- c(mu,sigma)

i <- 1L
ns <- c(25,50,100,250,500,seq(1000,100000,1000))
dat <- matrix(ncol=15,nrow=length(ns))
colnames(dat) <- c("err(mu)","err(var)",
                   "bias(mu)","bias(var)",
                   "bias(mu).boot","bias(var).boot",
                   "pe(mu)","pe(var)",
                   "E{pe(mu)}","E{pe(var)}",
                   "var(mu).boot","var(var).boot",
                   "var(mu).fim","var(var).fim",
                   "n")
for (n in ns)
{
    xs <- rnorm(n,mu,sqrt(sigma))

    theta.hat <- mle_normal(xs)
    theta.b <- raply(20000,
    {
        point(mle_normal(sample(xs,n,replace=T)))
    })

    err <- point(theta.hat) - theta
    bias.actual <- c(0,-sigma/n)
    bias.hat <- colMeans(theta.b) - point(theta.hat)
    pe.obs <- pe(point(theta.hat),theta)
    pe.expected <- bias.actual/theta
    var.theta.boot <- diag(cov(theta.b))
    var.theta.fim <- diag(ginv(algebraic.mle::normal_fisher_info(xs)(point(theta.hat))))

    dat[i,] <- c(err,
                 bias.actual,
                 bias.hat,
                 pe.obs,
                 pe.expected,
                 var.theta.boot,
                 var.theta.fim,
                 n)
    print(dat[i,])
    i <- i + 1L
}

df <- as_tibble(dat)

ggplot(df[3:10,]) +
    geom_line(aes(x=n,y=`var(var).boot`,color="`var(var).boot`")) +
    geom_line(aes(x=n,y=`var(var).fim`,color="`var(var).fim`")) +
    xlab("sample size") +
    ylab("estimate") +
    labs(title="Estimate of variance of variance of MLE for normal",
         subtitle="using Bootstrap and Fisher information methods")

ggplot(df[3:10,]) +
    geom_line(aes(x=n,y=`var(mu).boot`,color="`var(mu).boot`")) +
    geom_line(aes(x=n,y=`var(mu).fim`,color="`var(mu).fim`")) +
    xlab("sample size") +
    ylab("estimate") +
    labs(title="Estimate of variance of mean of MLE for normal",
         subtitle="using Bootstrap and Fisher information methods")

ggplot(df[3:10,]) +
    geom_line(aes(x=n,y=`E{pe(var)}`,color="E{pe(var)}")) +
    #geom_line(aes(x=n,y=`pe(var)`,color="pe(var)")) +
    xlab("sample size") +
    ylab("percent error") +
    labs(title="Percent error of MLE estimator for normal distribution")

ggplot(df) +
    geom_line(aes(x=n,y=`bias(mu)`,color="`bias(mu)`")) +
    geom_line(aes(x=n,y=`bias(mu).boot`,color="`bias(mu).boot`")) +
    xlab("sample size") +
    ylab("bias") +
    labs(title="Bias of MLE estimator of mean for normal distribution")

ggplot(df[10:105,]) +
    geom_line(aes(x=n,y=`bias(var)`,color="`bias(var)")) +
    geom_line(aes(x=n,y=`bias(var).boot`,color="`bias(var).boot`")) +
    xlab("sample size") +
    ylab("bias") +
    labs(title="Bias of MLE estimator of variance for normal distribution")

ggplot(df) +
    geom_line(aes(x=n,y=`err(mu)`,color="err(mu)")) +
    geom_line(aes(x=n,y=`err(var)`,color="err(var)")) +
    xlab("sample size") +
    ylab("error") +
    labs(title="Error of MLE estimator for normal distribution")



ggplot(df[10:105,]) +
    geom_line(aes(x=n,y=`bias(var)`,color="`bias(var)")) +
    #geom_line(aes(x=n,y=`bias(var).boot`,color="`bias(var).boot`")) +
    geom_ma(aes(x=n,y=`bias(var).boot`,color="`ma(bias(var).boot)`")) +
    #geom_smooth(aes(x=n,y=`bias(var).boot`,color="`smooth(bias(var).boot)`")) +
    xlab("sample size") +
    ylab("bias") +
    labs(title="Bias of MLE estimator of variance for normal distribution")

