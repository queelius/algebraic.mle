
# store sequence of steps in gradient ascent/newton raphson and plot the points
# overlay it with loglike
library(tidyverse)
library(md.tools)
library(stats)

theta <- c(100,2)
n <- 17
data <- rweibull(n,shape=theta[1],scale=theta[2])
loglik <- weibull_shape_scale_loglike(data)
scr <- weibull_shape_scale_score(data)
nfo <- weibull_shape_scale_fim(data)

theta0 <- c(5,15)

sup.weibull <- function(theta) {
    all(theta > 0)
}

theta.start <- sim_anneal(
    f=loglik,
    x0=theta0,
    options=list(
        t_init=100,
        t_end=1e-4,
        alpha=0.99,
        iter_per_temp=200,
        sup=sup.weibull,
        debug=FALSE,
        trace=TRUE))

logliks <- apply(theta.start$path,1,loglik)
plot(logliks,type="l",xlab="iteration",ylab="log-likelihood")

theta.mle <- mle_weibull_shape_scale(
    data,
    k0=theta.start$argmax[1],
    eps=1e-10)

theta.nr <- mle_newton_raphson(
    ll=loglik,
    theta0=theta.start$argmax,
    score=scr,
    info=nfo,
    options=list(
        sup=sup.weibull,
        rel_tol=1e-12,
        eta=.1,
        trace=TRUE))

trace.ll <- apply(theta.nr$trace,1,loglik)
plot(trace.ll,type="l")

point(theta.nr)
mle_local_search

theta.optim <- mle_optim(optim(
    par=theta.start$argmax,
    fn=loglik,
    gr=scr,      
    hessian=TRUE,
    control=list(fnscale=-1,reltol=1e-16, maxit=2000000)))



#####################################
# TODO: bias correction, using bootstrap: theta.bias_cor = theta.hat + boot_estimate(data,mle_solver)$bias.hat
# we knwo true parameter since it's a sim, so evaluate performacne of bias corrected theta.bias_cor
# versus pure MLE estimate theta.hat.
# since we have bootstrap of MLE, we don't necessarily need to use inverse of FIM
# for variance-covariance.


library(algebraic.mle)

boot_estimate(xs,mle_normal)
bias(theta.hat,theta)


library(plyr)

pe <- function(x,y) { abs(x-y)/x*100 }

mu <- 5
sigma <- 3
theta <- c(mu,sigma)

i <- 1L
ns <- seq(1000,50000,1000)
dat <- matrix(ncol=6,nrow=length(ns))
colnames(dat) <- c("boot1","boot2","fim1","fim2","pe1","pe2")
for (n in ns)
{
    print(i)
    xs <- rnorm(n,mu,sqrt(sigma))
    theta.hat <- mle_normal(xs)

    theta.b <- raply(10000,
    {
        xs.b <- sample(xs,n,replace=T)
        theta.hat <- point(mle_normal(xs.b))
    })

    cov.boot <- cov(theta.b)
    fim <- algebraic.mle::normal_fisher_info(xs)(point(theta.hat))
    cov.fim <- ginv(algebraic.mle::normal_fisher_info(xs)(point(theta.hat)))

    dat[i,] <- c(diag(cov.boot),
                 diag(cov.fim),
                 pe(diag(cov.fim),diag(cov.boot)))
    i <- i + 1L
}

df <- as.tibble(dat)
df$n <- ns
df.tmp <- df[8:nrow(df),]

ggplot(df.tmp) +
    geom_line(aes(x=n,y=boot1),color="blue") +
    geom_line(aes(x=n,y=boot2),color="blue") +
    geom_line(aes(x=n,y=fim1),color="green") +
    geom_line(aes(x=n,y=fim2),color="green") +
    xlab("sample size") +
    ylab("bias") +
    labs(title="title",subtitle="subtitle",caption="caption")















## PIs
n <- 100
theta <- c(4,2)
x <- rnorm(n,mean=theta[1],sd=sqrt(theta[2]))
head(x,n=4)
hist(x)
theta.hat <- mle_normal_mu_var(x)
summary(theta.hat)
point(theta.hat)
fim(theta.hat)
vcov(theta.hat)
confint(theta.hat)
bias(theta.hat,theta)
bias(theta.hat)
mse(theta.hat)        # estimate of MSE
mse(theta.hat,theta)  # true MSE

mle_solver <- function(data, ind)
    point(mle_normal_mu_var(data[ind]))
R <- 100000 # number of bootstrap replicates

theta.boot <- mle_boot(mle_solver, x, R, parallel="multicore", ncpus=4)
bias(theta.boot)
bias(theta.hat)

samplr <- function(n=1,theta) rnorm(n,theta[1],theta[2])
pred(x=theta.hat, samp=samplr)
