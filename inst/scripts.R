
# store sequence of steps in gradient ascent/newton raphson and plot the points
# overlay it with loglike
library(tidyverse)
library(md.tools)
library(stats)

theta <- c(100,2)
n <- 17
data <- rweibull(n,shape=theta[1],scale=theta[2])
loglik <- weibull_loglike(data)
scr <- weibull_score(data)
nfo <- weibull_fim(data)

theta0 <- c(5,15)

sup.weibull <- function(theta) {
    all(theta > 0)
}

theta.start <- sim_anneal(
    fn=loglik,
    par=theta0,
    control=list(
        t_init=150,
        t_end=1e-9,
        fnscale=-1,
        alpha=0.995,
        it_per_temp=225,
        REPORT=250L,
        maxit=1000000L,
        sup=sup.weibull,
        debug=1L,
        trace=TRUE))

theta.start$par
# Load required library
library(ggplot2)


m <- 2
# Convert the matrix to a data frame
data_df <- as_tibble(theta.start$trace_info)
data_df$best <- as.factor(data_df$best)
# Reshape data for ggplot2
library(reshape2)
long_data <- melt(data_df,
    id.vars = c("it", "value", "temp", "best"),
    variable.name = "parameter",
    value.name = "value_par")

# Convergence plot (1)
convergence_plot <- ggplot(data_df, aes(x = it, y = value, color = best)) +
  geom_line() +
  labs(title = "Convergence Plot",
       x = "Iteration",
       y = "Best Function Value") +
  scale_color_discrete(name = "Best Value", labels = c("No", "Yes"))

print(convergence_plot)

# Parameter traces plot (3)
parameter_traces_plot <- ggplot(long_data, aes(x = it, y = value_par, color = best)) +
  geom_line() +
  facet_wrap(~parameter, ncol = m, scales = "free_y") +
  labs(title = "Parameter Traces",
       x = "Iteration",
       y = "Parameter Value") +
  scale_color_discrete(name = "Best Value", labels = c("No", "Yes"))

print(parameter_traces_plot)


logliks <- -theta.start$trace_info[,"value"]
iters <- theta.start$trace_info[,"it"]
plot(logliks,type="l",xlab="iteration",ylab="log-likelihood")

best.trace <- theta.start$trace_info[theta.start$trace_info[,"best"] == 1,]
plot(-best.trace[,"value"],type="l",xlab="iteration",ylab="log-likelihood")

theta.mle <- mle_weibull(
    data,
    k0=theta.start$par[1],
    eps=1e-11)

theta.nr2 <- newton_raphson(
    par=theta.start$par,
    fn=loglik,
    gr=scr,
    hess=function(x) -nfo(x),
    #hess=function(x) numDeriv::hessian(loglik, x,
    #    method.args=list(eps=1e-5, d=0.00001, r=6)),
    control=list(
        eta=1,
        abs_tol=1e-6,
        rel_tol=1e-6,
        fnscale=-1,
        proj=function(x) pmax(x,1e-3),
        trace=FALSE,
        fn=loglik,
        r=.5,
        debug=1L,
        REPORT=1000L,
        trace_info_size_inc=10000L,
        inverted=FALSE,
        maxit=100))

point(theta.mle)
theta.optim$par
theta.nr2$par
theta.nr$par
theta.start$par

par <- c(1,2)
-nfo(par)
numDeriv::hessian(loglik, par, method.args=list(eps=1e-5, d=0.00001, r=6))
scr(par)
numDeriv::grad(loglik, par)
loglik(par)
sum(dweibull(data,shape=par[1],scale=par[2],log=TRUE))

theta.optim <- optim(
    par=theta.start$par,
    fn=loglik,
    gr=scr,      
    hessian=TRUE,
    control=list(fnscale=-1,reltol=1e-16, maxit=2000000))

confint(mle_numerical(theta.optim))
confint(mle_numerical(theta.nr2))
confint(theta.mle)




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
