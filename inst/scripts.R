
# store sequence of steps in gradient ascent/newton raphson and plot the points
# overlay it with loglike


theta <- c(100,2)
n <- 17
data <- rweibull(n,shape=theta[1],scale=theta[2])
loglik <- weibull_shape_scale_loglike(data)
scr <- weibull_shape_scale_score(data)
nfo <- weibull_shape_scale_fim(data)

theta0 <- c(5,15)
k0 <- theta0[1]

sup.weibull <- function(theta) {
    k <- theta[1]
    b <- theta[2]
    if (k <= 0) {
        return(FALSE)
    }
    if (b <= 0) {
        return(FALSE)
    }
    TRUE
}
proj.weibull <- function(theta) {
    k <- theta[1]
    b <- theta[2]
    if (k <= 0) {
        k <- 1e-3
    }
    if (b <= 0) {
        b <- 1e-3
    }
    c(k,b)
}

theta.sa <- sim_anneal(
    f=loglik,
    x0=theta0,
    t_init=100,
    alpha=.99,
    sup=sup.weibull,
    iter_per_temp=100,
    t_end=1e-12,
    trace=TRUE,
    debug=FALSE)

sa.optim <- optim(par=theta0,fn=loglik,method="SANN",control=list(maxit=2000000,
    fnscale=-1,temp=100,trace=0,reltol=1e-12, tol=1e-12, maxeval=20000000, maxtime=Inf))

z <- theta.sa$trace_fx
x <- theta.sa$trace_x[,1]
y <- theta.sa$trace_x[,2]
df <- tibble(x=x, y=y, z=z)

threshold <- quantile(df$z, 0.05) # removes the lowest 5% of data
df <- df[df$z > threshold, ]

# let's plot the z vs iteration
ggplot(df, aes(x=1:length(z), y=z)) +
  geom_line() +
  theme_minimal() +
  xlab("iteration") +
  ylab("loglike")

# plot the heatmap
library(ggplot2)
ggplot(df, aes(x=x, y=y, color=z)) +
  geom_point() +
  theme_minimal() +
    xlab("shape") +
    ylab("scale") +
    scale_color_gradient(low="blue", high="red")

theta.mle <- mle_weibull_shape_scale(
    data,
    k0=theta.sa$argmax[1],
    eps=1e-16)

theta.ga <- mle_gradient_ascent(
    ll=loglik,
    theta0=theta.sa$argmax,
    score=scr,
    debug=FALSE,
    proj=proj.weibull,
    sup=sup.weibull,
    eps=1e-16,
    eta=1,
    max_iter=1000L)

theta.nr <- mle_newton_raphson(
    ll=loglik,
    theta0=theta.sa$argmax,
    score=scr,
    info=nfo,
    proj=proj.weibull,
    sup=sup.weibull,
    eps=1e-20,
    eta=1,
    max_iter=1000L,
    debug=FALSE)

theta.optim <- mle_optim(optim(
    par=theta.sa$argmax,
    fn=loglik,
    gr=scr,      
    hessian=TRUE,
    control=list(fnscale=-1,reltol=1e-16, maxit=2000000)))

options(digits=20)
names(theta.sa$max) <- "sa"
l.mle <- loglike(theta.mle)
names(l.mle) <- "mle"
l.optim <- loglike(theta.optim)
names(l.optim) <- "optim"
l.nr <- loglike(theta.nr)
names(l.nr) <- "nr"
l.ga <- loglike(theta.ga)
names(l.ga) <- "ga"

sort(c(theta.sa$max,l.mle,l.optim,l.nr,l.ga))

##########
# ggplot #
##########

library(ggplot2)

delta <- 0.75
diff <- delta / 100
as <- seq(center[1]-delta/2,center[1]+delta/2,diff)
bs <- seq(center[2]-delta/2,center[2]+delta/2,diff)

df <- tibble(x=numeric(),y=numeric(),z=numeric())

N <- length(as) * length(bs)
t <- numeric(length=N)
u <- numeric(length=N)
v <- numeric(length=N)

i <- 1L
for (a in as)
{
    for (b in bs)
    {
        t[i] <- a
        u[i] <- b
        v[i] <- loglik(c(a,b))
        i <- i + 1L
    }
}

#sqrt(t(scr(center))%*%scr(center))

d <- scr(center)
mag <- sqrt(sum(d^2))
d.begin <- center
d.end <- center + .01*(d/mag)


d2 <- MASS::ginv(nfo(center)) %*% scr(center)
mag2 <- sqrt(sum(d2^2))
d2.begin <- center
d2.end <- center + .01*(d2/mag2)


df <- tibble(x=t,y=u,z=v)
ggplot(df, aes(x=x,y=y,z=z)) + geom_contour(bins=100,aes(colour=after_stat(level))) +
    geom_point(aes(x=center[1],y=center[2]),colour="red") +
    geom_segment(aes(x=d.begin[1],
                     y=d.begin[2],
                     xend=d.end[1],
                     yend=d.end[2]),
                 arrow = arrow(length = unit(0.2, "cm"))) +
    geom_segment(aes(x=d.begin[1],
                     y=d.begin[2],
                     xend=d2.end[1],
                     yend=d2.end[2]),
                 arrow = arrow(length = unit(0.2, "cm")))
##########
# plotly #
##########
library(plotly)

center <- point(good2)
delta <- 1.5
diff <- delta / 1000
as <- seq(center[1]-delta/2,center[1]+delta/2,diff)
bs <- seq(center[2]-delta/2,center[2]+delta/2,diff)
zs <- matrix(nrow=length(as),ncol=length(bs))

i <- 1L
j <- 1L
for (a in as)
{
    for (b in bs)
    {
        zs[i,j] <- loglik(c(a,b))
        j <- j + 1L
    }
    j <- 1L
    i <- i + 1L
}

plot_ly(x=as,y=bs,z=zs,type="surface")
#plot_ly(x=as,y=bs,z=zs,type="contour")




###


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















##



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
