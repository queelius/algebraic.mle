
# store sequence of steps in gradient ascent/newton raphson and plot the points
# overlay it with loglike


theta <- c(100,2)
data <- rweibull(n,shape=theta[1],scale=theta[2])
loglik <- weibull_loglike(data)
scr <- weibull_score(data)
nfo <- weibull_fisher_info(data)

theta0 <- c(5,15)
good <- algebraic.mle::mle_weibull(data,k0=100)
center <- point(good)
good2 <- algebraic.mle::mle_gradient_ascent(l=loglik,theta0=theta0,score=scr,debug=F)
bad <- algebraic.mle::mle_newton_raphson(l=loglik,theta0=theta0,score=scr,info=nfo,debug=T)



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
