## ------------------------------------------------------------------------
# The easy way:
binom.test(24, 29)

# The drawn-out way:
# Lower endpoint of rejection region:
qbinom(0.025, 29, 0.5) 
pbinom(9, 29, 0.5) # (9 is too big)
pbinom(8, 29, 0.5) # (8 is the lower critical value)
# 29-8 = 21 is the upper critical value

# Actual significance level = Pr(X <= 8 or X >= 21 | p = 1/2)
pbinom(8, 29, 0.5) + 1 - pbinom(20, 29, 0.5) 

# P-value for the data X = 24: 2*Pr(X >= 24 | p = 1/2)
2*(1 - pbinom(23, 29, 0.5)) 

## ------------------------------------------------------------------------
# The easy way:
binom.test(17, 25)

# The drawn-out way:
# Lower endpoint of rejection region:
qbinom(0.025, 25, 0.5) 
pbinom(8, 25, 0.5) # (8 is too big)
pbinom(7, 25, 0.5) # (7 is the lower critical value)
# 25-7 = 18 is the upper critical value

# Actual significance level = Pr(X <= 7 or X >= 18 | p = 1/2)
pbinom(7, 25, 0.5) + 1 - pbinom(17, 25, 0.5) 

# P-value for the data X = 17: 2*Pr(X >= 17 | p = 1/2)
2*(1 - pbinom(16, 25, 0.5)) 

## ------------------------------------------------------------------------
# The easy way:
binom.test(24,29)$conf.int

# The drawn-out way:
p <- seq(0, 1, by=0.001) 

# upper value:
min(p[pbinom(24, 29, p) <= 0.025])

# lower value:
max(p[1-pbinom(23, 29, p) <= 0.025]) 

## ------------------------------------------------------------------------
# The easy way:
binom.test(17,25)$conf.int 

# The drawn-out way:
p <- seq(0, 1, by=0.001) 

# upper value:
min(p[pbinom(17, 25, p) <= 0.025]) 

# lower value:
max(p[1-pbinom(16, 25, p) <= 0.025])

## ------------------------------------------------------------------------
# The easy way:
binom.test(0,15)$conf.int

# The direct way:
# lower limit = 0
# upper limit:
1-(0.025)^(1/15)

## ------------------------------------------------------------------------
# parameters

p <- c(0.01,0.05,seq(0.1,0.5,by=0.1))
n <- c(10,25,50,100,250,500,1000)

cvr <- matrix(NA,nrow=length(n),ncol=length(p))

# calculate coverage

for(i in 1:length(p)){
  for(j in 1:length(n)){
    prbs <- dbinom(0:n[j],n[j],p[i])
    z <- rep(NA,n[j]+1)
    for(k in 0:n[j]){
      tmp <- binom.test(k,n[j])$conf[1:2]
      z[k+1] <- ifelse(tmp[1]<p[i]&p[i]<tmp[2],1,0)
    }
    cvr[j,i] <- sum(z*prbs)
  }
}

# output

cvr <- data.frame(cvr)
names(cvr) <- as.character(p)
row.names(cvr) <- as.character(n)
print(round(100*cvr,1))

par(las=1,yaxs="i")
plot(range(p),range(cvr),type="n",xlab="p",ylab="coverage",xlim=c(0,0.5),ylim=c(0.94,1))
abline(h=0.95,lwd=2,col="red")
for(j in 1:length(n)) lines(p,cvr[j,],type="b",cex=0.5,pch=19)
axis(4,cvr[,length(p)],as.character(n),mgp=c(3,0.2,0),tcl=0,cex.axis=0.3)

# coverage - small n

n <- 25
x <- 7
binom.test(x,n)$conf

p <- seq(0.01,0.5,by=0.01)
cvr95 <- rep(NA,length(p))
cvr90 <- rep(NA,length(p))

for(i in 1:length(p)){
  prbs <- dbinom(0:n,n,p[i])
  z <- rep(NA,n+1)
  for(k in 0:n){
    tmp <- binom.test(k,n)$conf[1:2]
    z[k+1] <- ifelse(tmp[1]<p[i]&p[i]<tmp[2],1,0)
  }
  cvr95[i] <- sum(z*prbs)
  z <- rep(NA,n+1)
  for(k in 0:n){
    tmp <- binom.test(k,n,conf.level=0.90)$conf[1:2]
    z[k+1] <- ifelse(tmp[1]<p[i]&p[i]<tmp[2],1,0)
  }
  cvr90[i] <- sum(z*prbs)
}

par(las=1,yaxs="i")
plot(range(p),range(cvr95),type="n",xlab="p",ylab="coverage",xlim=c(0,0.5),ylim=c(0.9,1))
abline(h=0.95,lwd=2,col="red")
lines(p,cvr95,type="b",cex=0.25,pch=19)
lines(p,cvr90,type="b",cex=0.25,pch=19,col="blue")

# does it matter?

round(binom.test(x,n)$conf,2)
round(binom.test(x,n,conf.level=0.90)$conf,2)

## ------------------------------------------------------------------------
set.seed(1)
x <- 11; n1 <- 20; alpha1 <- 1; beta1 <- 1
y <- 5; n2 <- 20; alpha2 <- 1; beta2 <- 1
p1 <- rbeta(1000, x + alpha1, n1 - x + beta1)
p2 <- rbeta(1000, y + alpha2, n2 - y + beta2)
rd <- p1 - p2
quantile(rd, c(.025, .975))
mean(rd)
median(rd)
plot(density(rd))
abline(v=quantile(rd,c(0.025,0.975)),lwd=2,lty=2,col="blue")
abline(v=0,lwd=2,col="red")

