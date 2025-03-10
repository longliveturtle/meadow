## ------------------------------------------------------------------------
xA <- c(27.0,54.6,33.5,27.6,46.0,22.2,44.2,17.3,15.9,32.8)
xB <- c(17.4,20.5,13.9,14.8,27.9,10.6,33.7,15.4,25.0,24.1)

stripchart(list(xA,xB),vertical=T,pch=21,xlim=c(0.5,2.5))
segments(0.9,mean(xA),1.1,mean(xA),lwd=2,col="red")
segments(1.9,mean(xB),2.1,mean(xB),lwd=2,col="blue")
abline(h=mean(c(xA,xB)),lty=2)

t.test(xA,xB)
wilcox.test(xA,xB)

## ------------------------------------------------------------------------
sort(xA)
sort(xA,decreasing=T)
order(xA)
order(xA,decreasing=T)

xAnew <- xA
xAnew[2] <- 99.9
xA
xAnew

stripchart(list(xAnew,xB),vertical=T,pch=21,xlim=c(0.5,2.5))
segments(0.9,mean(xAnew),1.1,mean(xAnew),lwd=2,col="red")
segments(1.9,mean(xB),2.1,mean(xB),lwd=2,col="blue")
abline(h=mean(c(xAnew,xB)),lty=2)

t.test(xA,xB)$p.value
t.test(xAnew,xB)$p.value # why is the p-value now > 0.05?
wilcox.test(xA,xB)$p.value
wilcox.test(xAnew,xB)$p.value # why are the p-values the same?

## ------------------------------------------------------------------------
my.perm <- function(xa,xb){
  # combine the two vectots
  x <- c(xa,xb)
  # calculate the total length of observations
  n <- length(x)
  # calculate the length of the first vector
  na <- length(xa)
  # shuffle the order of the observations across grups
  z <- x[sample(1:n)]
  # return the t test statistic from the shuffled data
  return(t.test(z[1:na],z[-(1:na)],var.equal=T)$stat)
}

## ------------------------------------------------------------------------
# Let's try one with our xA and xB
my.perm(xA,xB)
# And onther one
my.perm(xA,xB)
# And onther one
my.perm(xA,xB)

# Let's do 1000 shuffles
nit <- 1000
# Create a vector of length nit (1000)
myres <- rep(NA,nit)
# Fill in the slots of that vector with the permutation test t statistics
for(j in 1:nit) myres[j] <- my.perm(xA,xB)

# Plot the data
hist(myres,breaks=21,col="lightgrey")
# Add the observed test statistic to the plot
tst <- t.test(xA,xB,var.equal=T)$stat 
abline(v=tst,col="red")
# Calculate the permutation p-value
mean(abs(myres)>abs(tst))

