######################################################################
# Computer notes                                Biostatistics 140.751
# Homework 1                                 Johns Hopkins University 
######################################################################
# Lines beginning with the symbol '#' are comments in R.  All other
# lines contain code.
#
# In R for Windows, you may wish to open this file from the menu bar
# (File:Display file); you can then easily copy commands into the
# command window.  (Use the mouse to highlight one or more lines; then
# right-click and select "Paste to console".)
######################################################################

###############################################
# 2d
###############################################

# This generates a somewhat fancy plot. Try to 
# figure what all these commands do!

func1=function(n){
z=2/(n-1)
return(z)}

func2=function(n){
z=(3*n-4)/((n-1)^2)
return(z)}

par(las=1)
n=10:100
plot(n,func2(n),pch=0,ylim=range(c(0,func2(n))),
     xlab="number of observations",ylab="variance",
     main=expression(paste("Comparison of the estimators ",S^2," and Q")))
points(n,func1(n),pch=15)
abline(h=0)
axis(1,seq(10,100,10))
axis(4,seq(0,0.3,0.05),rep("",7))
text(40,0.1,expression(bold(Q)))
text(15,0.1,expression(bold(S^2)))
text(90,0.3,expression(paste("We assumed that ",sigma^2,"=1")))


###############################################
# 3b
###############################################

# Here is the code for the iteration:
a=matrix(c(3/5,2/5,4/5,1/5),ncol=2,byrow=T);b=a;eps=0.1;cnt=0
while(eps>10^(-6)){
cnt=cnt+1;cat(cnt,"")
b.old=b;b=b%*%a
eps=max(abs(b-b.old))}

# Print the matrix
print(b)


###############################################
# 4c
###############################################

# Here is the code for the density and contour plots
mu=c(8,0)
sigma=matrix(c(30,-2,-2,2),ncol=2)

## 3D density plot of the bivariate normal distribution 
x1=seq(-10,26,0.5)
x2=seq(-3,3,0.1)
x=expand.grid(x1,x2)

z=solve(sigma)%*%(t(x)-mu)
z=(t(x)-mu)*z
z=apply(z,2,sum)
z=1/(2*pi*sqrt(det(sigma)))*exp(-0.5*z)
z=matrix(z,ncol=length(x2))

persp(x1,x2,z,theta=30,phi=30,expand=0.5,
      xlab="X1",ylab="X2",zlab="",
      col = "lightblue",main="",ticktype="d",
      ltheta = 120) 

## contour plot of the bivariate normal distribution 
x1=seq(-10,26,0.2)
x2=seq(-3,3,0.04)
x=expand.grid(x1,x2)

z=solve(sigma)%*%(t(x)-mu)
z=(t(x)-mu)*z
z=apply(z,2,sum)
z=1/(2*pi*sqrt(det(sigma)))*exp(-0.5*z)
z=matrix(z,ncol=length(x2))

par(las=1)
image(x1,x2,z,col=grey(1-seq(0.1,0.8,length=256)),
      xlab=expression(X[1]),ylab=expression(X[2]))
contour(x1,x2,z,drawlabels=FALSE,add=TRUE,levels=10^(-6:-2))
title(expression(paste("The joint distribution of ",X[1]," and ",X[2]),sep=""))
text(8,0,expression(paste(mu," = (8,0)",sep="")),col="white",cex=1.5)

## 3D density plot of the conditional distribution 
x1=seq(-10,26,0.5)
x2=seq(-3,3,0.1)
z=matrix(ncol=length(x1),nrow=length(x2))
for (j in 1:length(x2)) z[j,]=dnorm(x1,mean=8-x2[j],sd=sqrt(28))

persp(x1,x2,t(z),theta=30,phi=30,expand=0.5,
      xlab="X1",ylab="X2",zlab="",
      col = "lightblue",main="",ticktype="d",
      ltheta = 120) 

## contour plot of the conditional distribution 
x1=seq(-10,26,0.2)
x2=seq(-3,3,0.04)
z=matrix(ncol=length(x1),nrow=length(x2))
for (j in 1:length(x2)) z[j,]=dnorm(x1,mean=8-x2[j],sd=sqrt(28))

par(las=1)
image(x1,x2,t(z),col=grey(1-seq(0.1,0.8,length=256)),
      xlab=expression(X[1]),ylab=expression(x[2]))
lines(8-x2,x2,col="white")
title(expression(paste("The distribution of ",X[1],"|",X[2]==x[2]),sep=""))
text(11,0,expression(paste(x[1]==8-x[2]),sep=""),col="white",cex=1.5)


##################
# End of hw1.R
##################
