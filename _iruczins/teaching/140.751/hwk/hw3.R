######################################################################
# Computer notes                                Biostatistics 140.751
# Homework 3                                 Johns Hopkins University 
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
# 4
###############################################

myrmvn=function(mu,sigma,hm=1,...){
n=length(mu)
if(sum((dim(sigma)-rep(n,2))^2)!=0) 
stop("Check the dimensions of mu and sigma!")
if(det(sigma)==0) stop("The covariance matrix is singular!")
a=t(chol(sigma))
z=matrix(rnorm(n*hm),nrow=n)
y=t(a%*%z+mu)
return(y)}

mymu=c(14,-7);mysigma=matrix(c(1,-0.75,-0.75,1),2)
postscript("fig.ps")
par(las=1)
plot(myrmvn(mymu,mysigma,hm=1000),
xlab=expression(y[1]),ylab=expression(y[2]))
dev.off()

###############################################
# 5a
###############################################

myrchisq=function(n,df,lambda){
x=matrix(rnorm(n*df),ncol=df)
x=x+sqrt(2*lambda/df)
x=x^2
x=apply(x,1,sum)
return(x)}

###############################################
# 5b
###############################################

mypchisq=function(q,n,df,lambda) sum(myrchisq(n,df,lambda)<=q)/n

##################
# End of hw3.R
##################
