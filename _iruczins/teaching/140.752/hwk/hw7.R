######################################################################
# Computer notes                                 Biostatistics 140.752
# Homework 7                                  Johns Hopkins University 
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
# 2
###############################################

# (a)
hm=101
f=seq(0,2,length=hm)
mypw=as.list(1:4)
for(p in 2:5){
  x=matrix(ncol=4,nrow=hm)
  for(k in 1:4){
    n=5*k
    cv=qf(0.95,p-1,n*p-p)
    x[,k]=1-pf(cv/(1+n*f),p-1,n*p-p)
  }
  mypw[[(p-1)]]=x
}
#
par(las=1,mfrow=c(2,2),xaxs="i",yaxs="i")
for (k in 1:4){
  plot(c(0,2),c(0,1),xlab=expression(sigma[1]^2/sigma^2),ylab="power",type="n")
  par(cex.axis=0.5)
  for (j in 1:4){ 
    lines(f,mypw[[k]][,j])
    axis(4,mypw[[k]][hm,j],j*5)
  }
  par(cex.axis=1.0)
  title(paste(k+1,"groups",sep=" "))
}

# (d)
a=4
b=5
s=1
sa=1
sb=1
hm=1000
myf1=myf2=myf3=rep(0,hm)
for (j in 1:hm){
  mua=sa*rnorm(a)
  mub=sb*rnorm(b)
  dat=data.frame(y=apply(expand.grid(mub,mua),1,sum)+s*rnorm(20),
                 x1=factor(rep(LETTERS[1:a],rep(b,a))),
                 x2=factor(rep(LETTERS[1:b],a)))
  myf1[j]=summary(aov(y~x1+x2,dat))[[1]][1,4]
  myf2[j]=summary(aov(y~x1,dat))[[1]][1,4]
}
#
myf3=a*(b-1)/(a-1)*(s+b*sa)*rchisq(hm,a-1)/
     (s*rchisq(hm,(a-1)*(b-1))+(s+a*sb)*rchisq(hm,b-1))

# comparison of the two simulation methods
par(mfrow=c(2,1),las=1)
hist(myf2,xlim=c(0,20),breaks=51,prob=T)
hist(myf3,xlim=c(0,20),breaks=51,prob=T)
t.test(log(myf2),log(myf3))

# power: correct model
cv1=qf(0.95,a-1,(a-1)*(b-1))
cv2=cv1/(1+b*sa/s)
round(sum(myf1>cv1)/hm,3)           # simulation 
round(1-pf(cv2,a-1,(a-1)*(b-1)),3)  # theoretical value

# power: false model
cv=qf(0.95,a-1,a*(b-1))
round(sum(myf2>cv)/hm,3)            # simulation 


###############################################
# 3
###############################################

# data manipulation
dat=read.table(file="pittsburgh.txt",header=T)
n=nrow(dat)
x=matrix(ncol=7,nrow=n-6)
for (j in 0:6) x[,7-j]=dat[j+(1:(n-6)),2]
y=dat[-(1:6),1]
n=length(y)

# MLE estimates under constraints
t1=solve(t(x)%*%x)
t2=t1%*%(t(x)%*%y)
a=t2-apply(t1,1,sum)*sum(t2)/sum(t1)
b=apply(t1,1,sum)/sum(t1)

hm=101
gms=seq(-0.001,0.002,length=hm)
L=rep(0,hm)
for (j in 1:hm){
  gm=gms[j]
  bhat=a+gm*b
  rhat=y-x%*%bhat
  shat=mean(rhat^2)
  L[j]=exp(-n/2*(1+log(2*pi*shat)))
}

par(las=1)
plot(gms,L,type="l",xlab=expression(gamma),ylab="Likelihood")

##################
# End of hw7.R
##################
