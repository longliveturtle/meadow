######################################################################
# Computer notes                                 Biostatistics 140.752
# Homework 5                                  Johns Hopkins University 
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
# 1
###############################################

Y=c(75,67,70,75,65,71,67,67,76,68,
    57,58,60,59,62,60,60,57,59,61,
    58,61,56,58,57,56,61,60,57,58,
    58,59,58,61,57,56,58,57,57,59,
    62,66,65,63,64,62,65,65,62,67)

X=rep(c("C","G","F","G+F","S"),rep(10,5))

fc=unique(X);k=length(fc);my.fit=rep(0,k);n=rep(0,k)

for(j in 1:k){
my.fit[j]=mean(subset(Y,X==fc[j]))
n[j]=length(subset(Y,X==fc[j]))}
mf=rep(my.fit,n);s2=sum((Y-mf)^2)/(sum(n)-k)

cv=rep(0,4);level=0.95
cv[1]=qt(1-(1-level)/2,sum(n)-k)
cv[2]=qt(1-(1-level)/(2*k),sum(n)-k)
cv[3]=qt(0.5*(1+level^(1/k)),sum(n)-k)
cv[4]=sqrt(k*qf(level,k,sum(n)-k))

z=as.list(1:k);names(z)=unique(X)
for (j in 1:k){
z2=matrix(0,ncol=3,nrow=4)
for (i in 1:4){
z2[i,]=my.fit[j]+c(-1,0,1)*cv[i]*sqrt(s2/n[j])}
z2=as.data.frame(z2)
names(z2)=c("lower","fitted","upper")
row.names(z2)=c("unadjusted","bonferroni","max.modulus","scheffe")
z[[j]]=round(z2,rnd)}

print(z)


###############################################
# 2a
###############################################

Y=c(
541,524,561,414,410,457,344,467,464,498,580,471,525,508,566,635,603,714,865,
640,649,540,464,547,460,566,577,631,574,534,571,554,577,628,487,644,640,704,
648,968,587,699,632,591,782,510,610,524)

X=c(
52.5,57.2,58.0,52.9,54.4,57.1,45.1,55.3,52.9,55.2,53.0,52.5,57.4,54.5,60.8,
58.6,57.2,54.0,72.4,67.7,66.3,60.2,51.1,51.7,55.1,54.4,54.8,57.9,56.3,49.3,
51.8,51.3,57.8,54.7,48.7,62.9,56.6,58.6,66.3,67.2,62.6,56.3,60.3,50.8,67.2,
57.1,62.3,59.3)

# the mylm() function from HW#2
mylm=function(x,y,hatmat=F,addint=F){
if(addint) x=cbind(1,x)
if(det(t(x)%*%x)==0) stop("The design matrix is rank deficient!")
beta=as.vector(solve(t(x)%*%x)%*%t(x)%*%y)
rss=as.vector((t((y-x%*%beta))%*%(y-x%*%beta)))
sigma=sqrt(rss/(dim(x)[1]-dim(x)[2]))
varbeta=sigma^2*solve(t(x)%*%x)
hat=x%*%solve(t(x)%*%x)%*%t(x)
fitted=as.vector(x%*%beta)
residuals=as.vector(y-x%*%beta)
if(hatmat){
z=list(beta=beta,sigma=sigma,varbeta=varbeta,fitted=fitted,
       residuals=residuals,hat=hat)}
else{
z=list(beta=beta,sigma=sigma,varbeta=varbeta,fitted=fitted,
       residuals=residuals)}
return(z)}

# this function calculates the F statistics for various 
# b0/b1 combinations
my.image=function(b0,b1,X,Y,A,ai=T){
n0=length(b0);n1=length(b1);n=length(Y);d=dim(A)[1]
if(ai) X=cbind(rep(1,n),X)
mlm=mylm(X,Y)
bh=mlm$beta;s2=(mlm$sigma)^2
m2=solve(A%*%solve(t(X)%*%X)%*%t(A))
z=matrix(ncol=n0,nrow=n1)
for (j in 1:n0){
for (k in 1:n1){
b=c(b0[j],b1[k])
m1=A%*%(bh-b)
z[j,k]=drop(t(m1)%*%m2%*%m1/(d*s2))}}
return(z)}

# this function calculates unadjusted CI
ci.lm=function(X,Y,level=0.95,ai=T){
n=length(Y);lv=1-(1-level)/2
if(ai) X=cbind(rep(1,n),X)
if(is.vector(X)) k=1
if(is.matrix(X)) k=dim(X)[2]
mlm=mylm(X,Y)
z=matrix(nrow=k,ncol=3)
bt=qt(lv,n-k);be=sqrt(diag(mlm$varbeta))
z[,2]=mlm$beta
z[,1]=mlm$beta-bt*be
z[,3]=mlm$beta+bt*be
z=as.data.frame(z)
names(z)=c("lower","fitted","upper")
row.names(z)=0:(k-1)
return(z)}

# this function calculates the Scheffe CI
ci.sch=function(X,Y,h,level=0.95,ai=T){
n=length(Y)
if(ai) X=cbind(rep(1,n),X)
if(is.vector(X)) k=1
if(is.matrix(X)) k=dim(X)[2]
mlm=mylm(X,Y)
bh=mlm$beta;s2=(mlm$sigma)^2
mp=drop(h%*%bh);bm=solve(t(X)%*%X)
cl=sqrt(2*qf(level,2,n-2)*s2*drop(h%*%bm%*%h))
z=rep(0,3);z[1]=mp-cl;z[2]=mp;z[3]=mp+cl
return(z)}

# the actual F statistics for the joint confidence region in the example
A=diag(2);b0=seq(-600,100,length=151);b1=seq(8,20,length=151)
z=my.image(b0,b1,X,Y,A);z=log(z);z[z<0]=0

# calculate the unadjusted and simultaneous scheffe confidence intervals
ci1=ci.lm(X,Y);h=c(0,1);ci2a=ci.sch(X,Y,h);h=c(1,0);ci2b=ci.sch(X,Y,h)

# the figure with unadjusted CI
par(las=1,pty="s")
image(b0,b1,z,col=grey(seq(0.7,0.9,0.01)),
xlab=expression(beta[0]),ylab=expression(beta[1]))
contour(b0,b1,z,drawlabels=F,add=T,levels=c(log(qf(0.95,2,46))))
points(ci1[1,2],ci1[2,2],cex=1,pch=15)
text(-170,15,expression(paste("( ",hat(beta)[0]," , ",hat(beta)[1]," )")))
abline(v=ci1[1,c(1,3)],lty=2);abline(h=ci1[2,c(1,3)],lty=2)

# the figure with Scheffe CI
par(las=1,pty="s")
image(b0,b1,z,col=grey(seq(0.7,0.9,0.01)),
xlab=expression(beta[0]),ylab=expression(beta[1]))
contour(b0,b1,z,drawlabels=F,add=T,levels=c(log(qf(0.95,2,46))))
points(ci1[1,2],ci1[2,2],cex=1,pch=15)
text(-170,15,expression(paste("( ",hat(beta)[0]," , ",hat(beta)[1]," )")))
abline(h=ci2a[c(1,3)],lty=2);abline(v=ci2b[c(1,3)],lty=2)


###############################################
# 2b
###############################################

# this function calculates the values for the confidence band
my.band=function(x,X,Y,level=0.95){
n=length(Y);k=length(x)
x=cbind(rep(1,k),x)
X=cbind(rep(1,n),X)
mlm=mylm(X,Y)
bh=mlm$beta;s2=(mlm$sigma)^2
bm=solve(t(X)%*%X)
fc=2*qf(level,2,n-2)*s2
z=matrix(ncol=3,nrow=k)
for (j in 1:k){
xc=x[j,]
cl=sqrt(fc*drop(xc%*%bm%*%xc))
z[j,2]=drop(xc%*%bh)
z[j,1]=drop(xc%*%bh)-cl
z[j,3]=drop(xc%*%bh)+cl}
return(z)}

# calculate the values for the auto data
x=seq(40,75,1);rgs=my.band(x,X,Y)

# generate the plot
plot(X,Y,type="n")
polygon(c(x,rev(x)),c(rgs[,1],rev(rgs[,3])),col=grey(0.9))
points(X,Y,pch=15)
abline(lsfit(X,Y),lty=2)


###############################################
# 3
###############################################

n2=20;s2sq=1;n1=n2*2^(-1:3);s1sq=s2sq*2^(-3:3);hm=10000
myt=matrix(0,length(n1),length(s1sq))

for (i in 1:length(n1)){
print(i)
for (j in 1:length(s1sq)){
x1=matrix(sqrt(s1sq[j])*rnorm(hm*n1[i]),nrow=hm)
x2=matrix(rnorm(n2*hm),nrow=hm)
m1=apply(x1,1,mean)
m2=apply(x2,1,mean)
var1=apply(x1,1,var)
var2=apply(x2,1,var)
mys=((n1[i]-1)*var1+(n2-1)*var2)/(n1[i]+n2-2)
z=abs((m1-m2)/sqrt(mys*(1/n1[i]+1/n2)))
myt[i,j]=length(subset(z,z>qt(0.975,n1[i]+n2-2)))}}

print(round(myt/hm,2))


##################
# End of hw5.R
##################
