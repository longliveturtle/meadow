######################################################################
# Computer notes                                Biostatistics 140.751
# Homework 4                                 Johns Hopkins University 
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
# 3
###############################################

# the plot
par(mfrow=c(2,2),las=1,xaxs="i")
for(n in seq(5,20,5)){
q=3;x=seq(0,1,0.01);tst=qf(0.95,q,4*n-4)
lambda1=2*n*x^2*9/2;lambda2=2*n*x^2*9/4;lambda3=2*n*x^2*5/2
plot(c(0,1),c(0,1),type="n",xlab=expression(c/sigma),ylab="power",
     main=paste("n = ",as.character(n)))
lines(x,1-pf(tst,q,4*n-4,lambda1),col="red")
lines(x,1-pf(tst,q,4*n-4,lambda2),col="green")
lines(x,1-pf(tst,q,4*n-4,lambda3),col="blue")
text(0.55,0.2,"case 1",col="red")
text(0.7,0.2,"case 2",col="green")
text(0.85,0.2,"case 3",col="blue")
abline(h=0);abline(h=0.05,lty=2);axis(2,0.05,cex.axis=0.7)}

# a table
myp=as.list(1:4)
for(j in 1:4){
n=5*j;q=3;x=seq(0.1,0.9,0.1);tst=qf(0.95,q,4*n-4)
lambda1=2*n*x^2*9/2;lambda2=2*n*x^2*9/4;lambda3=2*n*x^2*5/2
pwmat=matrix(ncol=3,nrow=length(x))
pwmat[,1]=1-pf(tst,q,4*n-4,lambda1)
pwmat[,2]=1-pf(tst,q,4*n-4,lambda2)
pwmat[,3]=1-pf(tst,q,4*n-4,lambda3)
pwmat=round(pwmat,5)
pwmat=as.data.frame(pwmat)
names(pwmat)=c("case 1","case 2","case 3")
row.names(pwmat)=seq(0.1,0.9,0.1)
myp[[j]]=pwmat}
names(myp)=paste("n = ",seq(5,20,5))

# the simulation
n=10
hm=1000
cds=seq(0.1,0.7,0.1)

pval=matrix(0,ncol=length(cds),nrow=hm)
for (k in 1:length(cds)){
print(k)
for (j in 1:hm){
y=rnorm(4*n)+cds[k]*rep(c(0,0,3,3),n)
x=rep(LETTERS[1:4],n)
test=aov(y~as.factor(x))
pval[j,k]=anova(test)$"Pr(>F)"[1]}}
pval[pval<=0.05]=1
pval[pval<1]=0
phat1=apply(pval,2,sum)/hm

pst1=matrix(0,ncol=3,nrow=length(cds))
pst1[,2]=phat1
pst1[,1]=phat1-1.96*sqrt(phat1*(1-phat1)/hm)
pst1[,3]=phat1+1.96*sqrt(phat1*(1-phat1)/hm)
# note: I cheat a bit (but it doesn't make a great difference)
pst1[pst1>1]=1
pst1[pst1<0]=0
# end cheating

pval2=matrix(0,ncol=length(cds),nrow=hm)
for (k in 1:length(cds)){
print(k)
for (j in 1:hm){
y=rnorm(4*n)+cds[k]*rep(c(0,1.5,1.5,3),n)
x=rep(LETTERS[1:4],n)
test=aov(y~as.factor(x))
pval2[j,k]=anova(test)$"Pr(>F)"[1]}}
pval2[pval2<=0.05]=1
pval2[pval2<1]=0
phat2=apply(pval2,2,sum)/hm

pst2=matrix(0,ncol=3,nrow=length(cds))
pst2[,2]=phat2
pst2[,1]=phat2-1.96*sqrt(phat2*(1-phat2)/hm)
pst2[,3]=phat2+1.96*sqrt(phat2*(1-phat2)/hm)
pst2[pst2>1]=1
pst2[pst2<0]=0

par(mfrow=c(1,1),las=1,xaxs="i")
n=10;q=3;x=seq(0,1,0.01);tst=qf(0.95,q,4*n-4)
lambda1=2*n*x^2*9/2;lambda2=2*n*x^2*9/4;lambda3=2*n*x^2*5/2
plot(c(0,1),c(0,1),type="n",xlab=expression(c/sigma),ylab="power",
     main=paste("n = ",as.character(n)))
lines(x,1-pf(tst,q,4*n-4,lambda1),col="red",lty=2)
lines(x,1-pf(tst,q,4*n-4,lambda2),col="blue",lty=2)
text(0.1,0.4,"case 1",col="red",cex=2)
text(0.5,0.4,"case 2",col="blue",cex=2)
abline(h=0);abline(h=0.05,lty=2);axis(2,0.05,cex.axis=0.7)
for (j in 1:length(cds)){
lines(rep(cds[j],2),pst1[j,c(1,3)],lwd=1.5)
lines(cds[j]+c(-0.005,0.005),rep(pst1[j,2],2),lwd=1.5)
lines(cds[j]+c(-0.0025,0.0025),rep(pst1[j,1],2),lwd=1.5)
lines(cds[j]+c(-0.0025,0.0025),rep(pst1[j,3],2),lwd=1.5)}
for (j in 1:length(cds)){
lines(rep(cds[j],2),pst2[j,c(1,3)],lwd=1.5)
lines(cds[j]+c(-0.005,0.005),rep(pst2[j,2],2),lwd=1.5)
lines(cds[j]+c(-0.0025,0.0025),rep(pst2[j,1],2),lwd=1.5)
lines(cds[j]+c(-0.0025,0.0025),rep(pst2[j,3],2),lwd=1.5)}

##################
# End of hw4.R
##################
