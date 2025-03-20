# install.packages("mda")

library(mda)

# load and plot the tree data 
data(trees)
par(las=1,pty="s");pairs(trees,panel=panel.smooth,main="trees data")

# fit a MARS model
fit1 = mars(trees[,-3], trees[3])

# from the MARS helpfile
showcuts = function(obj)
{
  tmp = obj$cuts[obj$sel, ]
  dimnames(tmp) = list(NULL, names(trees)[-3])
  tmp
}

# output
# > showcuts(fit1)
#      Girth Height
# [1,]     0      0
# [2,]    12      0
# [3,]    12      0
# [4,]     0     76

# plot the fit (all of the below...)

x1=seq(8,21,length=21);x2=seq(60,90,length=21);x=expand.grid(x1,x2)
xprd=matrix(predict(fit1,x),ncol=21)

# this function projects 3D points onto a 2D plot
trans3d=function(x,y,z,pmat){
  tr=cbind(x,y,z,1) %*% pmat
  list(x = tr[,1]/tr[,4], y= tr[,2]/tr[,4])}

par(fig=c(0,0.65,0,1),pty="m")
persp(x1,x2,xprd,theta=30,phi=30,expand=0.8,col="lightblue",
xlab="girth",ylab="height",zlab="volume") -> theplot
points(trans3d(trees[,1],trees[,2],trees[,3],theplot),col=2,pch=16)
par(fig=c(0.65,1,0.5,1),new=T,las=1,cex.axis=0.6)
# this is from the MARS help file
Xp = matrix(sapply(trees[1:2], mean), nrow(trees), 2, byrow=TRUE)
xr = sapply(trees, range)
Xp1 = Xp; Xp1[,1] = seq(xr[1,1], xr[2,1], len=nrow(trees))
Xf = predict(fit1, Xp1)
plot(Xp1[ ,1], Xf, xlab=names(trees)[1], ylab="", type="l")
par(fig=c(0.65,1,0,0.5),new=T)
xr = sapply(trees, range)
Xp1 = Xp; Xp1[,2] = seq(xr[1,2], xr[2,2], len=nrow(trees))
Xf = predict(fit1, Xp1)
plot(Xp1[ ,2], Xf, xlab=names(trees)[2], ylab="", type="l")




