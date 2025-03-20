######################################################################
# Computer notes                                Biostatistics 140.751
# Homework 2                                 Johns Hopkins University 
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

mylm=function(x,y,hatmat=F,addint=T){
  if(addint) x=cbind(1,x)
  if(det(t(x)%*%x)==0) stop("The design matrix is rank deficient!")
  xtxinv=solve(t(x)%*%x)
  beta=as.vector(xtxinv%*%t(x)%*%y)
  fitted=x%*%beta
  residuals=as.vector(y-fitted)
  rss=sum(residuals^2)
  sigma=sqrt(rss/(nrow(x)-ncol(x)))
  varbeta=sigma^2*xtxinv
  hat=x%*%xtxinv%*%t(x)
  if(hatmat){
  z=list(beta=beta,sigma=sigma,varbeta=varbeta,fitted=fitted,
         residuals=residuals,hat=hat)}
  else{
  z=list(beta=beta,sigma=sigma,varbeta=varbeta,fitted=fitted,
         residuals=residuals)}
  return(z)}


##################
# End of hw2.R
##################
