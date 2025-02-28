aickin.alpha <- function(n, d=diag(1, nrow=nrow(n), ncol=ncol(n)), epsilon=1e-7, level=0.95)
{
  # Version 1.0 (March 2013)
  #
	# This function computes the alpha coefficient described in:
  # Maximum Likelihood Estimation of Agreement in the Constant Predictive Probability Model, and Its Relation to Cohen-s Kappa
	# Mikel Aickin
	# Biometrics 46, 293-302, June 1990.
  #
  #
  # Function developed by 
  # Lawrence Joseph and Patrick Bélisle
  # Division of Clinical Epidemiology
  # Montreal General Hospital
  # Montreal, Qc, Can
  #
  # patrick.belisle@rimuhc.ca
  # http://www.medicine.mcgill.ca/epidemiology/Joseph/PBelisle/Aickin-Alpha-Agreement-R.html
  #
  # Please refer to our webpage for details on each argument.
  
	
	if (any(dim(n) != dim(d))) stop("n and d must be of equal dimensions.")
	if (diff(dim(n)) != 0) stop("n and d must be square matrices.")

  m <- nrow(n)
  n <- n + 1/(m^2)
  J <- rep(1, m)
  ssize <- sum(n)
  A <- sum(n*d)
  p0 <- A/ssize
  rows.tot <- as.vector(n%*%matrix(J, ncol=1))
  cols.tot <- as.vector(matrix(J, nrow=1)%*%n)
  pr <- rows.tot/ssize
  pc <- cols.tot/ssize
  s <- sum(matrix(pc, nrow=m, ncol=m, byrow=T)*pr*d)
  alpha <- (p0 - s) / (1 - s)
  continue <- T
  
  while (continue)
  {
    previous.alpha <- alpha
    pr.denominator <- ssize * (1 - alpha + alpha * as.vector(d%*%matrix(pc, ncol=1)) / s)
    pr <- rows.tot/pr.denominator
    pr[1] <- 1- sum(pr[-1])
    pc.denominator <- ssize * (1 - alpha + alpha * as.vector(matrix(pr, nrow=1)%*%d) / s)
    pc <- cols.tot/pc.denominator
    pc[1] <- 1- sum(pc[-1])
    s <- sum(matrix(pc, nrow=m, ncol=m, byrow=T)*pr*d)
    alpha <- (p0 - s) / (1 - s)

    continue <- abs(alpha-previous.alpha) > epsilon
  }
  
  prdiff <- pr[-1] - pr[1] # m-1 x 1
  pcdiff <- pc[-1] - pc[1] # m-1 x 1
  
  d2L.da2 <- - ssize * (1-s)/(1-alpha)/((1-alpha)*s+alpha)
  R <- alpha/s/((1-alpha)*s+alpha)
  U <- 1/alpha - 1
  d2L.dadpri <- - A * pcdiff * ((ssize/A)^2) # m-1 x 1
  d2L.dadpcj <- - A * prdiff * ((ssize/A)^2) # m-1 x 1
  
  d2L.dpridprj <- -sum(n[1,])/(pr[1]^2) + A*(2*s*U+1)*R*R*matrix(pcdiff, ncol=1) %*% matrix(pcdiff, nrow=1) # m-1 x m-1
  d2L.dpri2 <- -rows.tot[-1]/(pr[-1]^2) - rows.tot[1]/(pr[1]^2) + A*(2*s*U + 1)*R*R*(pcdiff^2)
  diag(d2L.dpridprj) <- d2L.dpri2
  
  d2L.dpridpcj <- -A*R + A*(2*s*U+1)*R*R* matrix(pcdiff, ncol=1) %*% matrix(prdiff, nrow=1) # m-1 x m-1
  d2L.dpridpci <- -2*A*R + A*(2*s*U+1)*R*R*prdiff*pcdiff
  diag(d2L.dpridpcj) <- d2L.dpridpci
  
  d2L.dpcidpcj <- -sum(n[,1])/(pc[1]^2) + A*(2*s*U+1)*R*R*matrix(prdiff, ncol=1) %*% matrix(prdiff, nrow=1) # m-1 x m-1
  d2L.dpci2 <- -cols.tot[-1]/(pc[-1]^2) - cols.tot[1]/(pc[1]^2) + A*(2*s*U + 1)*R*R*(prdiff^2)
  diag(d2L.dpcidpcj) <- d2L.dpci2

	#		The matrix of second derivatives is composed of vectors and matrices
	#		computed above. The variance of alpha-s estimator is the top-left component of its inverse.
	#		The matrix is symetric and constructed as follows:
  #
	#					 alpha				  Pr1  Pr2 ... Prq				    Pc1  Pc2 ... Pcq
  #
	#					 ------------		-----------------------    	-------------------
	#		alpha  - d2L.da2  -   -   t(d2L.dadpri)     -    	-  t(d2L.dadpcj)  -
	#					 ------------   -----------------------    	-------------------
  #
	#		Pr1                   ----------------------   	  ------------------
	#		Pr2                   -   d2L.dpridprj     -      -  d2L.dpridpcj  -
	#		...                   -                    -      -                -
	#		Prq                   ----------------------      ------------------
  #
	#		Pc1                                             	-------------------
	#		Pc2                                             	-  d2L.dpcidpcj   -
	#		...                                             	-                 -
	#		Pcq                                             	-------------------
  
  matrix.second.derivatives.top <- matrix(c(d2L.da2, d2L.dadpri, d2L.dadpcj), nrow=1)
  matrix.second.derivatives.middle <- cbind(matrix(d2L.dadpri, ncol=1), d2L.dpridprj, d2L.dpridpcj)
  matrix.second.derivatives.bottom <- cbind(matrix(d2L.dadpcj, ncol=1), t(d2L.dpridpcj), d2L.dpcidpcj)
  matrix.second.derivatives <- rbind(matrix.second.derivatives.top, matrix.second.derivatives.middle, matrix.second.derivatives.bottom)
  
  parms.cov.matrix <- solve(-matrix.second.derivatives)
  alpha.sd <- sqrt(parms.cov.matrix[1])
  z <- qnorm((1+level)/2)
  lcl <- alpha - z*alpha.sd
  ucl <- alpha + z*alpha.sd
	
	list(alpha=alpha, lcl=lcl, ucl=ucl)
}
