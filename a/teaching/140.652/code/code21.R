## ------------------------------------------------------------------------
x <- rbind(c(18, 2), c(11, 9))
x
fisher.test(x) 

## ------------------------------------------------------------------------
x <- rbind(c(9,9), c(20, 62))
x
fisher.test(x) 

## ------------------------------------------------------------------------
# Simulate this
rhyper(1, 12, 18, 10) 

# Simulate this 100 times
rhyper(100, 12, 18, 10)

# Probability of getting no white balls in the sample
dhyper(0, 12, 18, 10)

# Probability of getting exactly 2 white balls in the sample
dhyper(2, 12, 18, 10)

# Probability of getting 2 or fewer white balls in the sample
phyper(2, 12, 18, 10)

# 2.5th and 97.5th %iles of number of white balls drawn
qhyper(c(0.025, 0.975), 12, 18, 10)

## ------------------------------------------------------------------------
# The data
x <- rbind(c(122, 117, 19, 244),
           c(1781, 1351, 288,3301),
           c(353, 269, 60, 713))
x

## ------------------------------------------------------------------------
fisher <- function(tab, n.sim=1000, return.all=FALSE, prnt=FALSE){
  bot0 <- sum(lgamma(tab+1))

  bot <- 1:n.sim
  a <- list(rep(row(tab),tab), rep(col(tab),tab))
  for(i in 1:n.sim) {
    a[[1]] <- sample(a[[1]])
    bot[i] <- sum(lgamma(table(a)+1))
    if(prnt) { if(i == round(i/10)*10) cat(i,"\n") }
  }
  if(return.all) return(list(bot0,bot))
  mean(bot0 <= bot)
}

# Apply Fisher's exact test, using 1000 simulated tables
set.seed(5)
fisher(x) 

# Built-in function (faster!)
fisher.test(x,simulate=TRUE,B=1000)

