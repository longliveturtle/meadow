######################################################################
# Computer notes                                 Biostatistics 140.752
# Homework 8                                  Johns Hopkins University 
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

library(nlme)
library(lattice)
trellis.device(col=F)
data(Oxboys)

# a
plot(Oxboys)

# b
oxb.lm1=lm(height~age,Oxboys)
bwplot(Subject~resid(oxb.lm1),Oxboys)
summary(oxb.lm1)

# c
oxb.lmList=lmList(height~age,Oxboys)
plot(oxb.lmList,Subject~resid(.))
oxb.lmList

# d
plot(intervals(oxb.lmList))

# e
oxb.lme1=lme(height~age,Oxboys,random=~age|Subject)
plot(oxb.lme1,Subject~resid(.))

# f
plot(oxb.lme1)
qqnorm(oxb.lme1)

# g
plot(augPred(oxb.lme1))

# h
plot(oxb.lme1,resid(.)~age|Subject)

# i
oxb.lmList2=update(oxb.lmList,height~age+I(age^2))
plot(intervals(oxb.lmList2))

# j
oxb.lme2=lme(oxb.lmList2)
summary(oxb.lme2)
plot(augPred(oxb.lme2))

# k
intervals(oxb.lme2)
plot(oxb.lme2,resid(.)~age|Subject)



###############################################
# 1
###############################################

data(Alfalfa)
attach(Alfalfa)

# a
plot(Alfalfa)
interaction.plot(Date,Variety,Yield)

# b
alf.lme1=lme(Yield~Date*Variety,data=Alfalfa,random=~1|Block/Variety)
summary(alf.lme1)
intervals(alf.lme1)

# c
anova(alf.lme1)

# d
summary(aov(Yield~Variety*Date+Error(Block/Variety),data=Alfalfa))

# e
alf.lme2=update(alf.lme1,Yield~Date)
summary(alf.lme2)
intervals(alf.lme2)

# f
plot(alf.lme2)
qqnorm(alf.lme2)


##################
# End of hw8.R
##################
