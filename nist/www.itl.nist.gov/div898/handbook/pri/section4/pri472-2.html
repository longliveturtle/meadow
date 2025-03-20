R commands and output:

## Read data and save relevant variables.
fname = "CATAPULT.DAT"
m = matrix(scan(fname,skip=50),ncol=13,byrow=T)
distance = m[,1]
order = m[,12]

## Save variables as factors.
height = as.factor(m[,2])
start = as.factor(m[,3])
bands = as.factor(m[,4])
length = as.factor(m[,5])
stop = as.factor(m[,6])

## Save numeric variables and interactions.
h = m[,7]
s = m[,8]
b = m[,9]
l = m[,10]
e = m[,11]
hs = h*s
hb = h*b
hl = h*l
he = h*e
sb = s*b
sl = s*l
se = s*e
bl = b*l
be = b*e
le = l*e

df = data.frame(distance,height,start,bands,length,stop,order,
     h,s,b,l,e,hs,hb,hl,he,sb,sl,se,bl,be,le)
df[order(df$order),1:7]

##>    distance height start bands length stop order
##> 1     28.00   3.25     0    1      0   80     1
##> 17    99.00      4    10    2      2   62     2
##> 16   126.50   4.75    20    2      4   80     3
##> 14   126.50   4.75     0    2      4   45     4
##> 15    45.00   3.25    20    2      4   45     5
##> 2     35.00   4.75     0    1      0   45     6
##> 18    45.00      4    10    1      2   62     7
##> 4     28.25   4.75    20    1      0   80     8
##> 10    85.00   4.75     0    1      4   80     9
##> 3      8.00   3.25    20    1      0   45    10
##> 12    36.50   4.75    20    1      4   45    11
##> 9     33.00   3.25     0    1      4   45    12
##> 19    84.50      4    10    2      2   62    13
##> 8     28.50   4.75    20    2      0   45    14
##> 5     33.50   3.25     0    2      0   45    15
##> 7     36.00   3.25    20    2      0   80    16
##> 6     84.00   4.75     0    2      0   80    17
##> 11    45.00   3.25    20    1      4   80    18
##> 20    37.50      4    10    1      2   62    19
##> 13   106.00   3.25     0    2      4   80    20


## Generate four plots.
center = which(height[]==4)
par(mfrow=c(2,2),bg=rgb(1,1,0.8))
qqnorm(distance)
qqline(distance, col = 2)
boxplot(distance, horizontal=TRUE, main="Box Plot", xlab="Distance")
hist(distance, main="Histogram", xlab="Distance")
plot(order, distance, xlab="Actual Run Order", ylab="Distance",
     main="Run Order Plot")
points(df$order[center],df$distance[center],pch=19)
par(mfrow=c(1,1))
  
# Plots of responses versus factor columns
par(mfrow=c(2,3),bg=rgb(1,1,0.8))
boxplot(distance~height, data=df, main="Distance by Band Height",
        xlab="Height",ylab="Distance")

boxplot(distance~start, data=df, main="Distance by Start Angle",
        xlab="Start Angle",ylab="Distance")

boxplot(distance~bands, data=df, main="Distance by Number of Bands",
        xlab="Number of Bands",ylab="Distance")

boxplot(distance~length, data=df, main="Distance by Arm Length",
        xlab="Arm Length",ylab="Distance")

boxplot(distance~stop, data=df, main="Distance by Stop Angle",
        xlab="Stop Angle",ylab="Distance")
par(mfrow=c(1,1))


## Fit a model with up to second order interactions.
q = lm(distance~h+s+b+l+e+hs+hb+hl+he+sb+sl+se+bl+be+le,data=df)
qq = summary(q)
qq

##> Call:
##> lm(formula = distance ~ h + s + b + l + e + hs + hb + hl + he + 
##>     sb + sl + se + bl + be + le, data = df)

##> Residuals:
##>       1       2       3       4       5       6       7       8       9      10 
##> -0.7813 -0.7812 -0.7812 -0.7812 -3.7000 -3.7000 -3.7000 -3.7000 -0.7812 -0.7812 
##>      11      12      13      14      15      16      17      18      19      20 
##> -0.7813 -0.7812 -3.7000 -3.7000 -3.7000 -3.7000 22.0500  6.8750  7.5500 -0.6250 

##> Coefficients:
##>             Estimate Std. Error t value Pr(>|t|)    
##> (Intercept)  57.5375     2.9691  19.378 4.18e-05 ***
##> h            13.4844     3.3196   4.062  0.01532 *  
##> s           -11.0781     3.3196  -3.337  0.02891 *  
##> b            19.4125     2.9691   6.538  0.00283 ** 
##> l            20.1406     3.3196   6.067  0.00373 ** 
##> e            12.0469     3.3196   3.629  0.02218 *  
##> hs           -2.7656     3.3196  -0.833  0.45163    
##> hb            4.6406     3.3196   1.398  0.23467    
##> hl            4.7031     3.3196   1.417  0.22950    
##> he            0.1094     3.3196   0.033  0.97529    
##> sb           -3.1719     3.3196  -0.955  0.39343    
##> sl           -1.1094     3.3196  -0.334  0.75502    
##> se            2.6719     3.3196   0.805  0.46601    
##> bl            7.6094     3.3196   2.292  0.08365 .  
##> be            2.8281     3.3196   0.852  0.44225    
##> le            3.1406     3.3196   0.946  0.39768    
##> ---
##> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

##> Residual standard error: 13.28 on 4 degrees of freedom
##> Multiple R-squared: 0.9709,     Adjusted R-squared: 0.8619 
##> F-statistic: 8.905 on 15 and 4 DF,  p-value: 0.02375 

## Generate normal probability plot of the effects.
## Save parameters in a vector, but remove intercept and the two largest
## parameters (b & l).
qef = q$coef
ex = c(1,4,5)
qef = qef[-ex]

## Sort effects and save labels.
sef = qef[order(qef)]
qlab = names(sef)

## Leave off the two largest effects, b & l.
#large = c(1,2)
#sef = sef[-large]
#qlab = qlab[-large]

## Generate theoretical quantiles.
ip = ppoints(length(sef))
zp = qnorm(ip)

## Generate normal probability plot of all effects (excluding the
## intercept).  Bands and length are not shown.
par(bg=rgb(1,1,0.8))
plot(zp,sef, 
     ylab="Parameter Estimate", xlab="Theoretical Quantiles",
     main="Normal Probability Plot of Saturated Model Effects")
qqline(sef, col=2)
abline(h=0, col=4)
text(-1.5,12,"b & l not shown",pos=4)

## Add labels for largest 4 effects (two are not shown).
small2 = c(1:(length(sef)-3))
text(zp[1],sef[1],label=qlab[1],pos=4,cex=1, font=2)
text(zp[-small2],sef[-small2],label=qlab[-small2],pos=2,cex=1, font=2)
par(mfrow=c(1,1))


## Generate effect estimates.
q = lm(distance~h+s+b+l+e+bl,data=df)
summary(q)

##> Call:
##> lm(formula = distance ~ h + s + b + l + e + bl, data = df)

##> Residuals:
##>     Min      1Q  Median      3Q     Max 
##> -23.091  -5.195  -1.528   6.993  22.050 

##> Coefficients:
##>             Estimate Std. Error t value Pr(>|t|)    
##> (Intercept)   57.537      2.847  20.212 3.33e-11 ***
##> h             13.484      3.183   4.237  0.00097 ***
##> s            -11.078      3.183  -3.481  0.00406 ** 
##> b             19.412      2.847   6.819 1.23e-05 ***
##> l             20.141      3.183   6.328 2.62e-05 ***
##> e             12.047      3.183   3.785  0.00227 ** 
##> bl             7.609      3.183   2.391  0.03264 *  
##> ---
##> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

##> Residual standard error: 12.73 on 13 degrees of freedom
##> Multiple R-squared: 0.9131,     Adjusted R-squared: 0.873 
##> F-statistic: 22.78 on 6 and 13 DF,  p-value: 3.452e-06

## Generate anova table.
anova(q)

##> Analysis of Variance Table

##> Response: distance
##>           Df Sum Sq Mean Sq F value    Pr(>F)    
##> h          1 2909.3  2909.3 17.9499 0.0009708 ***
##> s          1 1963.6  1963.6 12.1153 0.0040616 ** 
##> b          1 7536.9  7536.9 46.5023 1.226e-05 ***
##> l          1 6490.3  6490.3 40.0449 2.625e-05 ***
##> e          1 2322.0  2322.0 14.3268 0.0022708 ** 
##> bl         1  926.4   926.4  5.7161 0.0326399 *  
##> Residuals 13 2107.0   162.1                      
##> ---
##> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Perform lack-of-fit test.
lof = factor(paste(height,start,bands,length,stop,bands))
inner.model = lm(distance ~ h+s+b+l+e+bl, data = df)
outer.model = lm(distance ~ lof)
anova(inner.model, outer.model)

##> Analysis of Variance Table

##> Model 1: distance ~ h + s + b + l + e + bl
##> Model 2: distance ~ lof
##>   Res.Df     RSS Df Sum of Sq      F Pr(>F)
##> 1     13 2106.99                           
##> 2      2  133.25 11   1973.74 2.6931 0.3018


## Generate four plots.
center = which(height[]==4)
par(mfrow=c(2,2),bg=rgb(1,1,0.8))
qqnorm(q$residuals)
qqline(q$residuals, col = 2)
abline(h=0)
boxplot(q$residuals, horizontal=TRUE, main="Box Plot", xlab="Residual")
hist(q$residuals, main="Histogram", xlab="Residual")
plot(order, q$residuals, xlab="Actual Run Order", ylab="Residual",
     main="Run Order Plot")
points(df$order[center],q$residuals[center],pch=19)
par(mfrow=c(1,1))


## Plot residuals versus predicted response.
par(mfrow=c(1,1),bg=rgb(1,1,0.8))
plot(predict(q),q$residuals,ylab="Residual",
     xlab="Predicted Distance", col=4)
abline(h=0, col=2)
par(mfrow=c(1,1))


## Plots of residuals versus the factor variables
par(mfrow=c(2,3),bg=rgb(1,1,0.8))
plot(q$residuals~h, data=df, main="Residuals by Band Height",
        xlab="Height",ylab="Residual", xlim=c(-1.5,1.5), col=4)

plot(q$residuals~s, data=df, main="Residuals by Start Angle",
        xlab="Start Angle",ylab="Residual", xlim=c(-1.5,1.5), col=4)

plot(q$residuals~b, data=df, main="Residuals by Number of Bands",
        xlab="Number of Bands",ylab="Residual", xlim=c(-1.5,1.5), col=4)

plot(q$residuals~l, data=df, main="Residuals by Arm Length",
        xlab="Arm Length",ylab="Residual", xlim=c(-1.5,1.5), col=4)

plot(q$residuals~e, data=df, main="Residuals by Stop Angle",
        xlab="Stop Angle",ylab="Residual", xlim=c(-1.5,1.5), col=4)
par(mfrow=c(1,1))


## Add log(distance) to the data frame.
logdist = log(distance)
df = data.frame(df,logdist)

## Fit a model with up to second order interactions.
z = lm(logdist~h+s+b+l+e+hs+hb+hl+he+sb+sl+se+bl+be+le,data=df)
summary(z)

##> Call:
##> lm(formula = logdist ~ h + s + b + l + e + hs + hb + hl + he + 
##>     sb + sl + se + bl + be + le, data = df)

##> Residuals:
##>        1        2        3        4        5        6        7        8 
##> -0.05182 -0.05182 -0.05182 -0.05182 -0.07753 -0.07753 -0.07753 -0.07753 
##>        9       10       11       12       13       14       15       16 
##> -0.05182 -0.05182 -0.05182 -0.05182 -0.07753 -0.07753 -0.07753 -0.07753 
##>       17       18       19       20 
##>  0.38930  0.29844  0.23093  0.11612 

##> Coefficients:
##>             Estimate Std. Error t value Pr(>|t|)    
##> (Intercept)  3.85702    0.06865  56.186 6.01e-07 ***
##> h            0.25735    0.07675   3.353  0.02849 *  
##> s           -0.24174    0.07675  -3.150  0.03452 *  
##> b            0.34880    0.06865   5.081  0.00708 ** 
##> l            0.39437    0.07675   5.138  0.00680 ** 
##> e            0.26273    0.07675   3.423  0.02670 *  
##> hs          -0.02582    0.07675  -0.336  0.75348    
##> hb          -0.02035    0.07675  -0.265  0.80403    
##> hl          -0.01396    0.07675  -0.182  0.86457    
##> he          -0.04873    0.07675  -0.635  0.55999    
##> sb           0.00853    0.07675   0.111  0.91686    
##> sl           0.06775    0.07675   0.883  0.42724    
##> se           0.07955    0.07675   1.036  0.35855    
##> bl           0.01499    0.07675   0.195  0.85472    
##> be          -0.01152    0.07675  -0.150  0.88794    
##> le          -0.01120    0.07675  -0.146  0.89108    
##> ---
##> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

##> Residual standard error: 0.307 on 4 degrees of freedom
##> Multiple R-squared: 0.9564,     Adjusted R-squared: 0.7927 
##> F-statistic: 5.845 on 15 and 4 DF,  p-value: 0.0502 


## Fit model with significant effects.
z = lm(logdist~h+s+b+l+e,data=df)
summary(z)

##> Call:
##> lm(formula = logdist ~ h + s + b + l + e, data = df)

##> Residuals:
##>      Min       1Q   Median       3Q      Max 
##> -0.27259 -0.13143 -0.03024  0.12220  0.38930 

##> Coefficients:
##>             Estimate Std. Error t value Pr(>|t|)    
##> (Intercept)  3.85702    0.04702  82.035  < 2e-16 ***
##> h            0.25735    0.05257   4.896 0.000236 ***
##> s           -0.24174    0.05257  -4.599 0.000413 ***
##> b            0.34880    0.04702   7.419 3.26e-06 ***
##> l            0.39437    0.05257   7.502 2.87e-06 ***
##> e            0.26273    0.05257   4.998 0.000195 ***
##> ---
##> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1 

##> Residual standard error: 0.2103 on 14 degrees of freedom
##> Multiple R-squared: 0.9284,     Adjusted R-squared: 0.9028 
##> F-statistic: 36.28 on 5 and 14 DF,  p-value: 1.565e-07

## Generate anova table 
anova(z)

##> Analysis of Variance Table

##> Response: logdist
##>           Df  Sum Sq Mean Sq F value    Pr(>F)    
##> h          1 1.05968 1.05968  23.968 0.0002362 ***
##> s          1 0.93505 0.93505  21.149 0.0004134 ***
##> b          1 2.43324 2.43324  55.036 3.258e-06 ***
##> l          1 2.48839 2.48839  56.284 2.869e-06 ***
##> e          1 1.10443 1.10443  24.980 0.0001952 ***
##> Residuals 14 0.61896 0.04421                      
##> ---
##> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Perform lack-of-fit test.
lof = factor(paste(height,start,bands,length,stop,bands))
inner.model = lm(logdist ~ h+s+b+l+e, data = df)
outer.model = lm(logdist ~ lof)
anova(inner.model, outer.model)

##> Analysis of Variance Table

##> Model 1: logdist ~ h + s + b + l + e
##> Model 2: logdist ~ lof
##>   Res.Df     RSS Df Sum of Sq     F Pr(>F)
##> 1     14 0.61896                          
##> 2      2 0.02916 12   0.58980 3.371 0.2514


## Generate four plots.
center = which(height[]==4)
par(mfrow=c(2,2),bg=rgb(1,1,0.8))
qqnorm(z$residuals)
qqline(z$residuals, col = 2)
abline(h=0)
boxplot(z$residuals, horizontal=TRUE, main="Box Plot", xlab="Residual")
hist(z$residuals, main="Histogram", xlab="Residual")
plot(order, z$residuals, xlab="Actual Run Order", ylab="Residual",
     main="Run Order Plot")
points(df$order[center],z$residuals[center],pch=19)
par(mfrow=c(1,1))


## Plot residuals versus predicted response.
par(mfrow=c(1,1),bg=rgb(1,1,0.8))
plot(predict(z),z$residuals,ylab="Residual",
     xlab="Predicted log(Distance)", col=4)
abline(h=0,col=2)
par(mfrow=c(1,1))

## Plot of residuals versus the factor variables
par(mfrow=c(2,3),bg=rgb(1,1,0.8))
plot(z$residuals~h, data=df, main="Residuals by Band Height",
        xlab="Height",ylab="Residual", xlim=c(-1.5,1.5), col=4)

plot(z$residuals~s, data=df, main="Residuals by Start Angle",
        xlab="Start Angle",ylab="Residual", xlim=c(-1.5,1.5), col=4)

plot(z$residuals~b, data=df, main="Residuals by Number of Bands",
        xlab="Number of Bands",ylab="Residual", xlim=c(-1.5,1.5), col=4)

plot(z$residuals~l, data=df, main="Residuals by Arm Length",
        xlab="Arm Length",ylab="Residual", xlim=c(-1.5,1.5), col=4)

plot(z$residuals~e, data=df, main="Residuals by Stop Angle",
        xlab="Stop Angle",ylab="Residual", xlim=c(-1.5,1.5), col=4)
par(mfrow=c(1,1))


## Rearrange data so that factors and levels are in single columns.
group = rep(1:5,each=length(df$logdist))
newd = rep(logdist,5)
level = c(m[,7],m[,8],m[,9],m[,10],m[,11])
dflong = data.frame(group,level,newd)
dflong = dflong[order(group,level),]

## Generate means by factor and level.
gmn = aggregate(x=dflong$newd,by=list(dflong$group,dflong$level),FUN="mean")
cgroup = rep(c("Height","Start","Bands","Length","Stop"),3)
cgroup = cgroup[-8]
dfp = data.frame(cgroup,gmn)
names(dfp)=c("cgroup","group","level","tmean")

## Attach lattice library and generate main effects plot.
library(lattice)
xyplot(tmean~level|cgroup,data=dfp,layout=c(5,1),xlim=c(-2,2),
       ylab="log(Distance)",xlab="Factor Levels", type="b",
panel = function(x, y, ...){
panel.xyplot(x, y, ...)
panel.abline(h = mean(logdist), lty = 2, col = 2)})
par(mfrow=c(1,1))
