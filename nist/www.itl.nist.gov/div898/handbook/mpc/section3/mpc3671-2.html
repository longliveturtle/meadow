R commands and output:

## 2017/06: Correct line that reads data file to remove
##          the "skip=1" (this was actually skipping the
##          a line of good data).  This results in minor
##          changes to the output as well.

## Read data and save variables.
m = matrix(scan("LOADCELL.DAT", skip=25),ncol=2,byrow=T)
x = m[,1]
y = m[,2]
x2 = x*x

## Generate quadratic regression curve and print results.
z = lm(y ~ x + x2)
zz = summary(z)
zz

> Call:
> lm(formula = y ~ x + x2)

> Residuals:
>        Min         1Q     Median         3Q        Max 
> -9.969e-05 -1.467e-05  5.364e-06  2.030e-05  5.595e-05

> Coefficients:
>               Estimate Std. Error   t value Pr(>|t|)    
> (Intercept) -1.850e-05  2.081e-05    -0.889    0.378    
> x            1.001e-01  3.681e-06 27190.614   <2e-16 ***
> x2           7.030e-06  1.460e-07    48.137   <2e-16 ***
> ---
> Signif. codes:  0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1

> Residual standard error: 3.409e-05 on 51 degrees of freedom
> Multiple R-squared:      1,	Adjusted R-squared:      1 
> F-statistic: 7.639e+09 on 2 and 51 DF,  p-value: < 2.2e-16

## Print covariance matrix.
v = vcov(z)
v

>               (Intercept)             x            x2
> (Intercept)  4.330796e-10 -7.027501e-11  2.498281e-12
> x           -7.027501e-11  1.355352e-11 -5.244666e-13
> x2           2.498281e-12 -5.244666e-13  2.133052e-14

## Save coefficients and variances.
a = z$coef[1]
b = z$coef[2]
c = z$coef[3]
sa2 = v[1,1]
sb2 = v[2,2]
sc2 = v[3,3]
sy2 = zz$sigma**2

## Generate additional points on the x-axis for plotting.
xnew = seq(2,21,by=.1)
xnew2 = xnew^2
xp = data.frame(x=xnew,x2=xnew2)

## Predict response for given X' values.
yp = predict(z,new=xp)

## Generate calibration curve plot.
plot(x,y, main="Quadratic Calibration Curve for Load Cell 32066",
     xlab="Load, psi", ylab="Readings")
lines(xp$x,yp)
llab = paste("Y = ", round(a,7), " + ", round(b,7),
             "*X + ", round(c,7), "*X*X",sep="")
text(8,2,llab)
segments(x0=0,y0=1.344,x1=13.417,y1=1.344,lty=2,col="blue")
text(6,1.4, "Y'=1.344",col="blue")
segments(x0=13.417,y0=1.344,x1=13.417,y1=0,lty=2,col="blue")
text(15,0.6, "X'=13.417",col="blue")

## The equation for the quadratic calibration curve is:
## f = sqrt(-b + (b^2 - 4*c*(a-y)))/(2*c)
## The partial derivatives of f with respect to Y is:

dfdy = 1/sqrt(b^2 - 4*c*(a-y))

## The other partial derivatives are:

dfda = -1/sqrt(b^2 - 4*c*(a-y))

dfdb = (-1 + b/sqrt(b^2 - 4*c*(a-y)))/(2*c)

dfdc = (-4*a + 4*y)/(sqrt(b^2 - 4*c*(a-y))*(4*c)) - 
       (-b + sqrt(b^2 - 4*c*(a-y)))/(2*c^2)

## The standard deviation of X' is defined from propagation of error. 
u = sqrt(dfdy^2*sy2 + dfda^2*sa2 + dfdb^2*sb2 + dfdc^2*sc2)

## Plot uncertainty versus X'.
plot(y,u,type="n",xlab="Scale for Instrument Response",
     ylab="psi",
     main="Standard deviation of calibrated value X' for a given response Y'")
lines(spline(y,u))

## Print the covariance matrix for the loadcell data.
v

>               (Intercept)             x            x2
> (Intercept)  4.330796e-10 -7.027501e-11  2.498281e-12
> x           -7.027501e-11  1.355352e-11 -5.244666e-13
> x2           2.498281e-12 -5.244666e-13  2.133052e-14

## Save covariances.
sab = v[1,2]
sac = v[1,3]
sbc = v[2,3]

## Compute updated uncertainty.
unew = sqrt(u^2 + 2*dfda*dfdb*sab + 2*dfda*dfdc*sac + 2*dfdb*dfdc*sbc)

## Plot predicted value versus X'.
plot(y,unew,type="n",xlab="Scale for Instrument Response",
     ylab="psi",
     main="Standard deviation of calibrated value X' for a given response Y'")
lines(spline(y,unew))
