######################################################################
##
## R commands and output for "Chi-Square Test Example"

## Attach the "nortest" library that contains the chi-square test.
library(nortest)
n=1000

## Generate normal random numbers and perform the chi-square test.
y1 = rnorm(n, mean = 0, sd = 1)
pearson.test(y1)

## >         Pearson chi-square normality test
## >
## > data:  y1 
## > P = 32.256, p-value = 0.3087

## Generate double exponential random numbers and perform
## the chi-square test.
y2 = ifelse(runif(n) > 0.5, 1, -1) * rexp(n) 
pearson.test(y2)

## >         Pearson chi-square normality test
## >
## > data:  y2 
## > P = 91.776, p-value = 1.935e-08

## Generate t random numbers and perform the chi-square test.
y3 = rt(n, 3)
pearson.test(y3)

## >         Pearson chi-square normality test
## >
## > data:  y3 
## > P = 101.488, p-value = 5.647e-10

## Generate lognormal random numbers and perform the chi-square test.
y4 = rlnorm(n, meanlog = 0, sdlog = 1)
z = pearson.test(y4)
z

## >         Pearson chi-square normality test
## >
## > data:  y4 
## > P = 1085.104, p-value < 2.2e-16

## Compute critical value.
qchisq(.05,z$n.classes-3,lower.tail=FALSE)

## > [1] 42.55697

######################################################################
##
## R commands and output for "Application Example"

## Input data from 1.4.2.9.1. "Background and Data"
## Measurements of fatigue life (thousands of cycles until rupture) of
## rectangular strips of 6061-T6 aluminum sheeting, subjected to
## periodic loading with maximum stress of 21,000 psi (pounds per
## square inch), as reported by Birnbaum and Saunders (1958)

x = c( 370, 1016, 1235, 1419, 1567, 1820,
       706, 1018, 1238, 1420, 1578, 1868,
       716, 1020, 1252, 1420, 1594, 1881,
       746, 1055, 1258, 1450, 1602, 1890,
       785, 1085, 1262, 1452, 1604, 1893,
       797, 1102, 1269, 1475, 1608, 1895,
       844, 1102, 1270, 1478, 1630, 1910,
       855, 1108, 1290, 1481, 1642, 1923,
       858, 1115, 1293, 1485, 1674, 1940,
       886, 1120, 1300, 1502, 1730, 1945,
       886, 1134, 1310, 1505, 1750, 2023,
       930, 1140, 1313, 1513, 1750, 2100,
       960, 1199, 1315, 1522, 1763, 2130,
       988, 1200, 1330, 1522, 1768, 2215,
       990, 1200, 1355, 1530, 1781, 2268,
      1000, 1203, 1390, 1540, 1782, 2440,
      1010, 1222, 1416, 1560, 1792)

## Maximum likelihood estimation of the parameters of the 3-parameter
## Weibull distribution based on the individual measured values

require(bbmle)
negloglik.wei = function (xi, beta, eta, x) {
    -sum(dweibull(x-xi, shape=beta, scale=eta, log=TRUE)) }
x.wei = mle2(negloglik.wei, data=list(x=x),
             optimizer="optim", method="L-BFGS-B",
             lower=c(xi=-Inf, beta=0, eta=0),
             start=list(xi=100, beta=2, eta=15))
xiHAT = coef(x.wei)[1]
betaHAT = coef(x.wei)[2]
etaHAT = coef(x.wei)[3]

## Definition of 10 bins so that each contains approximately 10
## observations, and counting the number of observations that fall
## into each bin

pVALUEs = seq(from=0.1, to=0.9, by=0.1)
qVALUEs = c(0, xiHAT + qweibull(pVALUEs, shape=betaHAT, scale=etaHAT), Inf)
counts = hist(x, breaks=qVALUEs, plot=FALSE)$counts

## Maximum likelihood estimation of the parameters of the 3-parameter
## Weibull distribution based on the bin counts (grouped data)

negloglik = function (xbe, counts, breaks)
{
    xi = xbe[1]; beta = xbe[2]; eta = xbe[3]
    nc = length(counts)
    s = 0
    for (ic in 1:nc)
    {
        s = s + counts[ic] *
            log(pweibull(breaks[ic+1]-xi, shape=beta, scale=eta) -
                pweibull(breaks[ic]-xi, shape=beta, scale=eta))
    }
    return(-1*s) 
}

wei.optim = optim(par=c(xi=xiHAT, beta=betaHAT, eta=etaHAT),
                  fn=negloglik, method="L-BFGS-B",
                  lower=c(xi=-Inf, beta=0, eta=0),
                  counts=counts, breaks=qVALUEs)
xiTILDE = wei.optim$par[1]
betaTILDE = wei.optim$par[2]
etaTILDE = wei.optim$par[3]

## Comparing the maximum likelihood estimates of the Weibull
## parameters derived from the individual observations and derived
## from the bin counts

cbind(Individual=c(xiHAT, betaHAT, etaHAT),
      Grouped=c(xiTILDE, betaTILDE, etaTILDE))
## >       Individual     Grouped
## > xi    181.156435  316.430889
## > beta    3.427996    3.096627
## > eta  1356.439166 1212.094501

## Computation of the expected numbers of observations in each bin
## according to the 3-parameter Weibull model fitted to the bin counts
## by the method of maximum likelihood

observed = counts
nc = length(counts)
expectedTILDE = numeric(nc)
breaks = qVALUEs
for (ic in 1:nc)
{
    expectedTILDE[ic] = nx *
        (pweibull(breaks[ic+1]-xiTILDE, shape=betaTILDE, scale=etaTILDE) -
         pweibull(breaks[ic]-xiTILDE, shape=betaTILDE, scale=etaTILDE))
}

cbind(observed, expectedTILDE)
## > observed expectedTILDE
## >        9      9.225443
## >       12     10.507247
## >        8     10.566971
## >       14     10.499642
## >        7     10.387510
## >       11     10.252879
## >       13     10.103506
## >        4      9.942504
## >       14      9.775073
## >        9      9.739226

## Test statistic and corresponding p-value

x2.TILDE = sum((observed-expectedTILDE)^2/expectedTILDE); x2.TILDE
## > 9.431575
1-pchisq(x2.TILDE, df=nc-1-3)
## > 0.1507216

######################################################################
