R commands and output:

## Load required libraries
require(car)
require(lme4)

## Read data
z = read.table(header=TRUE, as.is=TRUE, text="
batch value
B1 74 
B1 76 
B1 75 
B2 68 
B2 71 
B2 72
B3 75 
B3 77 
B3 77 
B4 72 
B4 74 
B4 73 
B5 79
B5 81
B5 79")

## Analysis of variance
z.aov = car::Anova(lm(value~batch, data=z))
print(z.aov)

##> Anova Table (Type II tests)
##>
##> Response: value
##>           Sum Sq Df F value    Pr(>F)    
##> batch     147.73  4  20.518 8.246e-05 ***
##> Residuals  18.00 10                      
##> ---
##> Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

## Calculation of the p-value of the test of batch effect
1 - pf(z.aov$"F value"[1], df1=z.aov$Df[1], df2=z.aov$Df[2])

##> [1] 8.24639e-05

## REML Fit and Estimates of Variance Components
z.lmer = lmer(value~1 + (1|batch), REML=TRUE, data=z)
summary(z.lmer)

##> Linear mixed model fit by REML ['lmerMod']
##> Formula: value ~ 1 + (1 | batch)
##>    Data: z
##> 
##> REML criterion at convergence: 62.8
##> 
##> Scaled residuals: 
##>      Min       1Q   Median       3Q      Max 
##> -1.90384 -0.53153  0.00484  0.61386  1.16817 
##> 
##> Random effects:
##>  Groups   Name        Variance Std.Dev.
##>  batch    (Intercept) 11.71    3.422   
##>  Residual              1.80    1.342   
##> Number of obs: 15, groups:  batch, 5
##> 
##> Fixed effects:
##>             Estimate Std. Error t value
##> (Intercept)   74.867      1.569   47.71

## Compute bootstrap confidence intervals
## (This may take a minute or two.)
confint.merMod(z.lmer, level=0.95, method="boot", nsim=5000,
               boot.type="perc", quiet=FALSE)

##> 21 message(s): boundary (singular) fit: see ?isSingular
##>
##>                  2.5 %    97.5 %
##> .sig01       0.9911935  5.782427
##> .sigma       0.7703256  1.912063
##> (Intercept) 71.7913451 78.000213

## The "boundary (singular) fit" message arises when an estmate
## of the variance hits 0. These messages can safely be ignored.
## Occasionally, the model will fail to converge for a bootstrap 
## sample and a warning will be printed. These messages can also 
## safely be ignored. 


