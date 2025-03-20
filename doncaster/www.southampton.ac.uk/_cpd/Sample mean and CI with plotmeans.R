# Plot of a single sample mean with 95% CI.
##################################################################
library(gplots) # *** requires installation of package 'gplots' ***
mu <- 0         # the null hypothesis for the value of the population mean
alpha <- 0.05   # critical threshold for Type-I error
y <- c(-3.4076, -8.5488, -1.8395, -5.6626, -2.8127, -4.6421, 
       -1.6908, -4.7169, -0.6191, 1.2085, 5.1147, 7.2654, 4.4188, 
       6.3951, 0.5986, -4.5055, 13.322, 17.0426, 14.66, 16.1283, 
       8.5927, 2.2129, 1.8934, 2.0203) # the data
##################################################################
N <- length(y) ; DF <- N - 1
x <- c(rep("1", N))        # single sample x = 1
SE <- sqrt(var(y)/N)       # SE = sqrt(v/N)
CL <- SE*qt(1-alpha/2, DF) # confidence limits either side of mean
setcex <- 1.8              # set the font size for labels
par(cex = setcex, mar = c(4, 4, 1, 2) + 0.1) # define plot margins
# Plot the mean and CI using 'plotmeans' from the 'gplots' library
plotmeans(y ~ x, p = 1-alpha, las = 1, xlab = "", ylab = "", n.label = FALSE, xaxt = "n")
# Add  reference line, axis labels, and legend
windowsFonts(A = windowsFont("Times"),
             B = windowsFont("Arial"),
             C = windowsFont("Cambria"))
abline(h = mu, lty = 3) # reference line for null hypothesis
if (mu > (mean(y)+qt(1-alpha/2,DF))*0.9) {ps = 1} else {ps = 3} # text below/above line
text(1.3,mu,substitute(paste(italic("H"),""[0],":"~italic("\u03bc")~"= ",v), list(v=mu)),
     pos=ps, cex=setcex/2, family = "A")
mtext("x", font = 3, side = 1, line = 0.5, 
      las = 1, cex = setcex, family = "A")
mtext(substitute(paste("Response    mean and ",v,"% CI"), list(v=100*(1-alpha))), 
      side = 2, line = 2.5, las = 0, cex = setcex, family = "B")
mtext(expression("    "~italic("y")~"                  "), 
      side = 2, line = 2.5, las = 0, cex = setcex, family = "A")
mtext(expression(~
      bold("Fig. 1.")~" Sample from a population with unknown"~italic("\u03bc .")),
      outer = TRUE, side = 1, line = -1.5, cex = setcex, family = "C")
# Report statistics
writeLines(sprintf("Sample mean = %.2f g.",mean(y))) # mean (with trailing zero)
t <- (mean(y)-mu)/SE ; P <- 2*(1-pt(abs(t),DF)) # Student's t test
tcrit <- qt(1-alpha/2, DF) # critical t at alpha for a two-tailed test
if (abs(t) > tcrit) {result <- "Reject"} else {result <- "Fail to reject"}
writeLines(c(result, " H0: mu = ", mu, sprintf(" (t = %.3f, ",t),
             sprintf("DF = %.0f, ",DF), sprintf("P = %.3f; ",P),
             sprintf("t[alpha=%.2f] ",alpha), sprintf("= %.3f).",tcrit)
             ),sep = "")
#
rm(list = ls())
par(par(no.readonly = TRUE))