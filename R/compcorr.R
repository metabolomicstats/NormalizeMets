#Compare two correlation coefficients using Fisher's Z-transformation: 
#A function in the DiffCorr package
#
#@n1 sample size under condition 1
#@r1 correlation coefficient under condition 1
#@n2 sample size under condition 2
#@r2 correlation coefficient under condition 1
#
compcorr<-function (n1, r1, n2, r2) 
{
  num1a <- which(r1 >= 0.99)
  num2a <- which(r2 >= 0.99)
  r1[num1a] <- 0.99
  r2[num2a] <- 0.99
  num1b <- which(r1 <= -0.99)
  num2b <- which(r2 <= -0.99)
  r1[num1b] <- -0.99
  r2[num2b] <- -0.99
  z1 <- atanh(r1)
  z2 <- atanh(r2)
  dz <- (z1 - z2)/sqrt(1/(n1 - 3) + (1/(n2 - 3)))
  pv <- 2 * (1 - pnorm(abs(dz)))
  return(list(diff = dz, pval = pv))
}
