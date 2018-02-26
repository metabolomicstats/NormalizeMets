#Correlation Test function in the DiffCorr package for R
#
# @n the number of samples
# @r the correlation coefficient
# @method "pearson" and "spearman" can be used.
cor2.test<-function (n, r, method = c("pearson", "kendall", "spearman")) 
{
  method <- match.arg(method)
  if (method == "pearson" || method == "spearman") {
    t <- abs(r) * sqrt((n - 2)/(1 - r^2))
    df <- n - 2
    p <- pt(t, df, lower.tail = FALSE) * 2
    c(p)
  }
  else {
    z <- abs(r)/sqrt((4 * n + 10)/(9 * n * (n - 1)))
    p <- pnorm(z, lower.tail = FALSE) * 2
    c(p)
  }
}
