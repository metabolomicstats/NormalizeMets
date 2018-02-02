#loplot function in statTarget, edited as commented


loplot_fn<-function (x, z, i,dirout.loplot) 
{
  cn <- colnames(x)
  qcid <- grep("QC", cn)
  RSD30_CV = paste(rownames(x)[i], "_", i, ".pdf", sep = "")
  pdf(paste(dirout.loplot, RSD30_CV, sep = "/"), width = 6, 
      height = 6)
  graphics::layout(matrix(1:2, nrow = 2))
  numY <- 1:dim(x)[2]
  graphics::plot(numY, x[i, ], pch = 19, col = "yellow", ylab = c("Intensity"), 
                 xlab = c("Injection Order"), main = "Raw Peak")
  points(qcid, x[i, qcid], pch = 19, col = "blue")
  legend("top", c("Sample", "QC"), col = c("yellow", "blue"), 
         lty = 1, pch = 19, bty = "n", cex = 0.75, horiz = TRUE)
  loe <- loess(x[i, qcid] ~ qcid)
  points(numY, predict(loe, numY), type = "l", col = rgb(0, 
                                                         0, 0, 0.3), lwd = 4)
  graphics::plot(numY, z[i, ], pch = 19, col = "yellow", ylab = c("Intensity"), 
                 xlab = c("Injection Order"), main = "Corrected Peak")
  points(qcid, z[i, qcid], pch = 19, col = "blue")
  legend("top", c("Sample", "QC"), col = c("yellow", "blue"), 
         lty = 1, pch = 19, bty = "n", cex = 0.75, horiz = TRUE)
  dev.off()
}
