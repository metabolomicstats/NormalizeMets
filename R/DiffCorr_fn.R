#edited comp.2.cc.fdr function 
comp.2.cc.fdr_fn<-function (#output.file = "res.txt", 
                            data1, #metabolites by samples 
                            #must have rownames
                            data2, 
                            method = "pearson", 
          p.adjust.methods = "BH", #Not using the local method or the get.lfdr function
          #must be one of # c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
          #   "fdr", "none")
          threshold = 0.05 #using threshold=1 to include everything
          ) 
{
  #Compute metabolite-wise correlation
  cc1 <- cor(t(data1), method = method) 
  cc2 <- cor(t(data2), method = method)
  #Extract just the lower.tri
  ccc1 <- as.vector(cc1[lower.tri(cc1)])
  ccc2 <- as.vector(cc2[lower.tri(cc2)])
  n1 <- ncol(data1)
  n2 <- ncol(data2)
  n <- nrow(data1)
  N <- n * (n - 1)/2
  
  #######################
  p1 <- rep(1, N)
  p2 <- rep(1, N)
  pdiff <- rep(1, N)
  diff <- rep(1, N)
  #######################
  
  mol.names <- rownames(cc1)
  p1 <- cor2.test(n1, ccc1)
  p2 <- cor2.test(n2, ccc2)
  pdiff <- compcorr(n1, ccc1, n2, ccc2)$pval
  diff <- ccc1 - ccc2
  pdiff[(is.na(pdiff)) == TRUE] <- 1
#   if (p.adjust.methods == "local") {
#     p1.lfdr <- get.lfdr(p1)$lfdr
#     p2.lfdr <- get.lfdr(p2)$lfdr
#     pdiff.lfdr <- get.lfdr(pdiff)$lfdr
#   }
# #  else (p.adjust.methods == "BH" | p.adjust.methods == "bh") {
#   else {  
    p1.lfdr <- p.adjust(p1, method = p.adjust.methods)
    p2.lfdr <- p.adjust(p2, method = p.adjust.methods)
    pdiff.lfdr <- p.adjust(pdiff, method = p.adjust.methods)
#  }
  # else {
  #   p1.lfdr <- rep("not adjusted", N)
  #   p2.lfdr <- rep("not adjusted", N)
  #   pdiff.lfdr <- rep("not adjusted", N)
  # }
  myindex <- which((lower.tri(cc1)) == TRUE, arr.ind = TRUE)
  mol.names1 <- mol.names[myindex[, 2]]
  mol.names2 <- mol.names[myindex[, 1]]
  
  
  fin.ind <- pdiff.lfdr < threshold
  res <- data.frame(mol.names1[fin.ind], mol.names2[fin.ind], ccc1[fin.ind], 
               p1[fin.ind], ccc2[fin.ind], p2[fin.ind], pdiff[fin.ind], 
               diff[fin.ind], p1.lfdr[fin.ind], p2.lfdr[fin.ind], pdiff.lfdr[fin.ind])
  
  head <- c("metabolite X", "metabolite Y", "r1", "p1", "r2", "p2", 
            "p (difference)", "(r1-r2)", "p-adjusted (in cond. 1)", "p-adjusted (in cond. 2)", 
            "p-adjusted (difference)")
  colnames(res) <- head
  # write.table(res, file = output.file, row.names = FALSE, col.names = FALSE, 
  #             sep = "\t", quote = FALSE)
  return(res)
}